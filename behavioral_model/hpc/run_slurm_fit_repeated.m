function run_slurm_fit_repeated
% One condition x repeat per array task; subjects run on local workers.
root = fileparts(mfilename('fullpath'));
modelRoot = fileparts(root);
addpath(modelRoot);
cd(root);
codeRoot = fileparts(modelRoot);
inputFile = fullfile(codeRoot,'figure_data','behavioral_results','bhv_data_85.mat');
assert(isfile(inputFile),'Missing model input: %s',inputFile);

names = {'ots_data','pls_data','otn_data','pln_data'};
taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
nRepeats = str2double(getenv('FIT_NREPEAT'));
assert(isfinite(nRepeats) && nRepeats >= 10 && nRepeats == floor(nRepeats), ...
    'FIT_NREPEAT must be an integer >= 10.');
assert(isfinite(taskID) && taskID >= 1 && taskID <= 4*nRepeats && ...
    taskID == floor(taskID), 'Invalid array task index.');
ci = mod(taskID-1,4)+1;
repeatIndex = floor((taskID-1)/4)+1;
ncpu = str2double(getenv('SLURM_CPUS_PER_TASK'));
assert(isfinite(ncpu) && ncpu >= 2, 'At least 2 CPUs required.');
tag = getenv('RUN_TAG');
assert(~isempty(regexp(tag,'^[A-Za-z0-9_-]+$','once')), 'Invalid RUN_TAG.');
limit = str2double(getenv('FIT_NSUB'));
assert(isfinite(limit) && limit >= 1 && limit == floor(limit), 'Invalid FIT_NSUB.');
assert(~isempty(which('particleswarm')) && ~isempty(which('patternsearch')), ...
    'Global Optimization Toolbox functions not found.');
S = load(inputFile,names{ci});
data = S.(names{ci});
assert(isstruct(data) && isvector(data), 'Data must be a struct vector.');
n = min(limit,numel(data));
out = fullfile(codeRoot,'analysis_outputs','model_fitting','full',tag,names{ci},sprintf('repeat_%03d',repeatIndex));
if ~exist(out,'dir'), mkdir(out); end
% Exact content hashes prevent resuming after input/code changes.
signature = [sha256(inputFile) ':' sha256(fullfile(modelRoot,'main_model.m')) ':' ...
    sha256('run_slurm_fit_repeated.m')];
pending = [];
for s = 1:n
    f = fullfile(out,sprintf('subject_%03d.mat',s));
    if exist(f,'file')
        old = load(f,'record');
        assert(strcmp(old.record.signature,signature), ...
            'Inputs changed. Use a new RUN_TAG.');
        assert(old.record.subjectIndex == s && old.record.conditionIndex == ci && ...
            old.record.repeatIndex == repeatIndex, ...
            'Checkpoint identity mismatch.');
    else
        pending(end+1) = s; %#ok<AGROW>
    end
end
fprintf('Condition %s: %d selected, %d pending.\n',names{ci},n,numel(pending));
if ~isempty(pending)
    jobDir = fullfile(out,['pool_' getenv('SLURM_JOB_ID')]);
    if ~exist(jobDir,'dir'), mkdir(jobDir); end
    c = parcluster('local');
    c.NumWorkers = min(ncpu-1,numel(pending));
    c.NumThreads = 1;
    c.JobStorageLocation = jobDir;
    pool = parpool(c,c.NumWorkers);
    cleanup = onCleanup(@() delete(pool));
    ok = false(size(pending));
    parfor j = 1:numel(pending)
        s = pending(j);
        try
            assert(s <= 100000, 'Subject index exceeds seed allocation.');
            seed = 20260916 + ((repeatIndex-1)*4+ci-1)*100000 + s;
            assert(seed <= 2^32-1, 'Repeat count exceeds seed allocation.');
            rng(seed,'twister');
            ticFit = tic;
            result = main_model(data(s));
            assert(numel(result)==1 && numel(result.optParams)==9 && ...
                all(isfinite(result.optParams)) && isfinite(result.nll) && ...
                isfinite(result.bic), 'Invalid fit output.');
            record = struct('result',result,'subjectIndex',s, ...
                'conditionIndex',ci,'condition',names{ci},'repeatIndex',repeatIndex,'seed',seed, ...
                'signature',signature,'elapsedSeconds',toc(ticFit), ...
                'matlabVersion',version,'parameterNames', ...
                {{'v0','b1_ratio','b2_ratio','k1','k2','sigma','lambda','tau','t0'}});
            saveRecord(fullfile(out,sprintf('subject_%03d.mat',s)),record);
            fprintf('SAVED %s subject %d\n',names{ci},s);
            ok(j) = true;
        catch ME
            writeError(fullfile(out,sprintf('error_%03d.txt',s)),ME);
            fprintf(2,'FAILED %s subject %d: %s\n',names{ci},s,ME.message);
        end
    end
    clear cleanup
    assert(all(ok), 'Some fits failed; inspect error files and resume.');
end
for s = 1:n
    z = load(fullfile(out,sprintf('subject_%03d.mat',s)),'record');
    results(s,1) = z.record.result; %#ok<AGROW>
end
subjectIndex = (1:n)';
parameterNames = z.record.parameterNames;
condition = names{ci};
save(fullfile(out,'results_all.mat'),'results','subjectIndex', ...
    'parameterNames','condition','signature','repeatIndex','-v7.3');
fprintf('COMPLETE %s, %d subjects.\n',condition,n);
end

function saveRecord(path,record)
tmp = [path '.tmp.mat'];
save(tmp,'record','-v7.3');
[ok,msg] = movefile(tmp,path);
assert(ok,msg);
end

function writeError(path,ME)
fid = fopen(path,'w');
if fid >= 0
    cl = onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',getReport(ME,'extended','hyperlinks','off'));
end
end

function hex = sha256(path)
fid = fopen(path,'rb');
assert(fid >= 0, 'Cannot open file for hashing.');
cl = onCleanup(@() fclose(fid));
md = java.security.MessageDigest.getInstance('SHA-256');
while true
    b = fread(fid,1048576,'*uint8');
    if isempty(b), break; end
    md.update(typecast(b,'int8'));
end
d = typecast(md.digest(),'uint8');
hex = lower(reshape(dec2hex(d,2).',1,[]));
end
