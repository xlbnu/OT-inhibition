function run_slurm_compare16
% One condition x one model x one repeat per array task.
% Subjects run on local workers.
root = fileparts(mfilename('fullpath'));
modelRoot = fileparts(root);
addpath(modelRoot);
cd(root);
inputFile = fullfile(root,'bhv_data_85.mat');
if ~isfile(inputFile)
    inputFile = fullfile(fileparts(modelRoot),'figure_data','behavioral_results','bhv_data_85.mat');
end
assert(isfile(inputFile),'Provide the required input MAT in hpc or figure_data.');
names = {'ots_data','pls_data','otn_data','pln_data'};
modelNames = {'1_FullModel','2_No_Lambda', ...
              '3_No_t0_tau','4_No_Lambda_t0_tau'};

nCondition = numel(names);
nModel = numel(modelNames);
tasksPerRepeat = nCondition * nModel;

taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
nRepeat = str2double(getenv('FIT_NREPEAT'));

assert(isfinite(taskID) && taskID >= 1 && taskID == floor(taskID), ...
    'Invalid SLURM_ARRAY_TASK_ID.');

assert(isfinite(nRepeat) && nRepeat >= 1 && ...
       nRepeat == floor(nRepeat), ...
    'Invalid FIT_NREPEAT.');

assert(taskID <= tasksPerRepeat * nRepeat, ...
    'Array task ID exceeds condition x model x repeat count.');

slot = mod(taskID - 1, tasksPerRepeat);

ci = mod(slot, nCondition) + 1;
modelType = floor(slot / nCondition) + 1;
repeatIndex = floor((taskID - 1) / tasksPerRepeat) + 1;

conditionName = names{ci};
modelName = modelNames{modelType};

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
assert(numel(data)>=limit, 'Fewer subjects than requested.');
n = limit;
out = fullfile(root, 'outputs', tag, conditionName, modelName, sprintf('repeat_%03d', repeatIndex));
if ~exist(out,'dir'), mkdir(out); end
% Exact content hashes prevent resuming after input/code changes.
signature = [sha256(inputFile) ':' sha256(fullfile(modelRoot,'main_model_compare_one.m')) ':' ...
    sha256('run_slurm_compare16.m')];
pending = [];
for s = 1:n
    f = fullfile(out,sprintf('subject_%03d.mat',s));
    if exist(f,'file')
        old = load(f,'record');
        assert(strcmp(old.record.signature,signature), ...
            'Inputs changed. Use a new RUN_TAG.');
        assert(old.record.subjectIndex == s && old.record.conditionIndex == ci && ...
            old.record.repeatIndex == repeatIndex && old.record.modelType == modelType, ...
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
            seed = 20260916 + ...
                1000000 * repeatIndex + ...
                10000 * modelType + ...
                100 * ci + s;
            assert(seed <= 2^32-1, 'Repeat count exceeds seed allocation.');
            rng(seed,'twister');
            ticFit = tic;
            wrapped = main_model_compare_one(data(s),modelType);
            result = wrapped.modelFit;
            assert(isscalar(result) && numel(result.optParams_FullArray)==9 && ...
                all(isfinite(result.optParams_FullArray)) && isfinite(result.nll) && ...
                isfinite(result.bic), 'Invalid fit output.');
            record = struct('result',result,'subjectIndex',s, ...
                'conditionIndex',ci,'condition',names{ci},'repeatIndex',repeatIndex,'seed',seed, ...
                'modelType',modelType,'modelName',modelName, ...
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
results = repmat(struct( ...
    'optParams_FullArray', [], ...
    'nll', [], ...
    'bic', [], ...
    'modelName', [], ...
    'nParams', [], ...
    'nTrials', [], ...
    'modelType', [], ...
    'b0', []), n, 1);
for s = 1:n
    z = load(fullfile(out,sprintf('subject_%03d.mat',s)),'record');
    results(s,1) = z.record.result; 
end
subjectIndex = (1:n)';
parameterNames = z.record.parameterNames;
condition = names{ci};
save(fullfile(out,'results_all.mat'),'results','subjectIndex', ...
    'parameterNames','condition','modelType','modelName','signature','repeatIndex','-v7.3');
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
