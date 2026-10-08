function run_slurm_recovery

root = fileparts(mfilename('fullpath'));
modelRoot = fileparts(root);
addpath(modelRoot);
cd(root);
codeRoot = fileparts(modelRoot);
inputFile = fullfile(codeRoot,'analysis_outputs','model_prediction_results','prediction_data_85.mat');
if ~isfile(inputFile)
    inputFile = fullfile(codeRoot,'figure_data','model_results','prediction_data_85.mat');
end
assert(isfile(inputFile),'Missing model input: %s',inputFile);

names = {'ots_preds','pls_preds','otn_preds','pln_preds'};

taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));
nRepeat = str2double(getenv('FIT_NREPEAT'));

assert(isfinite(taskID) && taskID >= 1, ...
    'Invalid SLURM_ARRAY_TASK_ID.');
assert(isfinite(nRepeat) && nRepeat >= 1 && nRepeat == floor(nRepeat), ...
    'Invalid FIT_NREPEAT.');
assert(taskID <= 4*nRepeat, ...
    'Array task ID exceeds 4 conditions x FIT_NREPEAT.');

conditionIndex = mod(taskID - 1, 4) + 1;
repeatIndex = floor((taskID - 1) / 4) + 1;
conditionName = names{conditionIndex};

ncpu = str2double(getenv('SLURM_CPUS_PER_TASK'));
assert(isfinite(ncpu) && ncpu >= 2, ...
    'At least 2 CPUs are required.');

tag = getenv('RUN_TAG');
assert(~isempty(regexp(tag, '^[A-Za-z0-9_-]+$', 'once')), ...
    'Invalid RUN_TAG.');

nSubject = str2double(getenv('FIT_NSUB'));
assert(isfinite(nSubject) && nSubject >= 1 && nSubject == floor(nSubject), ...
    'Invalid FIT_NSUB.');

assert(~isempty(which('particleswarm')) && ...
       ~isempty(which('patternsearch')), ...
       'Global Optimization Toolbox is unavailable.');

S = load(inputFile, conditionName);
data = S.(conditionName);

assert(isstruct(data) && isvector(data), ...
    'Condition data must be a struct vector.');

n = min(nSubject, numel(data));

outDir = fullfile(codeRoot,'analysis_outputs','model_recovery',tag, conditionName, ...
    sprintf('repeat_%03d', repeatIndex));

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

signature = [ ...
    sha256(inputFile) ':' ...
    sha256(fullfile(modelRoot,'model_recovery.m')) ':' ...
    sha256('run_slurm_recovery.m') ...
];

pending = [];

for s = 1:n
    resultFile = fullfile(outDir, sprintf('subject_%03d.mat', s));

    if exist(resultFile, 'file')
        old = load(resultFile, 'record');

        assert(strcmp(old.record.signature, signature), ...
            'Inputs or code changed. Use a new RUN_TAG.');

        assert(old.record.subjectIndex == s && ...
               old.record.conditionIndex == conditionIndex && ...
               old.record.repeatIndex == repeatIndex, ...
            'Checkpoint identity mismatch.');
    else
        pending(end+1) = s; %#ok<AGROW>
    end
end

fprintf('Condition %s, repeat %d: %d selected, %d pending.\n', ...
    conditionName, repeatIndex, n, numel(pending));

if ~isempty(pending)

    jobDir = fullfile(outDir, ...
        ['pool_' getenv('SLURM_JOB_ID')]);

    if ~exist(jobDir, 'dir')
        mkdir(jobDir);
    end

    cluster = parcluster('local');
    cluster.NumWorkers = min(ncpu - 1, numel(pending));
    cluster.NumThreads = 1;
    cluster.JobStorageLocation = jobDir;

    pool = parpool(cluster, cluster.NumWorkers);
    cleanupPool = onCleanup(@() delete(pool));

    ok = false(size(pending));

    parfor j = 1:numel(pending)

        s = pending(j);

        try
            seed = 20260916 + 1000000*repeatIndex + ...
                10000*conditionIndex + s;

            rng(seed, 'twister');

            t0 = tic;

            result = model_recovery(data(s));

            assert(isscalar(result), ...
                'Unexpected result size.');

            assert(numel(result.optParams) == 9 && ...
                   all(isfinite(result.optParams)) && ...
                   isfinite(result.nll) && ...
                   isfinite(result.bic), ...
                'Invalid recovery result.');

            record = struct( ...
                'result', result, ...
                'subjectIndex', s, ...
                'conditionIndex', conditionIndex, ...
                'condition', conditionName, ...
                'repeatIndex', repeatIndex, ...
                'seed', seed, ...
                'signature', signature, ...
                'elapsedSeconds', toc(t0), ...
                'matlabVersion', version, ...
                'parameterNames', {{ ...
                    'v0','b1_ratio','b2_ratio','k1','k2', ...
                    'sigma','lambda','tau','t0'}});

            saveRecord( ...
                fullfile(outDir, sprintf('subject_%03d.mat', s)), ...
                record);

            fprintf('SAVED %s repeat %d subject %d\n', ...
                conditionName, repeatIndex, s);

            ok(j) = true;

        catch ME
            writeError( ...
                fullfile(outDir, sprintf('error_%03d.txt', s)), ME);

            fprintf(2, ...
                'FAILED %s repeat %d subject %d: %s\n', ...
                conditionName, repeatIndex, s, ME.message);
        end
    end

    clear cleanupPool

    assert(all(ok), ...
        'Some fits failed. Inspect error files and resubmit the same task.');
end

results = repmat(struct( ...
    'optParams', [], ...
    'nll', [], ...
    'bic', []), n, 1);

for s = 1:n
    z = load(fullfile(outDir, sprintf('subject_%03d.mat', s)), ...
        'record');
    results(s,1) = z.record.result;
end

subjectIndex = (1:n)';
parameterNames = z.record.parameterNames;
condition = conditionName;

save(fullfile(outDir, 'results_all.mat'), ...
    'results', 'subjectIndex', 'parameterNames', ...
    'condition', 'conditionIndex', 'repeatIndex', ...
    'signature', '-v7.3');

fprintf('COMPLETE %s repeat %d, %d subjects.\n', ...
    conditionName, repeatIndex, n);

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
