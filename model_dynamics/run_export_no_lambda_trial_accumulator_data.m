%% no-mutual-inhibition model: export MEG trial paths
% Two parameter modes are supported:
%   hierarchical85 - fit EB distributions to minb2_85 and extract rows 38:56
%   direct19       - use independently fitted minb2_MEG_19 point estimates
%
% Both modes use the current 9-slot vector
% [v0 b1_ratio b2_ratio k1 k2 sigma lambda tau t0], with parameter 7
% (lambda) equal to zero.

if ~exist('noLambdaOverride','var')
    noLambdaOverride = struct();
end
clearvars -except noLambdaOverride;
clc;

codeDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(codeDir);
addpath(codeDir);
parameterFile85 = fullfile(projectDir,'figure_data','model_results','model_fit_parameter_85_NOLambda.mat');
parameterFile19 = fullfile(projectDir,'figure_data','model_results','model_fit_parameter_MEG19_NOLambda.mat');
behaviorFile85 = fullfile(projectDir,'figure_data','behavioral_results','bhv_data_85.mat');
behaviorFile19 = fullfile(projectDir,'figure_data','behavioral_results','bhv_data_MEG_19.mat');
parameterMode = 'hierarchical85';

hierarchyConfig = struct();
hierarchyConfig.nPosteriorDraws = 200;
hierarchyConfig.seed = 20260721;
hierarchyConfig.boundaryClip = 1e-5;
hierarchyConfig.correlationRegularization = 0.25;
hierarchyConfig.minimumSEFraction = 0.02;
hierarchyConfig.megExpectedRows = 38:56;
hierarchyConfig.parameterVariant = 'minb2';
hierarchyConfig.parameterValueField = 'optParams_FullArray';
hierarchyConfig.fixedLambda = true;

bankConfig = struct();
bankConfig.nParameterDraws = 30;
bankConfig.nRepPerDraw = 100;
bankConfig.nGrid = 101;
bankConfig.maxSteps = 1000;
bankConfig.b0 = 300;
bankConfig.delta = 0.01;
bankConfig.seed = 20260809;
bankConfig.megSubjectRows = [];
bankConfig.conditionRows = [];
bankConfig.useParallel = true;
bankConfig.nWorkers = 4;

accumulatorConfig = struct();
accumulatorConfig.minimumPaths = 30;
accumulatorConfig.minimumBranchProbability = 0.01;
accumulatorConfig.minimumReliableTrialsPerSubject = 1;

if isfield(noLambdaOverride,'parameterMode')
    parameterMode = char(noLambdaOverride.parameterMode);
end
if isfield(noLambdaOverride,'parameterFile85')
    parameterFile85 = char(noLambdaOverride.parameterFile85);
end
if isfield(noLambdaOverride,'parameterFile19')
    parameterFile19 = char(noLambdaOverride.parameterFile19);
end
if isfield(noLambdaOverride,'behaviorFile85')
    behaviorFile85 = char(noLambdaOverride.behaviorFile85);
end
if isfield(noLambdaOverride,'behaviorFile19')
    behaviorFile19 = char(noLambdaOverride.behaviorFile19);
end
if isfield(noLambdaOverride,'hierarchy')
    hierarchyConfig = merge_struct_local( ...
        hierarchyConfig,noLambdaOverride.hierarchy);
end
if isfield(noLambdaOverride,'bank')
    bankConfig = merge_struct_local(bankConfig,noLambdaOverride.bank);
end
if isfield(noLambdaOverride,'accumulator')
    accumulatorConfig = merge_struct_local( ...
        accumulatorConfig,noLambdaOverride.accumulator);
end

assert(ismember(parameterMode,{'hierarchical85','direct19'}), ...
    'parameterMode must be hierarchical85 or direct19.');
assert(isfile(behaviorFile19),'Missing behavior file: %s',behaviorFile19);

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outputDir = fullfile(projectDir,'analysis_outputs','model_dynamics', ...
    ['no_lambda_trial_accumulator_' parameterMode '_' stamp]);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
hierarchyFile = fullfile(outputDir,'no_lambda_parameter_source.mat');
bankFile = fullfile(outputDir,'no_lambda_trial_trajectory_bank.mat');
accumulatorFile = fullfile(outputDir, ...
    'no_lambda_trial_accumulator_GLM_data.mat');

switch parameterMode
    case 'hierarchical85'
        assert(isfile(parameterFile85), ...
            'Missing parameter file: %s',parameterFile85);
        assert(isfile(behaviorFile85), ...
            'Missing behavior file: %s',behaviorFile85);
        fprintf('Fitting fixed-lambda=0 EB model to 85 subjects...\n');
        hierarchyArgs = namedargs2cell(hierarchyConfig);
        hierarchy = fit_hierarchical_eb_parameters_85( ...
            parameterFile85,behaviorFile85,behaviorFile19, ...
            hierarchyArgs{:});

    case 'direct19'
        assert(isfile(parameterFile19), ...
            'Missing parameter file: %s',parameterFile19);
        fprintf('Loading independently fitted MEG-19 no-lambda parameters...\n');
        hierarchy = make_direct_source_local( ...
            parameterFile19,behaviorFile19, ...
            bankConfig.nParameterDraws,38:56);
end

assert(all(abs(double(hierarchy.meg.parameters.draws(:,:,7,:))) ...
    < 1e-12,'all'),'No-lambda source contains a nonzero lambda draw.');
assert(isfield(hierarchy,'meta') && ...
    isfield(hierarchy.meta,'parameterNames') && ...
    numel(hierarchy.meta.parameterNames) == 9, ...
    'The no-lambda source is not a 9-parameter result.');
assert(size(hierarchy.meg.parameters.draws,3) == 9, ...
    'The no-lambda draw array does not have 9 parameter slots.');
hierarchy.meta.outputDirectory = outputDir;
hierarchy.meta.analysisRole = ...
    'refitted no-mutual-inhibition model; lambda fixed at zero';
save(hierarchyFile,'hierarchy','hierarchyConfig', ...
    'parameterMode','-v7.3');

fprintf('Generating stay/change trial trajectories with lambda=0...\n');
bankConfig.checkpointFile = bankFile;
bankArgs = namedargs2cell(bankConfig);
trialBankNoLambda = build_hierarchical_meg_trial_trajectory_bank( ...
    hierarchy,behaviorFile19,bankArgs{:});
trialBankNoLambda.meta.modelVariant = 'minb2-no-mutual-inhibition';
trialBankNoLambda.meta.lambda = 0;
trialBankNoLambda.meta.parameterMode = parameterMode;
trialBankNoLambda.meta.parameterSourceFile = hierarchyFile;
save(bankFile,'trialBankNoLambda','bankConfig','-v7.3');

fprintf('Building choice-aligned no-lambda accumulator data...\n');
accumulatorArgs = namedargs2cell(accumulatorConfig);
trialAccumulatorNoLambda = ...
    build_choice_aligned_trial_accumulator_data( ...
    trialBankNoLambda,accumulatorArgs{:});
trialAccumulatorNoLambda.meta.modelVariant = ...
    'minb2-no-mutual-inhibition';
trialAccumulatorNoLambda.meta.lambda = 0;
trialAccumulatorNoLambda.meta.parameterMode = parameterMode;
trialAccumulatorNoLambda.meta.trajectoryBankFile = bankFile;
save(accumulatorFile,'trialAccumulatorNoLambda', ...
    'accumulatorConfig','-v7.3');

fprintf('Saved no-lambda trajectory bank: %s\n',bankFile);
fprintf('Saved no-lambda trial GLM data: %s\n',accumulatorFile);

%% Local helpers
function hierarchy = make_direct_source_local( ...
        parameterFile,behaviorFile,nBlocks,rowsIn85)
    fitNames = {'ots_minb2','pls_minb2','otn_minb2','pln_minb2'};
    dataNames = {'ots_data','pls_data','otn_data','pln_data'};
    p = load(parameterFile,fitNames{:});
    d = load(behaviorFile,dataNames{:});
    nSubjects = numel(p.(fitNames{1}));
    nConditions = numel(fitNames);
    assert(nSubjects == 19, ...
        'The direct MEG parameter file must contain 19 subjects.');
    assert(numel(rowsIn85) == nSubjects);
    parameters = nan(nSubjects,nConditions,9);
    for c = 1:nConditions
        fits = p.(fitNames{c});
        assert(numel(fits) == nSubjects, ...
            'Condition %d does not contain 19 parameter rows.',c);
        for s = 1:nSubjects
            value = get_parameter_value_local(fits(s));
            assert(isnumeric(value) && numel(value) == 9);
            parameters(s,c,:) = reshape(double(value),1,1,9);
        end
    end
    assert(all(isfinite(parameters),'all'));
assert(all(abs(parameters(:,:,7)) < 1e-12,'all'), ...
    'Parameter 7 (lambda) must be exactly zero in every minb2 fit.');

    ids = extract_ids_local(d.ots_data);
    for c = 2:nConditions
        assert(isequal(extract_ids_local(d.(dataNames{c})),ids), ...
            'Subject order differs across MEG conditions.');
    end
    hierarchy.meta.created = datetime('now');
    hierarchy.meta.method = [ ...
        'independently fitted point estimates repeated in Monte Carlo ', ...
        'blocks; no hierarchical parameter uncertainty'];
    hierarchy.meta.parameterFile19 = parameterFile;
    hierarchy.meta.parameterNames = {'v0','b1_ratio','b2_ratio','k1','k2', ...
        'sigma','lambda','tau','t0'};
    hierarchy.meta.conditionNames = {'OT-social','PL-social', ...
        'OT-nonsocial','PL-nonsocial'};
    hierarchy.meta.fixedLambda = true;
    hierarchy.meg.rowsIn85 = rowsIn85(:)';
    hierarchy.meg.subjectIDs = ids;
    hierarchy.meg.parameters.mean = parameters;
    hierarchy.meg.parameters.draws = single(repmat( ...
        parameters,1,1,1,nBlocks));
end

function ids = extract_ids_local(data)
    ids = strings(numel(data),1);
    for s = 1:numel(data)
        raw = data(s).subsName;
        if iscell(raw) && isscalar(raw)
            raw = raw{1};
        end
        value = string(raw);
        ids(s) = value(1);
    end
end

function value = get_parameter_value_local(record)
    if isfield(record,'optParams_FullArray')
        value = record.optParams_FullArray;
    elseif isfield(record,'optParams')
        value = record.optParams;
    else
        error('No 9-slot parameter vector found in the minb2 record.');
    end
end

function destination = merge_struct_local(destination,source)
    names = fieldnames(source);
    for i = 1:numel(names)
        destination.(names{i}) = source.(names{i});
    end
end
