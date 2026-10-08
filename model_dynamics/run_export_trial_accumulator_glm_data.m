%% Export subject-by-condition-by-trial accumulator trajectories for GLM
% Required upstream result: hierarchy from run_hierarchical_eb_85_to_meg19.
% Full defaults generate 30 parameter draws x 100 paths = 3000 paths/trial.

if ~exist('trialAccumulatorOverride','var')
    trialAccumulatorOverride = struct();
end
clearvars -except trialAccumulatorOverride;
clc;

codeDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(codeDir);
addpath(codeDir);
behaviorFile19 = fullfile(projectDir,'figure_data','behavioral_results','bhv_data_MEG_19.mat');

hierarchyFile = find_latest_hierarchy_local(projectDir);
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

if isfield(trialAccumulatorOverride,'hierarchyFile')
    hierarchyFile = char(trialAccumulatorOverride.hierarchyFile);
end
if isfield(trialAccumulatorOverride,'behaviorFile19')
    behaviorFile19 = char(trialAccumulatorOverride.behaviorFile19);
end
if isfield(trialAccumulatorOverride,'bank')
    bankConfig = merge_struct_local( ...
        bankConfig,trialAccumulatorOverride.bank);
end
if isfield(trialAccumulatorOverride,'accumulator')
    accumulatorConfig = merge_struct_local( ...
        accumulatorConfig,trialAccumulatorOverride.accumulator);
end

assert(isfile(hierarchyFile),'Hierarchy file not found: %s',hierarchyFile);
assert(isfile(behaviorFile19),'Behavior file not found: %s',behaviorFile19);
loaded = load(hierarchyFile,'hierarchy');
assert(isfield(loaded,'hierarchy'), ...
    'The hierarchy file does not contain variable hierarchy.');
hierarchy = loaded.hierarchy;
assert(isfield(hierarchy,'meta') && ...
    isfield(hierarchy.meta,'parameterNames') && ...
    numel(hierarchy.meta.parameterNames) == 9, ...
    'The hierarchy file is not a 9-parameter result.');
assert(~isfield(hierarchy.meta,'fixedLambda') || ...
    ~hierarchy.meta.fixedLambda, ...
    'The full-model export received a fixed-lambda hierarchy.');
assert(size(hierarchy.meg.parameters.draws,3) == 9, ...
    'The hierarchy draw array does not have 9 parameter slots.');

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outputDir = fullfile(projectDir,'analysis_outputs','model_dynamics', ...
    ['trial_accumulator_GLM_' stamp]);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
bankFile = fullfile(outputDir,'MEG_trial_trajectory_bank.mat');
accumulatorFile = fullfile(outputDir, ...
    'MEG_trial_accumulator_GLM_data.mat');
bankConfig.checkpointFile = bankFile;

fprintf('Hierarchy: %s\n',hierarchyFile);
fprintf('Generating branch-conditioned trial trajectories...\n');
bankArgs = namedargs2cell(bankConfig);
trialBank = build_hierarchical_meg_trial_trajectory_bank( ...
    hierarchy,behaviorFile19,bankArgs{:});
trialBank.meta.hierarchyFile = hierarchyFile;
trialBank.meta.outputDirectory = outputDir;
save(bankFile,'trialBank','bankConfig','-v7.3');
fprintf('Saved trajectory bank: %s\n',bankFile);

fprintf('Building observed-choice-aligned accumulator data...\n');
accumulatorArgs = namedargs2cell(accumulatorConfig);
trialAccumulator = build_choice_aligned_trial_accumulator_data( ...
    trialBank,accumulatorArgs{:});
trialAccumulator.meta.hierarchyFile = hierarchyFile;
trialAccumulator.meta.trajectoryBankFile = bankFile;
trialAccumulator.meta.outputDirectory = outputDir;
save(accumulatorFile,'trialAccumulator','accumulatorConfig','-v7.3');
fprintf('Saved GLM data: %s\n',accumulatorFile);

%% Example: primary endpoint-aligned dynamics outcome
% s = 1; c = 1;
% Y = double(trialAccumulator.subject{s,c}.end.choiceAligned. ...
%     baselineCorrected.B0.difference);
% reliable = trialAccumulator.subject{s,c}.end.choiceAligned.isReliable;
% Y(~reliable) = NaN;

%% Local helpers
function hierarchyFile = find_latest_hierarchy_local(projectDir)
    candidates = dir(fullfile(projectDir,'analysis_outputs','model_dynamics', ...
        'hierarchical_eb_85_to_meg19_*', ...
        'hierarchical_parameters_85.mat'));
    assert(~isempty(candidates),[ ...
        'No hierarchical_parameters_85.mat was found under analysis_outputs/model_dynamics. ', ...
        'Run run_hierarchical_eb_85_to_meg19 first or set ', ...
        'trialAccumulatorOverride.hierarchyFile.']);
    [~,order] = sort([candidates.datenum],'descend');
    latest = candidates(order(1));
    hierarchyFile = fullfile(latest.folder,latest.name);
end

function destination = merge_struct_local(destination,source)
    names = fieldnames(source);
    for i = 1:numel(names)
        destination.(names{i}) = source.(names{i});
    end
end
