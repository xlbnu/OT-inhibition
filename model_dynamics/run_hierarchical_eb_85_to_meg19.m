%% joint hierarchical EB fit in 85 subjects and MEG-19 trajectory export
% The exact trial likelihood in the original fitting code is stochastic.
% This runnable workflow therefore uses a transparent two-stage
% empirical-Bayes approximation and preserves the unpooled estimates for
% sensitivity analysis.

if ~exist('hierarchicalTrajectoryOverride','var')
    hierarchicalTrajectoryOverride = struct();
end
clearvars -except hierarchicalTrajectoryOverride;
clc;

codeDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(codeDir);
addpath(codeDir);
parameterFile85 = fullfile(projectDir,'figure_data','model_results','model_fit_parameter_85.mat');
behaviorFile85 = fullfile(projectDir,'figure_data','behavioral_results','bhv_data_85.mat');
behaviorFile19 = fullfile(projectDir,'figure_data','behavioral_results','bhv_data_MEG_19.mat');

hierarchyConfig = struct();
hierarchyConfig.nPosteriorDraws = 200;
hierarchyConfig.seed = 20260721;
hierarchyConfig.boundaryClip = 1e-5;
hierarchyConfig.correlationRegularization = 0.25;
hierarchyConfig.minimumSEFraction = 0.02;
hierarchyConfig.megExpectedRows = 38:56;
hierarchyConfig.parameterVariant = 'minb1';
hierarchyConfig.parameterValueField = 'optParams';
hierarchyConfig.fixedLambda = false;

trajectoryConfig = struct();
trajectoryConfig.nParameterDraws = 30;
trajectoryConfig.nRepPerDraw = 100;
trajectoryConfig.nGrid = 101;
trajectoryConfig.maxSteps = 1000;
trajectoryConfig.b0 = 300;
trajectoryConfig.delta = 0.01;
trajectoryConfig.seed = 20260722;
trajectoryConfig.tieRule = 'original-self-priority';
trajectoryConfig.megSubjectRows = [];
trajectoryConfig.conditionRows = [];
trajectoryConfig.useParallel = true;
trajectoryConfig.nWorkers = 4;

if isfield(hierarchicalTrajectoryOverride,'parameterFile85')
    parameterFile85 = char( ...
        hierarchicalTrajectoryOverride.parameterFile85);
end
if isfield(hierarchicalTrajectoryOverride,'behaviorFile85')
    behaviorFile85 = char( ...
        hierarchicalTrajectoryOverride.behaviorFile85);
end
if isfield(hierarchicalTrajectoryOverride,'behaviorFile19')
    behaviorFile19 = char( ...
        hierarchicalTrajectoryOverride.behaviorFile19);
end
if isfield(hierarchicalTrajectoryOverride,'hierarchy')
    hierarchyConfig = merge_struct_local(hierarchyConfig, ...
        hierarchicalTrajectoryOverride.hierarchy);
end
if isfield(hierarchicalTrajectoryOverride,'trajectory')
    trajectoryConfig = merge_struct_local(trajectoryConfig, ...
        hierarchicalTrajectoryOverride.trajectory);
end

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
outputDir = fullfile(projectDir,'analysis_outputs','model_dynamics', ...
    ['hierarchical_eb_85_to_meg19_' stamp]);
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
hierarchyFile = fullfile(outputDir,'hierarchical_parameters_85.mat');
trajectoryFile = fullfile(outputDir, ...
    'hierarchical_MEG19_posterior_trajectories.mat');
trajectoryConfig.checkpointFile = trajectoryFile;

fprintf('Fitting hierarchical EB model to 85 subjects...\n');
hierarchyArgs = namedargs2cell(hierarchyConfig);
hierarchy = fit_hierarchical_eb_parameters_85( ...
    parameterFile85,behaviorFile85,behaviorFile19,hierarchyArgs{:});
hierarchy.meta.outputDirectory = outputDir;
save(hierarchyFile,'hierarchy','hierarchyConfig','-v7.3');
fprintf('Saved hierarchical parameters: %s\n',hierarchyFile);

fprintf('Generating MEG-19 posterior-predictive trajectories...\n');
trajectoryArgs = namedargs2cell(trajectoryConfig);
posteriorTrajectory = ...
    simulate_hierarchical_meg_posterior_trajectories( ...
    hierarchy,behaviorFile19,trajectoryArgs{:});
posteriorTrajectory.meta.hierarchyFile = hierarchyFile;
posteriorTrajectory.meta.outputDirectory = outputDir;
save(trajectoryFile,'posteriorTrajectory','trajectoryConfig','-v7.3');
fprintf('Saved posterior trajectories: %s\n',trajectoryFile);

%% Local helper
function destination = merge_struct_local(destination,source)
    names = fieldnames(source);
    for i = 1:numel(names)
        destination.(names{i}) = source.(names{i});
    end
end
