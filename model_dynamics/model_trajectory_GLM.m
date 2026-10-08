%% 
% Run time-resolved trajectory GLMs.
% Predictors are aligned through each bank's sourceTrialIndex. At every
% endpoint step, a trial enters a GLM only when its trajectory is reliable,
% finite, and its behavioral predictors are finite. Stay/change models add
% the corresponding observed-choice filter.

%%
codeRoot = fileparts(fileparts(mfilename('fullpath')));
outputRoot = fullfile(codeRoot,'analysis_outputs','model_dynamics');
% Set either input filename before running to select a specific exported run.
if ~exist('fullTrajectoryFile','var') || isempty(fullTrajectoryFile)
    fullTrajectoryFile = latest_export_local(outputRoot, ...
        'trial_accumulator_GLM_*','MEG_trial_accumulator_GLM_data.mat');
end
if ~exist('noLambdaTrajectoryFile','var') || isempty(noLambdaTrajectoryFile)
    noLambdaTrajectoryFile = latest_export_local(outputRoot, ...
        'no_lambda_trial_accumulator_hierarchical85_*','no_lambda_trial_accumulator_GLM_data.mat');
end
load(fullTrajectoryFile,'trialAccumulator');
load(noLambdaTrajectoryFile,'trialAccumulatorNoLambda');
load(fullfile(codeRoot,'figure_data','behavioral_results','bhv_data_MEG_19.mat'), ...
    'ots_data','pls_data','otn_data','pln_data');
if ~isfolder(outputRoot), mkdir(outputRoot); end

%%
data_tmp1=trialAccumulatorNoLambda.subject;
data_tmp=trialAccumulator.subject;
data_bhv={ots_data,pls_data,otn_data,pln_data};
beta_acc=nan(4,19,1000,3);beta_acc1=nan(4,19,1000,3);
for i_sub=1:length(data_tmp)
for i_condition = 1:4
Y = data_tmp{i_sub,i_condition}.end.choiceAligned.raw.B0.difference;
Y1 = data_tmp1{i_sub,i_condition}.end.choiceAligned.raw.B0.difference;
X1 = data_bhv{i_condition}(i_sub).lo_cplus(trialAccumulator.subject{i_sub,i_condition}.trial.sourceTrialIndex);
X2 = data_bhv{i_condition}(i_sub).lo_cminus(trialAccumulator.subject{i_sub,i_condition}.trial.sourceTrialIndex);


for i_step=1:size(Y,2)-1
    if sum(~isnan(Y(:,i_step)))>30
        valid = ~isnan(Y(:,i_step));
        b_tmp=glmfit(zscore([X1(valid),X2(valid)]),zscore(Y(valid,i_step)));
    else
        b_tmp=nan(3,1);
    end
    beta_acc(i_condition,i_sub,i_step,:)=b_tmp;

    if sum(~isnan(Y1(:,i_step)))>30
        valid = ~isnan(Y1(:,i_step));
        b_tmp=glmfit(zscore([X1(valid),X2(valid)]),zscore(Y1(valid,i_step)));
    else
        b_tmp=nan(3,1);
    end
    beta_acc1(i_condition,i_sub,i_step,:)=b_tmp;    
end
end
end

% Compute beta(self) & beta(other)
data_tmp1=trialAccumulatorNoLambda.subject;
data_tmp=trialAccumulator.subject;
data_bhv={ots_data,pls_data,otn_data,pln_data};
bs_acc=nan(4,19,1000,2);bs_acc1=nan(4,19,1000,2);
bc_acc=nan(4,19,1000,2);bc_acc1=nan(4,19,1000,2);
for i_sub=1:length(data_tmp)
for i_condition = 1:4
Y_s = data_tmp{i_sub,i_condition}.end.choiceAligned.raw.B0.self;
Y_o = data_tmp{i_sub,i_condition}.end.choiceAligned.raw.B0.other;
sourceRows = data_tmp{i_sub,i_condition}.trial.sourceTrialIndex;
X1 = data_bhv{i_condition}(i_sub).lo_s1(sourceRows);
X2 = data_bhv{i_condition}(i_sub).lo_o1(sourceRows);
ischange = double(data_tmp{i_sub,i_condition}.trial.ischange);

Y1_s = data_tmp1{i_sub,i_condition}.end.choiceAligned.raw.B0.self;
Y1_o = data_tmp1{i_sub,i_condition}.end.choiceAligned.raw.B0.other;
sourceRows1 = data_tmp1{i_sub,i_condition}.trial.sourceTrialIndex;
X12 = data_bhv{i_condition}(i_sub).lo_s1(sourceRows1);
X22 = data_bhv{i_condition}(i_sub).lo_o1(sourceRows1);
ischange1 = double(data_tmp1{i_sub,i_condition}.trial.ischange);


for i_step=1:size(Y_s,2)-1
    if sum(~isnan(Y_s(:,i_step)))>30
        valids = ~isnan(Y_s(:,i_step)) & ischange==0;
        bs_tmps=glmfit(zscore([X1(valids)]),(Y_s(valids,i_step)));
    else
        bs_tmps=nan(2,1);
    end
    if sum(~isnan(Y_o(:,i_step)))>30
        valids = ~isnan(Y_o(:,i_step)) & ischange==0;
        bs_tmpo=glmfit(zscore([X2(valids)]),(Y_o(valids,i_step)));
    else
        bs_tmpo=nan(2,1);
    end

    if sum(~isnan(Y_s(:,i_step)))>30
        valids = ~isnan(Y_s(:,i_step)) & ischange==1;
        bc_tmps=glmfit(zscore([X1(valids)]),(Y_s(valids,i_step)));
    else
        bc_tmps=nan(2,1);
    end
    if sum(~isnan(Y_o(:,i_step)))>30
        valids = ~isnan(Y_o(:,i_step)) & ischange==1;
        bc_tmpo=glmfit(zscore([X2(valids)]),(Y_o(valids,i_step)));
    else
        bc_tmpo=nan(2,1);
    end
    bs_acc(i_condition,i_sub,i_step,1)=bs_tmps(2);bs_acc(i_condition,i_sub,i_step,2)=bs_tmpo(2);
    bc_acc(i_condition,i_sub,i_step,1)=bc_tmps(2);bc_acc(i_condition,i_sub,i_step,2)=bc_tmpo(2);

    % no lambda

    if sum(~isnan(Y1_s(:,i_step)))>30
        valids = ~isnan(Y1_s(:,i_step)) & ischange1==0;
        bs_tmp1s=glmfit(zscore([X12(valids)]),(Y1_s(valids,i_step)));
    else
        bs_tmp1s=nan(2,1);
    end
    if sum(~isnan(Y1_o(:,i_step)))>30
        valids = ~isnan(Y1_o(:,i_step)) & ischange1==0;
        bs_tmp1o=glmfit(zscore([X22(valids)]),(Y1_o(valids,i_step)));
    else
        bs_tmp1o=nan(2,1);
    end

    if sum(~isnan(Y1_s(:,i_step)))>30
        valids = ~isnan(Y1_s(:,i_step)) & ischange1==1;
        bc_tmp1s=glmfit(zscore([X12(valids)]),(Y1_s(valids,i_step)));
    else
        bc_tmp1s=nan(2,1);
    end
    if sum(~isnan(Y1_o(:,i_step)))>30
        valids = ~isnan(Y1_o(:,i_step)) & ischange1==1;
        bc_tmp1o=glmfit(zscore([X22(valids)]),(Y1_o(valids,i_step)));
    else
        bc_tmp1o=nan(2,1);
    end

    bs_acc1(i_condition,i_sub,i_step,1)=bs_tmp1s(2);bs_acc1(i_condition,i_sub,i_step,2)=bs_tmp1o(2);       
    bc_acc1(i_condition,i_sub,i_step,1)=bc_tmp1s(2);bc_acc1(i_condition,i_sub,i_step,2)=bc_tmp1o(2);
end
end
end



%%
save(fullfile(outputRoot,'model_DDM_trajectory_GLM.mat'),'beta_acc','beta_acc1','bs_acc','bc_acc','bs_acc1','bc_acc1');


function filename = latest_export_local(outputRoot,pattern,name)
    candidates = dir(fullfile(outputRoot,pattern,name));
    assert(~isempty(candidates), ...
        'Missing %s. Run the corresponding trial-accumulator export first.',name);
    [~,index] = max([candidates.datenum]);
    filename = fullfile(candidates(index).folder,candidates(index).name);
end
