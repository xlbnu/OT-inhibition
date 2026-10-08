function results = main_model_compare_one(otn_data, modelType)
    assert(isscalar(modelType) && ismember(modelType,1:4), 'Invalid modelType.');
    assert(isstruct(otn_data) && isvector(otn_data) && ~isempty(otn_data), 'Invalid subjects.');
    % 1. Initialize settings
    numSubjects = length(otn_data);
    numSimulations = 100; % Number of Monte Carlo simulations 
    fixedParams = struct('b0', 300, 'delta', 0.01, 'plapse', 0.05./8);
    
    % Define four models
    modelNames = {'1_FullModel', '2_No_Lambda', '3_No_t0_tau', '4_No_Lambda_t0_tau'};
    numModels = 1;
    
    % Preallocate top-level result arrays
    results = struct('subjID', cell(1, numSubjects), ...
                     'modelFit', cell(1, numSubjects));  
                     
    totaltic = tic;
    
    % Use parfor for participant-level parallel processing
    for subjID = 1:numSubjects
        fprintf('开始处理被试 %d/%d...\n', subjID, numSubjects);
        subtic = tic;       
        subjectData = otn_data(subjID);
        fields = {'answer1','answer2','conf1','conf2','conflict','rtime2'};
        if isfield(subjectData,'inforate'), fields{end+1}='inforate';
        else, fields{end+1}='o1conf1'; end
        nRows = numel(subjectData.answer1);
        for j=1:numel(fields)
            value = subjectData.(fields{j});
            assert(isnumeric(value) || islogical(value), 'Invalid data type.');
            assert(isvector(value) && numel(value)==nRows, 'Trial shape mismatch.');
            subjectData.(fields{j}) = value(:);
        end
        
        % Keep data-extraction logic unchanged
        if isfield(subjectData,'inforate')
            valid = subjectData.answer1~=-1 & subjectData.answer2~=-1 & subjectData.conf1~=-1 & subjectData.conf2~=-1;
            conflictData = struct(...
                'selfConf', subjectData.conf1(subjectData.conflict==1 & valid==1), ...
                'otherConf', subjectData.inforate(subjectData.conflict==1 & valid==1), ...
                'initialChoice', subjectData.answer1(subjectData.conflict==1 & valid==1),...
                'finalChoice', subjectData.answer2(subjectData.conflict==1 & valid==1), ...
                'finalConf', subjectData.conf2(subjectData.conflict==1 & valid==1), ...
                'RT', subjectData.rtime2(subjectData.conflict==1 & valid==1).*1000);
        else
            valid = subjectData.answer1~=-1 & subjectData.answer2~=-1 & ~isnan(subjectData.conf1) & ~isnan(subjectData.conf2) & ~isnan(subjectData.o1conf1);
            conflictData = struct(...            
                'selfConf', subjectData.conf1(subjectData.conflict==1 & valid==1), ...
                'otherConf', subjectData.o1conf1(subjectData.conflict==1 & valid==1), ...
                'initialChoice', subjectData.answer1(subjectData.conflict==1 & valid==1),...
                'finalChoice', subjectData.answer2(subjectData.conflict==1 & valid==1), ...
                'finalConf', subjectData.conf2(subjectData.conflict==1 & valid==1), ...
                'RT', subjectData.rtime2(subjectData.conflict==1 & valid==1).*1000);            
        end
        
        assert(~isempty(conflictData.selfConf), 'No valid conflict trials.');
        checked = {'selfConf','otherConf','initialChoice','finalChoice','finalConf'};
        for j=1:numel(checked)
            assert(all(isfinite(conflictData.(checked{j}))), 'Nonfinite retained data.');
        end
        % One model per call.
        subjModelResults = struct('modelName', cell(1, numModels), ...
                                  'optParams_FullArray', cell(1, numModels), ...
                                  'nll', cell(1, numModels), ...
                                  'bic', cell(1, numModels));
                                  
        % Inner loop: Fit each current participant with 4 models
        for m = 1:numModels
            [optParams_FullArray, nll, bic] = fitModel(modelType, fixedParams, conflictData, numSimulations);
            
            subjModelResults(m).modelName = modelNames{modelType};
            subjModelResults(m).optParams_FullArray = optParams_FullArray; % Already restored to the original scale: 9 slot parameter array
            subjModelResults(m).nll = nll;
            subjModelResults(m).bic = bic;
            counts = [9 8 7 6];
            subjModelResults(m).nParams = counts(modelType);
            subjModelResults(m).nTrials = numel(conflictData.selfConf);
            subjModelResults(m).modelType = modelType;
            subjModelResults(m).b0 = fixedParams.b0;
        end
        
        % Aggregate and save
        resTemp = struct();
        resTemp.subjID = subjID;
        resTemp.modelFit = subjModelResults;
        results(subjID) = resTemp;
        
        subjTime = toc(subtic);
        fprintf('被试 %d 指定模型处理完成，耗时 %.3f 秒\n', subjID, subjTime);        
    end
    
    totalTime = toc(totaltic); 
    fprintf('\n所有被试处理完毕! 总耗时: %.3f 秒\n', totalTime);
end

%% Main optimization function (Based on modelType dynamically adjust the parameter space)
function [optParams_FullArray, nll, bic] = fitModel(modelType, fixedParams, conflictData, numSimulations)
    % Full base parameter settings
    zoom_all = [1 1 1 10 10 20 1./6 1000 1000];
    lb_all = [0.001, 0, 0, 0.001, 0.001, 0.1, 1e-5, 10, 10] ./ zoom_all; 
    ub_all = [0.99, 1, 1, 10, 10, 20, 0.15, 1000, 1000] ./ zoom_all; 
    
    % Select parameter indices to optimize by model type
    switch modelType
        case 1 % Full model (9 items)
            active_idx = [1, 2, 3, 4, 5, 6, 7, 8, 9];
        case 2 % without lambda (8 items)
            active_idx = [1, 2, 3, 4, 5, 6, 8, 9];
        case 3 % without t0, tau (7 items)
            active_idx = [1, 2, 3, 4, 5, 6, 7];
        case 4 % without lambda, t0, tau (6 items)
            active_idx = [1, 2, 3, 4, 5, 6];
    end
    
    numParams = length(active_idx);
    lb = lb_all(active_idx);
    ub = ub_all(active_idx);

    % Hybrid optimizer settings
    hybridOpts = optimoptions('patternsearch', 'Display', 'off', 'MaxIterations', 500, 'UseCompletePoll', true, 'UseParallel', false);
    pcoptions = optimoptions('particleswarm', 'SwarmSize', 120, 'MaxIterations', 500, ...
        'MaxStallIterations', 30, 'UseParallel', false, 'Display', 'off', ...
        'HybridFcn', {@patternsearch, hybridOpts});
        
    % Run particle-swarm optimization
    [optParams_active, pcNll] = particleswarm(...
        @(p) modelNLL(p, modelType, fixedParams, conflictData, numSimulations), ...
        numParams, lb, ub, pcoptions);
  
    nll = pcNll;
    
    % Reconstruct 1x9 parameter array on the original scale, For subsequent unified calls to the prediction function
    optParams_FullArray = zeros(1, 9);
    for i = 1:numParams
        orig_idx = active_idx(i);
        optParams_FullArray(orig_idx) = optParams_active(i) * zoom_all(orig_idx);
    end
    
    % Assign values to fixed parameters excluded from optimization
    if ~ismember(7, active_idx), optParams_FullArray(7) = 0; end     % lambda = 0
    if ~ismember(8, active_idx), optParams_FullArray(8) = 1e6; end   % tau = Very large value
    if ~ismember(9, active_idx), optParams_FullArray(9) = 1e6; end   % t0 = Very large value
    
    % Model evaluation (BIC Compute, Penalty based on the actual number of optimized parameters numParams Compute)
    bic = 2*nll + numParams * log(length(conflictData.selfConf));  
end

%% Compute model negative log-likelihood (Support dynamic parameter parsing)
function nll = modelNLL(active_params, modelType, fixedParams, data, numSimulations)
    zoom_all = [1 1 1 10 10 20 1./6 1000 1000];
    
    % Extract common active parameters
    v0       = active_params(1) * zoom_all(1);
    b1_ratio = active_params(2) * zoom_all(2);
    b2_ratio = active_params(3) * zoom_all(3);
    k1        = active_params(4) * zoom_all(4);
    k2        = active_params(5) * zoom_all(5);
    sigma    = active_params(6) * zoom_all(6);
    
    % Extract or set specific parameters by model type
    if modelType == 1 % Full model
        lambda = active_params(7) * zoom_all(7);
        tau    = active_params(8) * zoom_all(8);
        t0     = active_params(9) * zoom_all(9);
    elseif modelType == 2 % without lambda
        lambda = 0;
        tau    = active_params(7) * zoom_all(8);
        t0     = active_params(8) * zoom_all(9);
    elseif modelType == 3 % without t0, tau
        lambda = active_params(7) * zoom_all(7);
        tau    = 1e6;
        t0     = 1e6;
    elseif modelType == 4 % without lambda, t0, tau
        lambda = 0;
        tau    = 1e6;
        t0     = 1e6;
    end

    b0 = fixedParams.b0;
    delta = fixedParams.delta;
    plapse = fixedParams.plapse;
    numTrials = length(data.selfConf);

    % Step A and B: Extract, transform, and find unique stimulus combinations (Preserve the efficient aggregation logic)
    sti1_all = (data.selfConf) / 16 + 0.5 - 1./32 - 0.5;
    sti2_all = zeros(numTrials, 1);
    idx_high = data.otherConf >= 1;
    sti2_all(idx_high) = (data.otherConf(idx_high)) / 16 + 0.5 - 1./32 - 0.5;
    sti2_all(~idx_high) = data.otherConf(~idx_high) - 0.5;
    
    stim_pairs = [sti1_all, sti2_all];
    [unique_stims, ~, ic] = unique(stim_pairs, 'rows');
    numConditions = size(unique_stims, 1);
    
    condProbMatrix = zeros(2, 4, numConditions);
    validConditions = 0;
    simsPerCond = max(numSimulations, 200); 

    % Step C: Condition-wise simulation
    for c = 1:numConditions
        sti1_c = unique_stims(c, 1);
        sti2_c = unique_stims(c, 2);
        
        [choices, confs, ~, valids] = runDDM_vectorized(v0, b1_ratio, b2_ratio, k1, k2, sigma, ...
                                        lambda, tau, t0, delta, b0, plapse, sti1_c, sti2_c, simsPerCond);
        
        if sum(valids) == 0
            continue;
        end
        
        valid_choices = choices(valids);
        valid_confs = confs(valids);
        confLevels = discretize_conf_vectorized(valid_confs);
        
        counts = accumarray([valid_choices, confLevels], 1, [2, 4]);
        matrix_c = counts + 1e-5;
        condProbMatrix(:,:,c) = matrix_c / sum(matrix_c, 'all');
        validConditions = validConditions + 1;
    end

    if validConditions == 0
        nll = 1e10; return;
    end

    % Step D: Vectorized lookup calculation Log-Likelihood
    obsChoice_all = (data.finalChoice ~= data.initialChoice) + 1;
    obsConf_all = discretize_conf_vectorized((data.finalConf) / 16 - 1/32 + 0.5);
    
    logLikelihood = 0;
    for t = 1:numTrials
        cond_idx = ic(t); 
        p_matrix = condProbMatrix(:,:,cond_idx); 
        prob = p_matrix(obsChoice_all(t), obsConf_all(t));
        adjusted_prob = prob * (1 - plapse * 8) + plapse;
        logLikelihood = logLikelihood + log(adjusted_prob);
    end
    nll = -logLikelihood;
end

%% Vectorized λStructurally enhanced DDM model (Keep unchanged; no modifications needed)
function [choices, confs, rts, valids] = runDDM_vectorized(v0, b1_ratio, b2_ratio, k1, k2, sigma, ...
                                            lambda, tau, t0, delta, b0, plapse, sti1, sti2, numSimulations)
    b1 = b0 * max(0, b1_ratio);
    b2 = b0 * max(0, b2_ratio);
    v1 = v0 + k1 * sti1;
    v2 = v0 + k2 * sti2;
    x1 = repmat(b1, numSimulations, 1);
    x2 = repmat(b2, numSimulations, 1);
    
    choices = ones(numSimulations, 1);
    rts = zeros(numSimulations, 1);
    active = true(numSimulations, 1); 
    
    dt = 10;
    maxSteps = 1000;
    
    steps_vec = 1:maxSteps;
    boundaries = b0 * ones(1, maxSteps);
    idx_t0 = steps_vec >= t0;
    boundaries(idx_t0) = b0 * (1 - 0.9 * (steps_vec(idx_t0) - t0) ./ (steps_vec(idx_t0) - t0 + tau));
    
    for step = 1:maxSteps
        if ~any(active)
            break; 
        end
        n_active = sum(active);
        curr_bound = boundaries(step);
        
        noise1 = sigma * randn(n_active, 1);
        noise2 = sigma * randn(n_active, 1);
        
        cx1 = x1(active);
        cx2 = x2(active);
        
        inh1 = -lambda * cx2;
        inh2 = -lambda * cx1;   
        dx1 = inh1 + (v1 - delta * cx1) + noise1;
        dx2 = inh2 + (v2 - delta * cx2) + noise2;
        nx1 = max(0, cx1 + dx1);
        nx2 = max(0, cx2 + dx2);
        
        x1(active) = nx1;
        x2(active) = nx2;
        crossed1 = nx1 >= curr_bound;
        crossed2 = nx2 >= curr_bound & ~crossed1; 
        crossed_any = crossed1 | crossed2;
        if any(crossed_any)
            active_idx = find(active);
            just_crossed_global = active_idx(crossed_any);
            
            choices(active_idx(crossed1)) = 1;
            choices(active_idx(crossed2)) = 2;
            rts(just_crossed_global) = step * dt;
            
            active(just_crossed_global) = false;
        end
    end
    if any(active)
        choices(active) = randi([1, 2], sum(active), 1);
        rts(active) = maxSteps * dt;
    end
    
    valids = ~(isnan(x1) | isnan(x2) | isinf(x1) | isinf(x2));
    e1 = min(0.5 + x1./(2*b0), 1 - 8*plapse);
    e2 = min(0.5 + x2./(2*b0), 1 - 8*plapse);
    e1 = max(min(e1, 1-eps), eps);
    e2 = max(min(e2, 1-eps), eps);
    log_odds1 = log(e1) - log(1-e1);
    log_odds2 = log(e2) - log(1-e2);
    
    confs = zeros(numSimulations, 1);
    idx1 = (choices == 1);
    idx2 = (choices == 2);
    
    confs(idx1) = 1 ./ (1 + exp(-(log_odds1(idx1) - log_odds2(idx1))));
    confs(idx2) = 1 ./ (1 + exp(-(log_odds2(idx2) - log_odds1(idx2))));
end

%% Discretize confidence (Vectorized version)
function confLevels = discretize_conf_vectorized(confs)
    confs = max(min(confs, 1), 0);
    confLevels = ones(size(confs)); 
    confLevels(confs >= 0.625 & confs < 0.75) = 2;
    confLevels(confs >= 0.75 & confs < 0.875) = 3;
    confLevels(confs >= 0.875) = 4;
end
