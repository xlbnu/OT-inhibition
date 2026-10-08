function results = main_model(otn_data)
    % 1. Initialize settings
    assert(isstruct(otn_data) && isvector(otn_data) && ~isempty(otn_data), ...
        'Input must be a nonempty subject struct vector.');
    numSubjects = numel(otn_data);
    numSimulations = 100; % Number of Monte Carlo simulations (Now vectorized, Even when increased to 500 execution remains fast)
    fixedParams = struct('b0', 300, 'delta', 0.01, 'plapse', 0.05./8);
    
    % Preallocate result arrays, Ensure parfor compatibility
    results = struct('optParams', cell(numSubjects, 1), ...
                     'nll', cell(numSubjects, 1), ...
                     'bic', cell(numSubjects, 1));
                     
    totaltic = tic;
    
    % Subject-level parallelism and checkpoints belong to run_slurm_fit.
    for subjID = 1:numSubjects
        fprintf('开始处理被试 %d/%d...\n', subjID, numSubjects);
        subtic = tic;       
        subjectData = otn_data(subjID);
        fields = {'answer1','answer2','conf1','conf2','conflict','rtime2'};
        if isfield(subjectData,'inforate')
            fields{end+1} = 'inforate';
        else
            fields{end+1} = 'o1conf1';
        end
        nRows = [];
        for j = 1:numel(fields)
            name = fields{j};
            assert(isfield(subjectData,name), 'Missing field: %s', name);
            value = subjectData.(name);
            assert((isnumeric(value) || islogical(value)) && ...
                (isvector(value) || isempty(value)), 'Invalid vector: %s', name);
            if isempty(nRows), nRows = numel(value); end
            assert(numel(value) == nRows, 'Trial length mismatch: %s', name);
            subjectData.(name) = value(:);
        end
        
        % Select conflict trials (Initial choice ≠ Other-person choice)
        if isfield(subjectData,'inforate')
            valid=subjectData.answer1~=-1 & subjectData.answer2~=-1 & subjectData.conf1~=-1 & subjectData.conf2~=-1;
            conflictData = struct(...
                'selfConf', subjectData.conf1(subjectData.conflict==1 & valid==1), ...
                'otherConf', subjectData.inforate(subjectData.conflict==1 & valid==1), ...
                'initialChoice', subjectData.answer1(subjectData.conflict==1 & valid==1),...
                'finalChoice', subjectData.answer2(subjectData.conflict==1 & valid==1), ...
                'finalConf', subjectData.conf2(subjectData.conflict==1 & valid==1), ...
                'RT', subjectData.rtime2(subjectData.conflict==1 & valid==1).*1000);
        else
            valid=subjectData.answer1~=-1 & subjectData.answer2~=-1 & ~isnan(subjectData.conf1) & ~isnan(subjectData.conf2) & ~isnan(subjectData.o1conf1);
            conflictData = struct(...            
                'selfConf', subjectData.conf1(subjectData.conflict==1 & valid==1), ...
                'otherConf', subjectData.o1conf1(subjectData.conflict==1 & valid==1), ...
                'initialChoice', subjectData.answer1(subjectData.conflict==1 & valid==1),...
                'finalChoice', subjectData.answer2(subjectData.conflict==1 & valid==1), ...
                'finalConf', subjectData.conf2(subjectData.conflict==1 & valid==1), ...
                'RT', subjectData.rtime2(subjectData.conflict==1 & valid==1).*1000);            
        end
        
        assert(~isempty(conflictData.selfConf), 'No valid conflict trials.');
        names = {'selfConf','otherConf','initialChoice','finalChoice','finalConf'};
        for j = 1:numel(names)
            value = conflictData.(names{j});
            assert(all(isfinite(value)), 'Nonfinite retained data: %s', names{j});
        end

        initGuess = [0.1, 1e-5, 1e-5, 0.5, 0.5, 1, 0, 100, 100] ./ [1 1 1 10 10 20 1./6 1000 1000];
        
        % Call the optimization function
        [optParams, nll, bic] = fitModel(initGuess, fixedParams, conflictData, numSimulations);

        % Save results (Note: in parfor assign the complete structure element this way)
        resTemp = struct();
        resTemp.optParams = optParams .* [1 1 1 10 10 20 1./6 1000 1000];
        resTemp.nll = nll;
        resTemp.bic = bic;
        results(subjID) = resTemp;
        
        subjTime = toc(subtic);
        fprintf('被试 %d 处理完成，耗时 %.3f 秒\n', subjID, subjTime);        
    end
    
    totalTime = toc(totaltic);  % Get total elapsed time
    fprintf('\n所有被试处理完毕! 总耗时: %.3f 秒\n', totalTime);
end

%% Improved main optimization function (including PSO + Patternsearch hybrid optimization)
function [optParams, nll, bic] = fitModel(params, fixedParams, conflictData, numSimulations)
    stage1Params = params;
    
    % Add parameter bounds
    zoom_size = [1 1 1 10 10 20 1./6 1000 1000];
    lb = [0.001, 0, 0, 0.001, 0.001, 0.1, 1e-5, 10, 10] ./ zoom_size; 
    ub = [0.99, 1, 1, 10, 10, 20, 0.15, 1000, 1000] ./ zoom_size; 

    % Recommended 1 & 2: configuration patternsearch Hybrid optimizer, Handle stochastic noise
    hybridOpts = optimoptions('patternsearch', ...
        'Display', 'off', ...
        'MaxIterations', 500, ...
        'UseCompletePoll', true);

    % Recommended 3: Adjust PSO parameters, Disable internal parallelism, Reduce iteration count, Enable hybrid optimization
    pcoptions = optimoptions('particleswarm', ...
        'SwarmSize', 120, ...                % Reduce swarm size (10*D)
        'MaxIterations', 500, ...           % Reduce maximum generations, Let the local optimizer take over
        'MaxStallIterations', 30, ...       % Early stopping criterion
        'UseParallel', false, ...           % Disable lower-level parallelism, Leave resources for the outer parfor loop
        'Display', 'off', ...
        'HybridFcn', {@patternsearch, hybridOpts});

    % Run particle-swarm optimization
    [stage1Opt_pc, pcNll] = particleswarm(...
        @(p) modelNLL(p, fixedParams, conflictData, numSimulations), ...
        length(stage1Params), lb, ub, pcoptions);
  
    optParams = stage1Opt_pc;
    nll = pcNll;

    % Model evaluation (BIC Compute)
    % bic = nll + 0.5 * length(params) * log(length(conflictData.selfConf));
    nParams = numel(params);  % currently 9
    nTrials = numel(conflictData.selfConf);
    bic = 2*nll + nParams*log(nTrials);
end

%% Compute model negative log-likelihood (Optimized aggregation by condition)
function nll = modelNLL(params, fixedParams, data, numSimulations)
    % 1. Unpack parameters
    zoom_size = [1 1 1 10 10 20 1./6 1000 1000];
    v0       = params(1) * zoom_size(1);
    b1_ratio = params(2) * zoom_size(2);
    b2_ratio = params(3) * zoom_size(3);
    k1        = params(4) * zoom_size(4);
    k2        = params(5) * zoom_size(5);
    sigma    = params(6) * zoom_size(6);
    lambda   = params(7) * zoom_size(7);
    tau      = params(8) * zoom_size(8);
    t0       = params(9) * zoom_size(9);

    b0 = fixedParams.b0;
    delta = fixedParams.delta;
    plapse = fixedParams.plapse;

    numTrials = length(data.selfConf);

    % ==========================================
    % Step A: Extract and transform stimulus strengths for all trials
    % ==========================================
    sti1_all = (data.selfConf) / 16 + 0.5 - 1./32 - 0.5;
    
    sti2_all = zeros(numTrials, 1);
    idx_high = data.otherConf >= 1;
    sti2_all(idx_high) = (data.otherConf(idx_high)) / 16 + 0.5 - 1./32 - 0.5;
    sti2_all(~idx_high) = data.otherConf(~idx_high) - 0.5;

    % ==========================================
    % Step B: Find unique stimulus combinations (Group by Condition)
    % ==========================================
    stim_pairs = [sti1_all, sti2_all];
    % unique_stims: Matrix of unique combinations;  ic: Indices mapping original trials to unique combinations
    [unique_stims, ~, ic] = unique(stim_pairs, 'rows');
    numConditions = size(unique_stims, 1);

    % Preallocate a 3D matrix, Store the predicted distribution for each condition: 2(Choice) x 4(Conf) x numConditions
    condProbMatrix = zeros(2, 4, numConditions);
    validConditions = 0;

    % To ensure smooth condition distributions, Because the total simulation count has decreased, the simulations per condition can be increased appropriately
    % For example: ensure each condition has at least 200 particles running
    simsPerCond = max(numSimulations, 200); 

    % ==========================================
    % Step C: Iterate only over[unique stimulus conditions]Run simulations
    % ==========================================
    for c = 1:numConditions
        sti1_c = unique_stims(c, 1);
        sti2_c = unique_stims(c, 2);
        
        % Run vectorized DDM
        [choices, confs, ~, valids] = runDDM_vectorized(v0, b1_ratio, b2_ratio, k1 , k2, sigma, ...
                                        lambda, tau, t0, delta, b0, plapse, sti1_c, sti2_c, simsPerCond);
        
        if sum(valids) == 0
            continue;
        end
        
        valid_choices = choices(valids);
        valid_confs = confs(valids);
        confLevels = discretize_conf_vectorized(valid_confs);
        
        % Count the joint-distribution matrix
        counts = accumarray([valid_choices, confLevels], 1, [2, 4]);
        
        % Add a small baseline value to prevent log(0), and normalize probabilities
        matrix_c = counts + 1e-5;
        condProbMatrix(:,:,c) = matrix_c / sum(matrix_c(:));
        validConditions = validConditions + 1;
    end

    % If simulation fails for all conditions(extreme parameters), return a large penalty
    if validConditions == 0
        nll = 1e10;
        return;
    end

    % ==========================================
    % Step D: Vectorized lookup calculation Log-Likelihood
    % ==========================================
    % Preprocess all observed choices and confidence ratings
    obsChoice_all = (data.finalChoice ~= data.initialChoice) + 1;
    obsConf_all = discretize_conf_vectorized((data.finalConf) / 16 - 1/32 + 0.5);
    
    logLikelihood = 0;
    
    % This loop only accumulates lookup values, No model simulations are performed, Very fast
    for t = 1:numTrials
        cond_idx = ic(t); % Find the condition-matrix index for the current trial
        p_matrix = condProbMatrix(:,:,cond_idx); % Retrieve the predicted probability matrix
        
        % Look up the predicted probability of the observed trial outcome
        prob = p_matrix(obsChoice_all(t), obsConf_all(t));
        
        % Add lapse rate (lapse rate)
        adjusted_prob = prob * (1 - plapse * 8) + plapse;
        
        logLikelihood = logLikelihood + log(adjusted_prob);
    end

    nll = -logLikelihood;
end
%% Vectorized λStructurally enhanced DDM model
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
    
    % Precompute dynamic boundaries for all time steps
    steps_vec = 1:maxSteps;
    boundaries = b0 * ones(1, maxSteps);
    idx_t0 = steps_vec >= t0;
    boundaries(idx_t0) = b0 * (1 - 0.9 * (steps_vec(idx_t0) - t0) ./ (steps_vec(idx_t0) - t0 + tau));

    % Batch evidence accumulation
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




