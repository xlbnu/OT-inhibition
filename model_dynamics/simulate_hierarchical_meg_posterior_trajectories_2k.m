function posteriorTrajectory = ...
        simulate_hierarchical_meg_posterior_trajectories_2k( ...
        hierarchy,behaviorFile19,cfg)
%SIMULATE_HIERARCHICAL_MEG_POSTERIOR_TRAJECTORIES Propagate EB uncertainty.
% Posterior parameter draws are combined with DDM path simulations. The
% pooled mean uses all parameter-draw/path combinations. Pointwise 95%
% intervals describe variation across parameter draws and their simulated
% paths; they are approximate posterior-predictive intervals.

    arguments
        hierarchy (1,1) struct
        behaviorFile19 (1,:) char
        cfg.nParameterDraws (1,1) double = 30
        cfg.nRepPerDraw (1,1) double = 100
        cfg.nGrid (1,1) double = 101
        cfg.maxSteps (1,1) double = 1000
        cfg.b0 (1,1) double = 300
        cfg.delta (1,1) double = 0.01
        cfg.seed (1,1) double = 20260722
        cfg.tieRule (1,:) char = 'original-self-priority'
        cfg.megSubjectRows (1,:) double = []
        cfg.conditionRows (1,:) double = []
        cfg.checkpointFile (1,:) char = ''
        cfg.useParallel (1,1) logical = false
        cfg.nWorkers (1,1) double = 4
    end

    allParameterDraws = hierarchy.meg.parameters.draws;
    nAvailableSubjects = size(allParameterDraws,1);
    nAvailableConditions = size(allParameterDraws,2);
    nAvailableDraws = size(allParameterDraws,4);
    if isempty(cfg.megSubjectRows)
        cfg.megSubjectRows = 1:nAvailableSubjects;
    end
    if isempty(cfg.conditionRows)
        cfg.conditionRows = 1:nAvailableConditions;
    end
    validate_rows_local(cfg.megSubjectRows,nAvailableSubjects, ...
        'cfg.megSubjectRows');
    validate_rows_local(cfg.conditionRows,nAvailableConditions, ...
        'cfg.conditionRows');
    assert(cfg.nParameterDraws >= 2 && ...
        cfg.nParameterDraws <= nAvailableDraws, ...
        'cfg.nParameterDraws must be between 2 and %d.',nAvailableDraws);
    assert(cfg.nRepPerDraw >= 1 && cfg.nRepPerDraw == ...
        round(cfg.nRepPerDraw),'cfg.nRepPerDraw must be a positive integer.');

    drawRows = unique(round(linspace( ...
        1,nAvailableDraws,cfg.nParameterDraws)),'stable');
    assert(numel(drawRows) == cfg.nParameterDraws, ...
        'Could not select the requested number of unique parameter draws.');

    d = load(behaviorFile19, ...
        'ots_data','pls_data','otn_data','pln_data');
    dataSets = {d.ots_data,d.pls_data,d.otn_data,d.pln_data};
    isSocial = [true true false false];
    fullConditionNames = {'OT-social','PL-social', ...
        'OT-nonsocial','PL-nonsocial'};
    subjectRows = cfg.megSubjectRows;
    conditionRows = cfg.conditionRows;
    nSubjects = numel(subjectRows);
    nConditions = numel(conditionRows);

    posteriorTrajectory = initialize_output_local( ...
        nSubjects,nConditions,cfg);
    posteriorTrajectory.meta.created = datetime('now');
    posteriorTrajectory.meta.method = ...
        'two-k hierarchical empirical-Bayes posterior-predictive trajectories';
    posteriorTrajectory.meta.parameterNames = hierarchy.meta.parameterNames;
    posteriorTrajectory.meta.warning = hierarchy.meta.warning;
    posteriorTrajectory.meta.behaviorFile19 = behaviorFile19;
    posteriorTrajectory.meta.conditionNames = ...
        fullConditionNames(conditionRows);
    posteriorTrajectory.meta.conditionRows = conditionRows;
    posteriorTrajectory.meta.subjectRowsInMEG19 = subjectRows;
    posteriorTrajectory.meta.subjectRowsIn85 = ...
        hierarchy.meg.rowsIn85(subjectRows);
    posteriorTrajectory.meta.subjectIDs = ...
        hierarchy.meg.subjectIDs(subjectRows);
    posteriorTrajectory.meta.parameterDrawRows = drawRows;
    posteriorTrajectory.meta.pathsPerTrial = ...
        cfg.nParameterDraws*cfg.nRepPerDraw;
    posteriorTrajectory.config = cfg;

    nJobs = nSubjects*nConditions;
    jobSubject = repmat((1:nSubjects)',nConditions,1);
    jobCondition = repelem((1:nConditions)',nSubjects);
    jobSourceSubject = subjectRows(jobSubject);
    jobSourceCondition = conditionRows(jobCondition);
    jobResult = cell(nJobs,1);

    runParallel = cfg.useParallel && ...
        license('test','Distrib_Computing_Toolbox') && ...
        exist('parpool','file') == 2;
    if cfg.useParallel && ~runParallel
        warning('Parallel toolbox unavailable; using serial trajectory simulation.');
    end
    if runParallel
        start_parallel_pool_local(cfg.nWorkers);
        parfor job = 1:nJobs
            jobResult{job} = simulate_one_subject_condition_local( ...
                allParameterDraws,drawRows,jobSourceSubject(job), ...
                jobSourceCondition(job),dataSets,isSocial,cfg);
        end
    else
        for job = 1:nJobs
            jobResult{job} = simulate_one_subject_condition_local( ...
                allParameterDraws,drawRows,jobSourceSubject(job), ...
                jobSourceCondition(job),dataSets,isSocial,cfg);
        end
    end

    for job = 1:nJobs
        subjectIndex = jobSubject(job);
        conditionIndex = jobCondition(job);
        result = jobResult{job};
        fprintf(['Hierarchical trajectory: condition %d/%d (%s), ', ...
            'MEG subject %d/%d\n'],conditionIndex,nConditions, ...
            fullConditionNames{jobSourceCondition(job)},subjectIndex,nSubjects);
        posteriorTrajectory = assign_curves_local( ...
            posteriorTrajectory,result.pooledResult,result.drawResult, ...
            subjectIndex,conditionIndex);
        posteriorTrajectory.qc.nEmpiricalTrials( ...
            subjectIndex,conditionIndex) = uint16(result.nTrials);
        posteriorTrajectory.qc.nRequestedPaths( ...
            subjectIndex,conditionIndex) = uint32( ...
            result.requestedPaths);
        posteriorTrajectory.qc.nValidPaths( ...
            subjectIndex,conditionIndex) = uint32(result.validPaths);
        posteriorTrajectory.qc.nUnresolvedPaths( ...
            subjectIndex,conditionIndex) = uint32( ...
            result.requestedPaths-result.validPaths);
        posteriorTrajectory.qc.validFraction( ...
            subjectIndex,conditionIndex) = single( ...
            result.validPaths/result.requestedPaths);
        posteriorTrajectory.completed(subjectIndex,conditionIndex) = true;
        if ~isempty(cfg.checkpointFile)
            save(cfg.checkpointFile,'posteriorTrajectory','-v7.3');
        end
    end
    posteriorTrajectory.meta.completed = datetime('now');
end

function result = simulate_one_subject_condition_local( ...
        allParameterDraws,drawRows,sourceSubject,sourceCondition, ...
        dataSets,isSocial,cfg)
    trials = extract_conflict_trials_local( ...
        dataSets{sourceCondition}(sourceSubject), ...
        isSocial(sourceCondition));
    nTrials = numel(trials.stiSelf);
    assert(nTrials > 0, ...
        'No valid conflict trials for MEG subject %d, condition %d.', ...
        sourceSubject,sourceCondition);

    pooled = [];
    drawResult = cell(cfg.nParameterDraws,1);
    for drawIndex = 1:cfg.nParameterDraws
        sourceDraw = drawRows(drawIndex);
        params = double(reshape(allParameterDraws( ...
            sourceSubject,sourceCondition,:,sourceDraw),1,[]));
        seed = cfg.seed + 10000000*sourceCondition + ...
            10000*sourceSubject + drawIndex;
        summary = simulate_batch_summary_local( ...
            params,trials.stiSelf,trials.stiOther, ...
            cfg,cfg.nRepPerDraw,seed);
        pooled = add_batch_summary_local(pooled,summary);
        drawResult{drawIndex} = finalize_summary_local(summary);
    end
    result.pooledResult = finalize_summary_local(pooled);
    result.drawResult = drawResult;
    result.nTrials = nTrials;
    result.requestedPaths = sum(pooled.requestedCount);
    result.validPaths = sum(pooled.validCount);
end

function start_parallel_pool_local(nWorkers)
    assert(nWorkers >= 1 && nWorkers == round(nWorkers), ...
        'cfg.nWorkers must be a positive integer when parallel mode is enabled.');
    pool = gcp('nocreate');
    if isempty(pool)
        parpool('Processes',nWorkers);
    elseif pool.NumWorkers ~= nWorkers
        warning(['Using the existing pool with %d workers; requested %d. ', ...
            'Delete the existing pool first to change its size.'], ...
            pool.NumWorkers,nWorkers);
    end
end

function output = initialize_output_local(nSubjects,nConditions,cfg)
    output.process101.progress = linspace(0,1,cfg.nGrid);
    output.end.step = -cfg.maxSteps:0;
    output.completed = false(nSubjects,nConditions);
    output.qc.nEmpiricalTrials = zeros( ...
        nSubjects,nConditions,'uint16');
    output.qc.nRequestedPaths = zeros( ...
        nSubjects,nConditions,'uint32');
    output.qc.nValidPaths = zeros( ...
        nSubjects,nConditions,'uint32');
    output.qc.nUnresolvedPaths = zeros( ...
        nSubjects,nConditions,'uint32');
    output.qc.validFraction = nan( ...
        nSubjects,nConditions,'single');

    branches = {'stay','change'};
    normalizations = {'B0','boundary'};
    accumulators = {'self','other'};
    statistics = {'mean','lower95','upper95'};
    for b = 1:numel(branches)
        branch = branches{b};
        for n = 1:numel(normalizations)
            normalization = normalizations{n};
            for a = 1:numel(accumulators)
                accumulator = accumulators{a};
                for q = 1:numel(statistics)
                    statistic = statistics{q};
                    output.process101.(branch).(normalization). ...
                        (accumulator).(statistic) = nan( ...
                        nSubjects,nConditions,cfg.nGrid,'single');
                    output.end.(branch).(normalization). ...
                        (accumulator).(statistic) = nan( ...
                        nSubjects,nConditions,cfg.maxSteps+1,'single');
                end
            end
        end
    end
    roles = {'winner','loser'};
    for r = 1:numel(roles)
        role = roles{r};
        for n = 1:numel(normalizations)
            normalization = normalizations{n};
            for q = 1:numel(statistics)
                statistic = statistics{q};
                output.process101.(role).(normalization). ...
                    (statistic) = nan( ...
                    nSubjects,nConditions,cfg.nGrid,'single');
                output.end.(role).(normalization).(statistic) = nan( ...
                    nSubjects,nConditions,cfg.maxSteps+1,'single');
            end
        end
    end
    output.process101.stay.weight = zeros( ...
        nSubjects,nConditions,'single');
    output.process101.change.weight = zeros( ...
        nSubjects,nConditions,'single');
    output.process101.stayProbability = nan( ...
        nSubjects,nConditions,'single');
    output.process101.changeProbability = nan( ...
        nSubjects,nConditions,'single');
    output.process101.nPaths = zeros( ...
        nSubjects,nConditions,'uint32');
    output.end.stay.weight = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'single');
    output.end.change.weight = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'single');
    output.end.stay.nTrials = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'uint16');
    output.end.change.nTrials = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'uint16');
    output.end.nTrials = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'uint16');
    output.end.nPaths = zeros( ...
        nSubjects,nConditions,cfg.maxSteps+1,'uint32');
end

function summary = simulate_batch_summary_local(params,stiSelfTrial, ...
        stiOtherTrial,cfg,nRep,seed)
    nTrials = numel(stiSelfTrial);
    nPaths = nTrials*nRep;
    nNativeSteps = cfg.maxSteps+1;
    nGroups = 2*nTrials;
    trialId = repelem((1:nTrials)',nRep);
    stiSelf = repelem(stiSelfTrial(:),nRep);
    stiOther = repelem(stiOtherTrial(:),nRep);

    assert(numel(params) == 9 && all(isfinite(params)), ...
        'The two-k trajectory simulator requires 9 finite parameters.');
    v0 = params(1);
    xSelf = cfg.b0*max(0,params(2))*ones(nPaths,1);
    xOther = cfg.b0*max(0,params(3))*ones(nPaths,1);
    kSelf = params(4);
    kOther = params(5);
    sigma = params(6);
    lambda = params(7);
    tau = params(8);
    t0 = params(9);
    vSelf = v0+kSelf*stiSelf;
    vOther = v0+kOther*stiOther;

    steps = 1:cfg.maxSteps;
    boundaries = cfg.b0*ones(1,cfg.maxSteps);
    collapse = steps >= t0;
    boundaries(collapse) = cfg.b0*(1-0.9* ...
        (steps(collapse)-t0)./(steps(collapse)-t0+tau));
    boundaryTrace = [cfg.b0 boundaries];

    traceSelf = nan(nPaths,nNativeSteps,'single');
    traceOther = nan(nPaths,nNativeSteps,'single');
    traceSelf(:,1) = single(xSelf);
    traceOther(:,1) = single(xOther);
    choice = zeros(nPaths,1,'uint8');
    decisionIndex = nan(nPaths,1);
    active = true(nPaths,1);
    rng(seed,'twister');
    for step = 1:cfg.maxSteps
        if ~any(active)
            break;
        end
        rows = find(active);
        xs = xSelf(rows);
        xo = xOther(rows);
        nextSelf = max(0,xs+vSelf(rows)-cfg.delta*xs- ...
            lambda*xo+sigma*randn(numel(rows),1));
        nextOther = max(0,xo+vOther(rows)-cfg.delta*xo- ...
            lambda*xs+sigma*randn(numel(rows),1));
        xSelf(rows) = nextSelf;
        xOther(rows) = nextOther;
        traceSelf(rows,step+1) = single(nextSelf);
        traceOther(rows,step+1) = single(nextOther);
        hitSelf = nextSelf >= boundaries(step);
        hitOther = nextOther >= boundaries(step);
        localChoice = zeros(numel(rows),1,'uint8');
        if strcmp(cfg.tieRule,'original-self-priority')
            localChoice(hitSelf) = 1;
            localChoice(hitOther & ~hitSelf) = 2;
        else
            localChoice(hitSelf & ~hitOther) = 1;
            localChoice(hitOther & ~hitSelf) = 2;
            both = hitSelf & hitOther;
            localChoice(both & nextSelf >= nextOther) = 1;
            localChoice(both & nextOther > nextSelf) = 2;
        end
        crossed = localChoice > 0;
        globalRows = rows(crossed);
        choice(globalRows) = localChoice(crossed);
        decisionIndex(globalRows) = step+1;
        active(globalRows) = false;
    end

    valid = choice > 0 & isfinite(decisionIndex);
    summary.requestedCount = nRep*ones(nTrials,1);
    summary.validCount = accumarray( ...
        trialId(valid),1,[nTrials 1],@sum,0);
    summary.branchCount = accumarray( ...
        [trialId(valid),double(choice(valid))],1,[nTrials 2],@sum,0);
    fields = {'normSumB0Self','normSumB0Other', ...
        'normSumBoundarySelf','normSumBoundaryOther'};
    for f = 1:numel(fields)
        summary.(fields{f}) = zeros(nGroups,cfg.nGrid);
    end
    fields = {'endSumB0Self','endSumB0Other', ...
        'endSumBoundarySelf','endSumBoundaryOther','endPathCount'};
    for f = 1:numel(fields)
        summary.(fields{f}) = zeros(nGroups,nNativeSteps);
    end

    target = linspace(0,1,cfg.nGrid);
    durations = unique(decisionIndex(valid))';
    for durationValue = durations
        last = round(durationValue);
        rows = find(valid & decisionIndex == durationValue);
        groupId = (trialId(rows)-1)*2+double(choice(rows));
        groupMatrix = sparse(groupId,1:numel(rows),1, ...
            nGroups,numel(rows));
        position = 1+target*(last-1);
        lo = floor(position);
        hi = ceil(position);
        fraction = position-lo;
        selfInterp = double(traceSelf(rows,lo)).*(1-fraction)+ ...
            double(traceSelf(rows,hi)).*fraction;
        otherInterp = double(traceOther(rows,lo)).*(1-fraction)+ ...
            double(traceOther(rows,hi)).*fraction;
        boundaryInterp = boundaryTrace(lo).*(1-fraction)+ ...
            boundaryTrace(hi).*fraction;
        summary.normSumB0Self = summary.normSumB0Self+ ...
            groupMatrix*(selfInterp/cfg.b0);
        summary.normSumB0Other = summary.normSumB0Other+ ...
            groupMatrix*(otherInterp/cfg.b0);
        summary.normSumBoundarySelf = ...
            summary.normSumBoundarySelf+ ...
            groupMatrix*(selfInterp./boundaryInterp);
        summary.normSumBoundaryOther = ...
            summary.normSumBoundaryOther+ ...
            groupMatrix*(otherInterp./boundaryInterp);

        source = 1:last;
        destination = (nNativeSteps-last+1):nNativeSteps;
        selfSegment = double(traceSelf(rows,source));
        otherSegment = double(traceOther(rows,source));
        summary.endSumB0Self(:,destination) = ...
            summary.endSumB0Self(:,destination)+ ...
            groupMatrix*(selfSegment/cfg.b0);
        summary.endSumB0Other(:,destination) = ...
            summary.endSumB0Other(:,destination)+ ...
            groupMatrix*(otherSegment/cfg.b0);
        summary.endSumBoundarySelf(:,destination) = ...
            summary.endSumBoundarySelf(:,destination)+ ...
            groupMatrix*(selfSegment./boundaryTrace(source));
        summary.endSumBoundaryOther(:,destination) = ...
            summary.endSumBoundaryOther(:,destination)+ ...
            groupMatrix*(otherSegment./boundaryTrace(source));
        groupCount = full(sum(groupMatrix,2));
        summary.endPathCount(:,destination) = ...
            summary.endPathCount(:,destination)+groupCount;
    end
end

function total = add_batch_summary_local(total,batch)
    if isempty(total)
        total = batch;
        return;
    end
    names = fieldnames(batch);
    for f = 1:numel(names)
        total.(names{f}) = total.(names{f})+batch.(names{f});
    end
end

function result = finalize_summary_local(summary)
    nTrials = size(summary.validCount,1);
    inverseValid = zeros(nTrials,1);
    validTrial = summary.validCount > 0;
    inverseValid(validTrial) = 1./summary.validCount(validTrial);
    branchNames = {'stay','change'};
    for branch = 1:2
        name = branchNames{branch};
        rows = branch:2:(2*nTrials);
        branchCount = summary.branchCount(:,branch);
        branchWeight = branchCount.*inverseValid;
        denominator101 = sum(branchWeight);
        result.process101.(name).B0.self = scalar_divide_local( ...
            sum(summary.normSumB0Self(rows,:).*inverseValid,1), ...
            denominator101);
        result.process101.(name).B0.other = scalar_divide_local( ...
            sum(summary.normSumB0Other(rows,:).*inverseValid,1), ...
            denominator101);
        result.process101.(name).boundary.self = scalar_divide_local( ...
            sum(summary.normSumBoundarySelf(rows,:).*inverseValid,1), ...
            denominator101);
        result.process101.(name).boundary.other = scalar_divide_local( ...
            sum(summary.normSumBoundaryOther(rows,:).*inverseValid,1), ...
            denominator101);
        result.process101.(name).weight = denominator101;
        result.process101.(name).probability = ...
            denominator101/sum(validTrial);
        result.process101.(name).nPaths = sum(branchCount);

        pathCount = summary.endPathCount(rows,:);
        denominatorEnd = sum(pathCount.*inverseValid,1);
        result.end.(name).B0.self = vector_divide_local( ...
            sum(summary.endSumB0Self(rows,:).*inverseValid,1), ...
            denominatorEnd);
        result.end.(name).B0.other = vector_divide_local( ...
            sum(summary.endSumB0Other(rows,:).*inverseValid,1), ...
            denominatorEnd);
        result.end.(name).boundary.self = vector_divide_local( ...
            sum(summary.endSumBoundarySelf(rows,:).*inverseValid,1), ...
            denominatorEnd);
        result.end.(name).boundary.other = vector_divide_local( ...
            sum(summary.endSumBoundaryOther(rows,:).*inverseValid,1), ...
            denominatorEnd);
        result.end.(name).weight = denominatorEnd;
        result.end.(name).nPaths = sum(pathCount,1);
        result.end.(name).nTrials = sum(pathCount > 0,1);
    end
    stayRows = 1:2:(2*nTrials);
    changeRows = 2:2:(2*nTrials);
    result.end.nTrials = sum( ...
        summary.endPathCount(stayRows,:)+ ...
        summary.endPathCount(changeRows,:) > 0,1);
end

function output = assign_curves_local(output,pooled,draws,s,c)
    branches = {'stay','change'};
    normalizations = {'B0','boundary'};
    accumulators = {'self','other'};
    nDraws = numel(draws);
    for b = 1:numel(branches)
        branch = branches{b};
        for n = 1:numel(normalizations)
            normalization = normalizations{n};
            for a = 1:numel(accumulators)
                accumulator = accumulators{a};
                mean101 = pooled.process101.(branch). ...
                    (normalization).(accumulator);
                meanEnd = pooled.end.(branch). ...
                    (normalization).(accumulator);
                draw101 = nan(nDraws,numel(mean101));
                drawEnd = nan(nDraws,numel(meanEnd));
                for q = 1:nDraws
                    draw101(q,:) = draws{q}.process101.(branch). ...
                        (normalization).(accumulator);
                    drawEnd(q,:) = draws{q}.end.(branch). ...
                        (normalization).(accumulator);
                end
                output = store_distribution_local(output, ...
                    {'process101',branch,normalization,accumulator}, ...
                    mean101,draw101,s,c);
                output = store_distribution_local(output, ...
                    {'end',branch,normalization,accumulator}, ...
                    meanEnd,drawEnd,s,c);
            end
        end
    end

    for n = 1:numel(normalizations)
        normalization = normalizations{n};
        pooledWinner101 = weighted_merge_local( ...
            pooled.process101.stay.(normalization).self, ...
            pooled.process101.stay.weight, ...
            pooled.process101.change.(normalization).other, ...
            pooled.process101.change.weight);
        pooledLoser101 = weighted_merge_local( ...
            pooled.process101.stay.(normalization).other, ...
            pooled.process101.stay.weight, ...
            pooled.process101.change.(normalization).self, ...
            pooled.process101.change.weight);
        pooledWinnerEnd = weighted_merge_local( ...
            pooled.end.stay.(normalization).self, ...
            pooled.end.stay.weight, ...
            pooled.end.change.(normalization).other, ...
            pooled.end.change.weight);
        pooledLoserEnd = weighted_merge_local( ...
            pooled.end.stay.(normalization).other, ...
            pooled.end.stay.weight, ...
            pooled.end.change.(normalization).self, ...
            pooled.end.change.weight);
        drawWinner101 = nan(nDraws,numel(pooledWinner101));
        drawLoser101 = nan(size(drawWinner101));
        drawWinnerEnd = nan(nDraws,numel(pooledWinnerEnd));
        drawLoserEnd = nan(size(drawWinnerEnd));
        for q = 1:nDraws
            drawWinner101(q,:) = weighted_merge_local( ...
                draws{q}.process101.stay.(normalization).self, ...
                draws{q}.process101.stay.weight, ...
                draws{q}.process101.change.(normalization).other, ...
                draws{q}.process101.change.weight);
            drawLoser101(q,:) = weighted_merge_local( ...
                draws{q}.process101.stay.(normalization).other, ...
                draws{q}.process101.stay.weight, ...
                draws{q}.process101.change.(normalization).self, ...
                draws{q}.process101.change.weight);
            drawWinnerEnd(q,:) = weighted_merge_local( ...
                draws{q}.end.stay.(normalization).self, ...
                draws{q}.end.stay.weight, ...
                draws{q}.end.change.(normalization).other, ...
                draws{q}.end.change.weight);
            drawLoserEnd(q,:) = weighted_merge_local( ...
                draws{q}.end.stay.(normalization).other, ...
                draws{q}.end.stay.weight, ...
                draws{q}.end.change.(normalization).self, ...
                draws{q}.end.change.weight);
        end
        output = store_distribution_local(output, ...
            {'process101','winner',normalization}, ...
            pooledWinner101,drawWinner101,s,c);
        output = store_distribution_local(output, ...
            {'process101','loser',normalization}, ...
            pooledLoser101,drawLoser101,s,c);
        output = store_distribution_local(output, ...
            {'end','winner',normalization}, ...
            pooledWinnerEnd,drawWinnerEnd,s,c);
        output = store_distribution_local(output, ...
            {'end','loser',normalization}, ...
            pooledLoserEnd,drawLoserEnd,s,c);
    end

    output.process101.stay.weight(s,c) = ...
        single(pooled.process101.stay.weight);
    output.process101.change.weight(s,c) = ...
        single(pooled.process101.change.weight);
    output.process101.stayProbability(s,c) = ...
        single(pooled.process101.stay.probability);
    output.process101.changeProbability(s,c) = ...
        single(pooled.process101.change.probability);
    output.process101.nPaths(s,c) = uint32( ...
        pooled.process101.stay.nPaths+ ...
        pooled.process101.change.nPaths);
    output.end.stay.weight(s,c,:) = single(reshape( ...
        pooled.end.stay.weight,1,1,[]));
    output.end.change.weight(s,c,:) = single(reshape( ...
        pooled.end.change.weight,1,1,[]));
    output.end.stay.nTrials(s,c,:) = uint16(reshape( ...
        pooled.end.stay.nTrials,1,1,[]));
    output.end.change.nTrials(s,c,:) = uint16(reshape( ...
        pooled.end.change.nTrials,1,1,[]));
    output.end.nTrials(s,c,:) = uint16(reshape( ...
        pooled.end.nTrials,1,1,[]));
    output.end.nPaths(s,c,:) = uint32(reshape( ...
        pooled.end.stay.nPaths+pooled.end.change.nPaths,1,1,[]));
end

function output = store_distribution_local(output,path,meanCurve, ...
        drawCurves,s,c)
    interval = pointwise_interval_local(drawCurves);
    target = output;
    for p = 1:numel(path)
        target = target.(path{p});
    end
    target.mean(s,c,:) = single(reshape(meanCurve,1,1,[]));
    target.lower95(s,c,:) = single(reshape(interval(1,:),1,1,[]));
    target.upper95(s,c,:) = single(reshape(interval(2,:),1,1,[]));
    if numel(path) == 3
        output.(path{1}).(path{2}).(path{3}) = target;
    else
        output.(path{1}).(path{2}).(path{3}).(path{4}) = target;
    end
end

function interval = pointwise_interval_local(draws)
    interval = nan(2,size(draws,2));
    for t = 1:size(draws,2)
        values = draws(isfinite(draws(:,t)),t);
        if ~isempty(values)
            interval(:,t) = prctile(values,[2.5 97.5]);
        end
    end
end

function merged = weighted_merge_local(valueA,weightA,valueB,weightB)
    valueA = double(valueA);
    valueB = double(valueB);
    weightA = double(weightA)+zeros(size(valueA));
    weightB = double(weightB)+zeros(size(valueB));
    validA = isfinite(valueA) & weightA > 0;
    validB = isfinite(valueB) & weightB > 0;
    weightA(~validA) = 0;
    weightB(~validB) = 0;
    valueA(~validA) = 0;
    valueB(~validB) = 0;
    denominator = weightA+weightB;
    merged = nan(size(valueA));
    valid = denominator > 0;
    merged(valid) = (weightA(valid).*valueA(valid)+ ...
        weightB(valid).*valueB(valid))./denominator(valid);
end

function value = scalar_divide_local(numerator,denominator)
    if denominator > 0
        value = numerator/denominator;
    else
        value = nan(size(numerator));
    end
end

function value = vector_divide_local(numerator,denominator)
    value = nan(size(numerator));
    valid = denominator > 0;
    value(valid) = numerator(valid)./denominator(valid);
end

function trials = extract_conflict_trials_local(data,isSocial)
    valid = data.conflict(:) == 1 & ...
        data.answer1(:) ~= -1 & data.answer2(:) ~= -1 & ...
        data.conf1(:) ~= -1 & data.conf2(:) ~= -1 & ...
        isfinite(data.conf1(:)) & isfinite(data.conf2(:));
    if isSocial
        other = data.o1conf1(:);
    else
        other = data.inforate(:);
    end
    valid = valid & isfinite(other);
    selfConfidence = double(data.conf1(valid));
    other = double(other(valid));
    trials.stiSelf = selfConfidence/16-1/32;
    if isSocial
        trials.stiOther = other/16-1/32;
    else
        trials.stiOther = other-0.5;
    end
end

function validate_rows_local(rows,nMaximum,name)
    assert(all(isfinite(rows) & rows == round(rows) & ...
        rows >= 1 & rows <= nMaximum), ...
        '%s contains an invalid row.',name);
    assert(numel(unique(rows)) == numel(rows), ...
        '%s contains duplicate rows.',name);
end
