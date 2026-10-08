function trialBank = build_hierarchical_meg_trial_trajectory_bank( ...
        hierarchy,behaviorFile19,options)
%BUILD_HIERARCHICAL_MEG_TRIAL_TRAJECTORY_BANK Trial-level path summaries.
% The bank is shared by permutation and balanced-trial mixed controls.
% Parameter uncertainty is integrated by pooling paths generated from
% multiple hierarchical EB parameter draws.

    arguments
        hierarchy (1,1) struct
        behaviorFile19 (1,:) char
        options.nParameterDraws (1,1) double = 30
        options.nRepPerDraw (1,1) double = 100
        options.nGrid (1,1) double = 101
        options.maxSteps (1,1) double = 1000
        options.b0 (1,1) double = 300
        options.delta (1,1) double = 0.01
        options.seed (1,1) double = 20260723
        options.megSubjectRows (1,:) double = []
        options.conditionRows (1,:) double = []
        options.checkpointFile (1,:) char = ''
        options.useParallel (1,1) logical = false
        options.nWorkers (1,1) double = 0
    end

    parameterDraws = hierarchy.meg.parameters.draws;
    nAvailableSubjects = size(parameterDraws,1);
    nAvailableConditions = size(parameterDraws,2);
    nAvailableDraws = size(parameterDraws,4);
    if isempty(options.megSubjectRows)
        options.megSubjectRows = 1:nAvailableSubjects;
    end
    if isempty(options.conditionRows)
        options.conditionRows = 1:nAvailableConditions;
    end
    validate_rows_local(options.megSubjectRows,nAvailableSubjects, ...
        'options.megSubjectRows');
    validate_rows_local(options.conditionRows,nAvailableConditions, ...
        'options.conditionRows');
    assert(options.nParameterDraws >= 1 && ...
        options.nParameterDraws <= nAvailableDraws, ...
        'nParameterDraws must be between 1 and %d.',nAvailableDraws);
    drawRows = unique(round(linspace(1,nAvailableDraws, ...
        options.nParameterDraws)),'stable');
    assert(numel(drawRows) == options.nParameterDraws);

    d = load(behaviorFile19, ...
        'ots_data','pls_data','otn_data','pln_data');
    dataSets = {d.ots_data,d.pls_data,d.otn_data,d.pln_data};
    isSocial = [true true false false];
    conditionNames = {'OT-social','PL-social', ...
        'OT-nonsocial','PL-nonsocial'};
    subjectRows = options.megSubjectRows;
    conditionRows = options.conditionRows;
    nSubjects = numel(subjectRows);
    nConditions = numel(conditionRows);

    trialBank.meta.created = datetime('now');
    trialBank.meta.type = 'hierarchical-MEG19-trial-trajectory-bank';
    trialBank.meta.parameterNames = hierarchy.meta.parameterNames;
    trialBank.meta.behaviorFile19 = behaviorFile19;
    trialBank.meta.subjectRowsInMEG19 = subjectRows;
    trialBank.meta.subjectRowsIn85 = ...
        hierarchy.meg.rowsIn85(subjectRows);
    trialBank.meta.subjectIDs = hierarchy.meg.subjectIDs(subjectRows);
    trialBank.meta.conditionRows = conditionRows;
    trialBank.meta.conditionNames = conditionNames(conditionRows);
    trialBank.meta.parameterDrawRows = drawRows;
    trialBank.meta.pathsPerTrial = ...
        options.nParameterDraws*options.nRepPerDraw;
    trialBank.meta.process101 = linspace(0,1,options.nGrid);
    trialBank.meta.endStep = -options.maxSteps:0;
    trialBank.config = options;
    trialBank.subject = cell(nSubjects,nConditions);
    trialBank.completed = false(nSubjects,nConditions);

    nJobs = nSubjects*nConditions;
    jobSubject = repmat((1:nSubjects)',nConditions,1);
    jobCondition = repelem((1:nConditions)',nSubjects);
    jobSourceSubject = subjectRows(jobSubject(:));
    jobSourceCondition = conditionRows(jobCondition(:));
    jobIsSocial = isSocial(jobSourceCondition(:));
    jobConditionName = conditionNames(jobSourceCondition(:));
    jobData = cell(nJobs,1);
    for job = 1:nJobs
        sourceSubject = jobSourceSubject(job);
        sourceCondition = jobSourceCondition(job);
        jobData{job} = dataSets{sourceCondition}(sourceSubject);
    end

    runParallel = options.useParallel && ...
        license('test','Distrib_Computing_Toolbox') && ...
        exist('parpool','file') == 2;
    if options.useParallel && ~runParallel
        warning('Parallel toolbox unavailable; using serial execution.');
    end
    if runParallel
        start_parallel_pool_local(options.nWorkers);
        jobResult = cell(nJobs,1);
        parfor job = 1:nJobs
            s = jobSubject(job);
            c = jobCondition(job);
            sourceSubject = jobSourceSubject(job);
            sourceCondition = jobSourceCondition(job);
            jobResult{job} = build_one_bank_local( ...
                parameterDraws,drawRows,sourceSubject,sourceCondition, ...
                jobData{job},jobIsSocial(job),options, ...
                jobConditionName{job},s,nSubjects);
        end
        for job = 1:nJobs
            s = jobSubject(job);
            c = jobCondition(job);
            trialBank.subject{s,c} = jobResult{job};
            trialBank.completed(s,c) = true;
        end
        if ~isempty(options.checkpointFile)
            save(options.checkpointFile,'trialBank','-v7.3');
        end
    else
        for job = 1:nJobs
            s = jobSubject(job);
            c = jobCondition(job);
            sourceSubject = jobSourceSubject(job);
            sourceCondition = jobSourceCondition(job);
            trialBank.subject{s,c} = build_one_bank_local( ...
                parameterDraws,drawRows,sourceSubject,sourceCondition, ...
                jobData{job},jobIsSocial(job),options, ...
                jobConditionName{job},s,nSubjects);
            trialBank.completed(s,c) = true;
            if ~isempty(options.checkpointFile)
                save(options.checkpointFile,'trialBank','-v7.3');
            end
        end
    end
    trialBank.meta.completed = datetime('now');
end

function bank = build_one_bank_local(parameterDraws,drawRows, ...
        sourceSubject,sourceCondition,data,isSocial,options, ...
        conditionName,displaySubject,nDisplaySubjects)
    fprintf('Trial bank: %s, subject %d/%d\n', ...
        conditionName,displaySubject,nDisplaySubjects);
    trials = extract_trials_local(data,isSocial);
    assert(~isempty(trials.stiSelf), ...
        'No valid conflict trials for subject %d, condition %d.', ...
        sourceSubject,sourceCondition);
    pooled = [];
    for drawIndex = 1:options.nParameterDraws
        sourceDraw = drawRows(drawIndex);
        params = double(reshape(parameterDraws( ...
            sourceSubject,sourceCondition,:,sourceDraw),1,[]));
        seed = options.seed+10000000*sourceCondition+ ...
            10000*sourceSubject+drawIndex;
        summary = simulate_trial_summary_local( ...
            params,trials,options,options.nRepPerDraw,seed);
        pooled = add_summary_local(pooled,summary);
    end
    bank = finalize_trial_bank_local(pooled,trials,options);
end

function start_parallel_pool_local(nWorkers)
    assert(nWorkers >= 0 && nWorkers == round(nWorkers), ...
        'options.nWorkers must be zero or a positive integer.');
    pool = gcp('nocreate');
    if isempty(pool)
        if nWorkers > 0
            parpool('Processes',nWorkers);
        else
            parpool('Processes');
        end
    elseif nWorkers > 0 && pool.NumWorkers ~= nWorkers
        warning(['Using the existing pool with %d workers; requested %d. ', ...
            'Delete the existing pool first to change its size.'], ...
            pool.NumWorkers,nWorkers);
    end
end

function summary = simulate_trial_summary_local(params,trials,options,nRep,seed)
    nTrials = numel(trials.stiSelf);
    nPaths = nTrials*nRep;
    nNativeSteps = options.maxSteps+1;
    nGroups = 2*nTrials;
    trialId = repelem((1:nTrials)',nRep);
    stiSelf = repelem(trials.stiSelf(:),nRep);
    stiOther = repelem(trials.stiOther(:),nRep);

    v0 = params(1);
    xSelf = options.b0*max(0,params(2))*ones(nPaths,1);
    xOther = options.b0*max(0,params(3))*ones(nPaths,1);
    assert(numel(params) == 9 && all(isfinite(params)), ...
        'The trajectory bank requires 9 finite parameters.');
    kSelf = params(4);
    kOther = params(5);
    sigma = params(6);
    lambda = params(7);
    tau = params(8);
    t0 = params(9);
    vSelf = v0+kSelf*stiSelf;
    vOther = v0+kOther*stiOther;

    steps = 1:options.maxSteps;
    boundaries = options.b0*ones(1,options.maxSteps);
    collapse = steps >= t0;
    boundaries(collapse) = options.b0*(1-0.9* ...
        (steps(collapse)-t0)./(steps(collapse)-t0+tau));
    boundaryTrace = [options.b0 boundaries];

    traceSelf = nan(nPaths,nNativeSteps,'single');
    traceOther = nan(nPaths,nNativeSteps,'single');
    traceSelf(:,1) = single(xSelf);
    traceOther(:,1) = single(xOther);
    choice = zeros(nPaths,1,'uint8');
    decisionIndex = nan(nPaths,1);
    active = true(nPaths,1);
    rng(seed,'twister');
    for step = 1:options.maxSteps
        if ~any(active)
            break;
        end
        rows = find(active);
        xs = xSelf(rows);
        xo = xOther(rows);
        nextSelf = max(0,xs+vSelf(rows)-options.delta*xs- ...
            lambda*xo+sigma*randn(numel(rows),1));
        nextOther = max(0,xo+vOther(rows)-options.delta*xo- ...
            lambda*xs+sigma*randn(numel(rows),1));
        xSelf(rows) = nextSelf;
        xOther(rows) = nextOther;
        traceSelf(rows,step+1) = single(nextSelf);
        traceOther(rows,step+1) = single(nextOther);
        hitSelf = nextSelf >= boundaries(step);
        hitOther = nextOther >= boundaries(step) & ~hitSelf;
        localChoice = zeros(numel(rows),1,'uint8');
        localChoice(hitSelf) = 1;
        localChoice(hitOther) = 2;
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
    names = {'normSumB0Self','normSumB0Other', ...
        'normSumBoundarySelf','normSumBoundaryOther', ...
        'normSumB0SelfBC','normSumB0OtherBC', ...
        'normSumBoundarySelfBC','normSumBoundaryOtherBC'};
    for f = 1:numel(names)
        summary.(names{f}) = zeros(nGroups,options.nGrid);
    end
    names = {'endSumB0Self','endSumB0Other', ...
        'endSumBoundarySelf','endSumBoundaryOther', ...
        'endSumB0SelfBC','endSumB0OtherBC', ...
        'endSumBoundarySelfBC','endSumBoundaryOtherBC', ...
        'endPathCount'};
    for f = 1:numel(names)
        summary.(names{f}) = zeros(nGroups,nNativeSteps);
    end

    target = linspace(0,1,options.nGrid);
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
        selfStartB0 = double(traceSelf(rows,1))/options.b0;
        otherStartB0 = double(traceOther(rows,1))/options.b0;
        summary.normSumB0Self = summary.normSumB0Self+ ...
            groupMatrix*(selfInterp/options.b0);
        summary.normSumB0Other = summary.normSumB0Other+ ...
            groupMatrix*(otherInterp/options.b0);
        summary.normSumBoundarySelf = ...
            summary.normSumBoundarySelf+ ...
            groupMatrix*(selfInterp./boundaryInterp);
        summary.normSumBoundaryOther = ...
            summary.normSumBoundaryOther+ ...
            groupMatrix*(otherInterp./boundaryInterp);
        summary.normSumB0SelfBC = summary.normSumB0SelfBC+ ...
            groupMatrix*(selfInterp/options.b0-selfStartB0);
        summary.normSumB0OtherBC = summary.normSumB0OtherBC+ ...
            groupMatrix*(otherInterp/options.b0-otherStartB0);
        summary.normSumBoundarySelfBC = ...
            summary.normSumBoundarySelfBC+ ...
            groupMatrix*(selfInterp./boundaryInterp-selfStartB0);
        summary.normSumBoundaryOtherBC = ...
            summary.normSumBoundaryOtherBC+ ...
            groupMatrix*(otherInterp./boundaryInterp-otherStartB0);

        source = 1:last;
        destination = (nNativeSteps-last+1):nNativeSteps;
        selfSegment = double(traceSelf(rows,source));
        otherSegment = double(traceOther(rows,source));
        summary.endSumB0Self(:,destination) = ...
            summary.endSumB0Self(:,destination)+ ...
            groupMatrix*(selfSegment/options.b0);
        summary.endSumB0Other(:,destination) = ...
            summary.endSumB0Other(:,destination)+ ...
            groupMatrix*(otherSegment/options.b0);
        summary.endSumBoundarySelf(:,destination) = ...
            summary.endSumBoundarySelf(:,destination)+ ...
            groupMatrix*(selfSegment./boundaryTrace(source));
        summary.endSumBoundaryOther(:,destination) = ...
            summary.endSumBoundaryOther(:,destination)+ ...
            groupMatrix*(otherSegment./boundaryTrace(source));
        summary.endSumB0SelfBC(:,destination) = ...
            summary.endSumB0SelfBC(:,destination)+ ...
            groupMatrix*(selfSegment/options.b0-selfStartB0);
        summary.endSumB0OtherBC(:,destination) = ...
            summary.endSumB0OtherBC(:,destination)+ ...
            groupMatrix*(otherSegment/options.b0-otherStartB0);
        summary.endSumBoundarySelfBC(:,destination) = ...
            summary.endSumBoundarySelfBC(:,destination)+ ...
            groupMatrix*(selfSegment./boundaryTrace(source)-selfStartB0);
        summary.endSumBoundaryOtherBC(:,destination) = ...
            summary.endSumBoundaryOtherBC(:,destination)+ ...
            groupMatrix*(otherSegment./boundaryTrace(source)-otherStartB0);
        count = full(sum(groupMatrix,2));
        summary.endPathCount(:,destination) = ...
            summary.endPathCount(:,destination)+count;
    end
end

function total = add_summary_local(total,current)
    if isempty(total)
        total = current;
        return;
    end
    names = fieldnames(current);
    for f = 1:numel(names)
        total.(names{f}) = total.(names{f})+current.(names{f});
    end
end

function bank = finalize_trial_bank_local(summary,trials,options)
    nTrials = numel(trials.stiSelf);
    branches = {'stay','change'};
    bank.observedChoice = trials.observedChoice;
    bank.sourceTrialIndex = trials.sourceTrialIndex;
    bank.confSelf1 = trials.confSelf1;
    bank.confOther1 = trials.confOther1;
    bank.confSelf2 = trials.confSelf2;
    bank.stiSelf = trials.stiSelf;
    bank.stiOther = trials.stiOther;
    bank.requestedPaths = summary.requestedCount;
    bank.validPaths = summary.validCount;
    bank.process101.progress = linspace(0,1,options.nGrid);
    bank.end.step = -options.maxSteps:0;
    for branch = 1:2
        name = branches{branch};
        rows = branch:2:(2*nTrials);
        counts = summary.branchCount(:,branch);
        bank.process101.(name).nPaths = counts;
        bank.process101.(name).probability = ...
            safe_divide_local(counts,summary.validCount);
        bank.process101.(name).B0.self = row_divide_local( ...
            summary.normSumB0Self(rows,:),counts);
        bank.process101.(name).B0.other = row_divide_local( ...
            summary.normSumB0Other(rows,:),counts);
        bank.process101.(name).boundary.self = row_divide_local( ...
            summary.normSumBoundarySelf(rows,:),counts);
        bank.process101.(name).boundary.other = row_divide_local( ...
            summary.normSumBoundaryOther(rows,:),counts);
        bank.process101.(name).baselineCorrected.B0.self = ...
            row_divide_local(summary.normSumB0SelfBC(rows,:),counts);
        bank.process101.(name).baselineCorrected.B0.other = ...
            row_divide_local(summary.normSumB0OtherBC(rows,:),counts);
        bank.process101.(name).baselineCorrected.boundary.self = ...
            row_divide_local( ...
            summary.normSumBoundarySelfBC(rows,:),counts);
        bank.process101.(name).baselineCorrected.boundary.other = ...
            row_divide_local( ...
            summary.normSumBoundaryOtherBC(rows,:),counts);

        endCounts = summary.endPathCount(rows,:);
        bank.end.(name).nPaths = endCounts;
        bank.end.(name).B0.self = row_divide_local( ...
            summary.endSumB0Self(rows,:),endCounts);
        bank.end.(name).B0.other = row_divide_local( ...
            summary.endSumB0Other(rows,:),endCounts);
        bank.end.(name).boundary.self = row_divide_local( ...
            summary.endSumBoundarySelf(rows,:),endCounts);
        bank.end.(name).boundary.other = row_divide_local( ...
            summary.endSumBoundaryOther(rows,:),endCounts);
        bank.end.(name).baselineCorrected.B0.self = ...
            row_divide_local(summary.endSumB0SelfBC(rows,:),endCounts);
        bank.end.(name).baselineCorrected.B0.other = ...
            row_divide_local(summary.endSumB0OtherBC(rows,:),endCounts);
        bank.end.(name).baselineCorrected.boundary.self = ...
            row_divide_local( ...
            summary.endSumBoundarySelfBC(rows,:),endCounts);
        bank.end.(name).baselineCorrected.boundary.other = ...
            row_divide_local( ...
            summary.endSumBoundaryOtherBC(rows,:),endCounts);
    end
end

function output = row_divide_local(numerator,denominator)
    if isvector(denominator)
        denominator = denominator(:);
    end
    output = nan(size(numerator),'single');
    expanded = denominator+zeros(size(numerator));
    valid = expanded > 0;
    output(valid) = single(numerator(valid)./expanded(valid));
end

function output = safe_divide_local(numerator,denominator)
    output = nan(size(numerator),'single');
    valid = denominator > 0;
    output(valid) = single(numerator(valid)./denominator(valid));
end

function trials = extract_trials_local(data,isSocial)
    if isSocial
        other = data.o1conf1(:);
    else
        other = data.inforate(:);
    end
    valid = data.conflict(:) == 1 & ...
        data.answer1(:) ~= -1 & data.answer2(:) ~= -1 & ...
        data.conf1(:) ~= -1 & data.conf2(:) ~= -1 & ...
        isfinite(data.conf1(:)) & isfinite(data.conf2(:)) & ...
        isfinite(other) & isfinite(data.ischange(:));
    trials.sourceTrialIndex = find(valid);
    trials.observedChoice = uint8(double(data.ischange(valid))+1);
    trials.confSelf1 = double(data.conf1(valid));
    trials.confOther1 = double(other(valid));
    trials.confSelf2 = double(data.conf2(valid));
    trials.stiSelf = trials.confSelf1/16-1/32;
    other = double(other(valid));
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
