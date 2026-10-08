function hierarchy = fit_hierarchical_eb_parameters_85( ...
        parameterFile85,behaviorFile85,behaviorFile19,cfg)
%FIT_HIERARCHICAL_EB_PARAMETERS_85 Two-stage hierarchical EB shrinkage.
% This function jointly models the 85 independently fitted parameter
% vectors on their bounded-logit scales. It is an empirical-Bayes
% approximation, not an exact posterior from the stochastic trial
% likelihood. Trial counts enter as precision weights.

    arguments
        parameterFile85 (1,:) char
        behaviorFile85 (1,:) char
        behaviorFile19 (1,:) char
        cfg.nPosteriorDraws (1,1) double = 200
        cfg.seed (1,1) double = 20260721
        cfg.boundaryClip (1,1) double = 1e-5
        cfg.correlationRegularization (1,1) double = 0.25
        cfg.minimumSEFraction (1,1) double = 0.02
        cfg.megExpectedRows (1,:) double = 38:56
        cfg.parameterVariant (1,:) char = 'minb1'
        cfg.parameterValueField (1,:) char = 'optParams'
        cfg.fixedLambda (1,1) logical = false
    end

    assert(cfg.nPosteriorDraws >= 20 && ...
        cfg.nPosteriorDraws == round(cfg.nPosteriorDraws), ...
        'cfg.nPosteriorDraws must be an integer of at least 20.');

    assert(ismember(cfg.parameterVariant,{'minb1','minb2'}), ...
        'cfg.parameterVariant must be minb1 or minb2.');
    fitNames = strcat({'ots_','pls_','otn_','pln_'}, ...
        cfg.parameterVariant);
    p = load(parameterFile85,fitNames{:});
    d85 = load(behaviorFile85, ...
        'ots_data','pls_data','otn_data','pln_data');
    d19 = load(behaviorFile19, ...
        'ots_data','pls_data','otn_data','pln_data');

    dataNames = {'ots_data','pls_data','otn_data','pln_data'};
    conditionNames = {'OT-social','PL-social', ...
        'OT-nonsocial','PL-nonsocial'};
    parameterNames = {'v0','b1_ratio','b2_ratio','k1','k2', ...
        'sigma','lambda','tau','t0'};
    lowerBound = [0.001 0 0 0.001 0.001 0.1 1e-5 10 10];
    upperBound = [0.99 1 1 10 10 20 0.15 1000 1000];
    drugCode = [1 0 1 0];
    contextCode = [1 1 0 0];

    fitSets = cellfun(@(x) p.(x),fitNames,'UniformOutput',false);
    dataSets85 = cellfun(@(x) d85.(x),dataNames,'UniformOutput',false);
    dataSets19 = cellfun(@(x) d19.(x),dataNames,'UniformOutput',false);
    nSubjects = numel(fitSets{1});
    nConditions = numel(fitSets);
    nParameters = numel(parameterNames);
    assert(nSubjects == 85, ...
        'The all-cohort parameter file must contain 85 subjects.');
    verify_set_sizes_local(fitSets,dataSets85,nSubjects);

    subjectIDs = extract_ids_local(dataSets85{1});
    megIDs = extract_ids_local(dataSets19{1});
    [isMatched,megRows] = ismember(megIDs,subjectIDs);
    assert(all(isMatched) && numel(unique(megRows)) == numel(megRows), ...
        'The 19 MEG subject IDs do not map uniquely into the 85 subjects.');
    assert(isequal(megRows(:)',cfg.megExpectedRows), ...
        'ID matching does not agree with cfg.megExpectedRows.');
    for c = 2:nConditions
        assert(isequal(extract_ids_local(dataSets85{c}),subjectIDs), ...
            'Subject order differs across the 85-subject conditions.');
        assert(isequal(extract_ids_local(dataSets19{c}),megIDs), ...
            'Subject order differs across the 19-subject conditions.');
    end

    cohort = [ones(37,1); 2*ones(19,1); 3*ones(29,1)];
    assert(numel(cohort) == nSubjects);
    observed = nan(nSubjects,nConditions,nParameters);
    trialCount = zeros(nSubjects,nConditions);
    for c = 1:nConditions
        parameterMatrix = nan(nSubjects,nParameters);
        for s = 1:nSubjects
            value = get_parameter_value_local( ...
                fitSets{c}(s),cfg.parameterValueField);
            assert(isnumeric(value) && numel(value) == nParameters, ...
                'Subject %d, condition %d does not have 9 parameters.', ...
                s,c);
            parameterMatrix(s,:) = reshape(double(value),1,[]);
            trialCount(s,c) = count_valid_conflict_trials_local( ...
                dataSets85{c}(s),c <= 2);
        end
        observed(:,c,:) = reshape( ...
            parameterMatrix,nSubjects,1,nParameters);
    end
    assert(all(isfinite(observed(:))), ...
        'At least one optParams vector contains a nonfinite value.');
    assert(all(trialCount(:) > 0), ...
        'At least one subject-condition has no valid conflict trials.');
    if cfg.fixedLambda
        assert(all(abs(observed(:,:,7)) < 1e-12,'all'), ...
            'cfg.fixedLambda is true, but parameter 7 (lambda) is not identically zero.');
    end

    latentObserved = bounded_logit_local( ...
        observed,lowerBound,upperBound,cfg.boundaryClip);
    rowSubject = repelem((1:nSubjects)',nConditions);
    rowCondition = repmat((1:nConditions)',nSubjects,1);
    rowDrug = drugCode(rowCondition)';
    rowContext = contextCode(rowCondition)';
    rowCohort = cohort(rowSubject);
    rowWeight = reshape(trialCount',[],1);
    rowWeight = rowWeight/mean(rowWeight);
    subjectLabels = categorical(subjectIDs(rowSubject));
    cohortLabels = categorical(rowCohort,1:3, ...
        {'cohort1','MEG','cohort3'});

    latentMean = nan(size(latentObserved));
    latentLower = nan(size(latentObserved));
    latentUpper = nan(size(latentObserved));
    latentSE = nan(size(latentObserved));
    parameterDraws = nan( ...
        nSubjects,nConditions,nParameters,cfg.nPosteriorDraws,'single');
    model = repmat(struct(),nParameters,1);

    rng(cfg.seed,'twister');
    formula = ['Value ~ Drug*Context + Cohort + ', ...
        '(1 + Drug | Subject)'];
    for parameter = 1:nParameters
        if cfg.fixedLambda && parameter == 7
            latentMean(:,:,parameter) = latentObserved(:,:,parameter);
            latentLower(:,:,parameter) = latentObserved(:,:,parameter);
            latentUpper(:,:,parameter) = latentObserved(:,:,parameter);
            latentSE(:,:,parameter) = 0;
            parameterDraws(:,:,parameter,:) = single(repmat( ...
                latentObserved(:,:,parameter),1,1,1, ...
                cfg.nPosteriorDraws));
            model(parameter).parameter = parameterNames{parameter};
            model(parameter).formula = 'fixed at zero; no LME fitted';
            model(parameter).fixedValue = 0;
            continue;
        end
        responseMatrix = latentObserved(:,:,parameter);
        response = reshape(responseMatrix',[],1);
        tbl = table(response,subjectLabels,rowDrug,rowContext, ...
            cohortLabels,rowWeight, ...
            'VariableNames',{'Value','Subject','Drug','Context', ...
            'Cohort','TrialWeight'});

        lme = fitlme(tbl,formula, ...
            'FitMethod','REML', ...
            'DummyVarCoding','effects', ...
            'CovariancePattern','Isotropic', ...
            'Weights',tbl.TrialWeight, ...
            'CheckHessian',true);
        [prediction,predictionCI] = predict(lme,tbl, ...
            'Conditional',true, ...
            'Prediction','curve', ...
            'DFMethod','none', ...
            'Simultaneous',false);
        predictionSE = (predictionCI(:,2)-predictionCI(:,1))/ ...
            (2*1.95996398454005);
        responseScale = std(response,0,'omitnan');
        predictionSE = max(predictionSE, ...
            cfg.minimumSEFraction*responseScale);

        meanMatrix = reshape(prediction,nConditions,nSubjects)';
        lowerMatrix = reshape(predictionCI(:,1), ...
            nConditions,nSubjects)';
        upperMatrix = reshape(predictionCI(:,2), ...
            nConditions,nSubjects)';
        seMatrix = reshape(predictionSE,nConditions,nSubjects)';
        latentMean(:,:,parameter) = meanMatrix;
        latentLower(:,:,parameter) = lowerMatrix;
        latentUpper(:,:,parameter) = upperMatrix;
        latentSE(:,:,parameter) = seMatrix;

        residualMatrix = responseMatrix-meanMatrix;
        residualCorrelation = corr(residualMatrix, ...
            'Rows','pairwise');
        residualCorrelation(~isfinite(residualCorrelation)) = 0;
        residualCorrelation(1:nConditions+1:end) = 1;
        residualCorrelation = (1-cfg.correlationRegularization)* ...
            residualCorrelation + ...
            cfg.correlationRegularization*eye(nConditions);
        residualCorrelation = nearest_correlation_local( ...
            residualCorrelation);

        for s = 1:nSubjects
            covariance = diag(seMatrix(s,:))* ...
                residualCorrelation*diag(seMatrix(s,:));
            factor = chol(nearest_spd_local(covariance),'lower');
            draws = meanMatrix(s,:)' + factor*randn( ...
                nConditions,cfg.nPosteriorDraws);
            parameterDraws(s,:,parameter,:) = single(reshape( ...
                draws,1,nConditions,1,cfg.nPosteriorDraws));
        end

        [psi,mse] = covarianceParameters(lme);
        model(parameter).parameter = parameterNames{parameter};
        model(parameter).formula = formula;
        model(parameter).coefficients = lme.Coefficients;
        model(parameter).coefficientCovariance = ...
            lme.CoefficientCovariance;
        model(parameter).randomCovariance = psi;
        model(parameter).residualVariance = mse;
        model(parameter).residualCorrelationAcrossConditions = ...
            residualCorrelation;
        model(parameter).modelCriterion = lme.ModelCriterion;
        model(parameter).rsquared = lme.Rsquared;
    end

    originalDraws = inverse_bounded_logit_local( ...
        double(parameterDraws),lowerBound,upperBound);
    originalMean = inverse_bounded_logit_local( ...
        latentMean,lowerBound,upperBound);
    originalLower = inverse_bounded_logit_local( ...
        latentLower,lowerBound,upperBound);
    originalUpper = inverse_bounded_logit_local( ...
        latentUpper,lowerBound,upperBound);
    if cfg.fixedLambda
        originalDraws(:,:,7,:) = 0;
        originalMean(:,:,7) = 0;
        originalLower(:,:,7) = 0;
        originalUpper(:,:,7) = 0;
    end

    hierarchy.meta.created = datetime('now');
    hierarchy.meta.method = ...
        'two-stage joint hierarchical empirical Bayes approximation';
    hierarchy.meta.warning = ['Not an exact posterior from the stochastic ', ...
        'trial likelihood; independently fitted optParams are modeled with ', ...
        'trial-count-weighted linear mixed models.'];
    hierarchy.meta.parameterFile85 = parameterFile85;
    hierarchy.meta.behaviorFile85 = behaviorFile85;
    hierarchy.meta.behaviorFile19 = behaviorFile19;
    hierarchy.meta.parameterNames = parameterNames;
    hierarchy.meta.parameterVariant = cfg.parameterVariant;
    hierarchy.meta.parameterValueField = cfg.parameterValueField;
    hierarchy.meta.fixedLambda = cfg.fixedLambda;
    hierarchy.meta.conditionNames = conditionNames;
    hierarchy.meta.conditionOrder = 1:nConditions;
    hierarchy.meta.drugCode = drugCode;
    hierarchy.meta.contextCode = contextCode;
    hierarchy.meta.parameterLowerBound = lowerBound;
    hierarchy.meta.parameterUpperBound = upperBound;
    hierarchy.meta.subjectIDs = subjectIDs;
    hierarchy.meta.cohort = cohort;
    hierarchy.meta.megIDs = megIDs;
    hierarchy.meta.megRowsIn85 = megRows;
    hierarchy.meta.dimensions = ...
        'subject x condition x parameter x posterior draw';
    hierarchy.config = cfg;
    hierarchy.trialCount = trialCount;
    hierarchy.unpooled.parameters = observed;
    hierarchy.latent.observed = latentObserved;
    hierarchy.latent.mean = latentMean;
    hierarchy.latent.lower95 = latentLower;
    hierarchy.latent.upper95 = latentUpper;
    hierarchy.latent.standardError = latentSE;
    hierarchy.parameters.mean = originalMean;
    hierarchy.parameters.lower95 = originalLower;
    hierarchy.parameters.upper95 = originalUpper;
    hierarchy.parameters.draws = single(originalDraws);
    hierarchy.meg.rowsIn85 = megRows;
    hierarchy.meg.subjectIDs = megIDs;
    hierarchy.meg.parameters.mean = originalMean(megRows,:,:);
    hierarchy.meg.parameters.lower95 = originalLower(megRows,:,:);
    hierarchy.meg.parameters.upper95 = originalUpper(megRows,:,:);
    hierarchy.meg.parameters.draws = ...
        single(originalDraws(megRows,:,:,:));
    hierarchy.model = model;
end

function verify_set_sizes_local(fitSets,dataSets,nSubjects)
    for c = 1:numel(fitSets)
        assert(numel(fitSets{c}) == nSubjects && ...
            numel(dataSets{c}) == nSubjects, ...
            'Condition %d does not contain %d subjects.',c,nSubjects);
    end
end

function value = get_parameter_value_local(record,preferredField)
    if isfield(record,preferredField)
        value = record.(preferredField);
        return;
    end
    % The compare pipeline has used both names for the restored 9-slot
    % vector. Accept the fallback, but keep the preferred field explicit in
    % the saved metadata.
    alternatives = {'optParams','optParams_FullArray'};
    for i = 1:numel(alternatives)
        if isfield(record,alternatives{i})
            value = record.(alternatives{i});
            return;
        end
    end
    error('No parameter vector field found in the fit record.');
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

function count = count_valid_conflict_trials_local(data,isSocial)
    valid = data.conflict(:) == 1 & ...
        data.answer1(:) ~= -1 & data.answer2(:) ~= -1 & ...
        data.conf1(:) ~= -1 & data.conf2(:) ~= -1 & ...
        isfinite(data.conf1(:)) & isfinite(data.conf2(:));
    if isSocial
        valid = valid & isfinite(data.o1conf1(:));
    else
        valid = valid & isfinite(data.inforate(:));
    end
    count = sum(valid);
end

function latent = bounded_logit_local(value,lower,upper,clipValue)
    lower = reshape(lower,1,1,[]);
    upper = reshape(upper,1,1,[]);
    proportion = (value-lower)./(upper-lower);
    proportion = min(max(proportion,clipValue),1-clipValue);
    latent = log(proportion./(1-proportion));
end

function value = inverse_bounded_logit_local(latent,lower,upper)
    shape = ones(1,ndims(latent));
    shape(3) = numel(lower);
    lower = reshape(lower,shape);
    upper = reshape(upper,shape);
    probability = 1./(1+exp(-latent));
    value = lower+(upper-lower).*probability;
end

function correlation = nearest_correlation_local(correlation)
    correlation = (correlation+correlation')/2;
    [vectors,values] = eig(correlation);
    values = diag(max(diag(values),1e-6));
    correlation = vectors*values*vectors';
    scale = sqrt(diag(correlation));
    correlation = correlation./(scale*scale');
    correlation = (correlation+correlation')/2;
end

function matrix = nearest_spd_local(matrix)
    matrix = (matrix+matrix')/2;
    [vectors,values] = eig(matrix);
    floorValue = max(max(diag(values))*1e-10,1e-10);
    values = diag(max(diag(values),floorValue));
    matrix = vectors*values*vectors';
    matrix = (matrix+matrix')/2;
end
