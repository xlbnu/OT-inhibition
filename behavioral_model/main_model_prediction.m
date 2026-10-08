%% Parameter settings
codeRoot = fileparts(fileparts(mfilename('fullpath')));
n_run = 500;                       % Number of prediction repetitions
lapse_rate = 0.05;                  % Probability of random choice and confidence
cohortID = [ones(37, 1); 2 * ones(19, 1); 3 * ones(29, 1)]; % 85 participants' cohort
output_dir = fullfile(codeRoot,'analysis_outputs','model_prediction_results');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Single model prediction: Preserve the original code's use of ot_valid to select behavioral data
ots_preds = generate_model_predictions(ots_minb1, ots_data(ot_valid));
pls_preds = generate_model_predictions(pls_minb1, pls_data(ot_valid));
otn_preds = generate_model_predictions(otn_minb1, otn_data(ot_valid));
pln_preds = generate_model_predictions(pln_minb1, pln_data(ot_valid));
save(fullfile(output_dir, 'prediction_data_85.mat')); % Save the current workspace as in the original code

%% Repeated prediction and mixed-effects models: Preserve the original code's use of complete behavioral data
beta_runs = cell(n_run, 1);
parfor j = 1:n_run
    os = generate_model_predictions(ots_minb1, ots_data);
    ps = generate_model_predictions(pls_minb1, pls_data);
    on = generate_model_predictions(otn_minb1, otn_data);
    pn = generate_model_predictions(pln_minb1, pln_data);
    beta_runs{j} = getPredictionGLMVar(os, ps, on, pn, cohortID);
end

% Retain the original script's repetition-specific beta variables, For continued use of existing analysis code.
os_isc1 = cellfun(@(b) b.os_isc1, beta_runs, 'UniformOutput', false);
ps_isc1 = cellfun(@(b) b.ps_isc1, beta_runs, 'UniformOutput', false);
on_isc1 = cellfun(@(b) b.on_isc1, beta_runs, 'UniformOutput', false);
pn_isc1 = cellfun(@(b) b.pn_isc1, beta_runs, 'UniformOutput', false);
os_conf1 = cellfun(@(b) b.os_conf1, beta_runs, 'UniformOutput', false);
ps_conf1 = cellfun(@(b) b.ps_conf1, beta_runs, 'UniformOutput', false);
on_conf1 = cellfun(@(b) b.on_conf1, beta_runs, 'UniformOutput', false);
pn_conf1 = cellfun(@(b) b.pn_conf1, beta_runs, 'UniformOutput', false);

%% By participant, remove outliers across repeated simulations for each coefficient, then average
% The three columns in each field are, in order:: intercept, evidence-sum coefficient, evidence-difference coefficient.
pred_beta = struct();
names = fieldnames(beta_runs{1});
for f = 1:numel(names)
    field = names{f};
    values = cellfun(@(b) b.(field), beta_runs, 'UniformOutput', false);
    values = cat(3, values{:});       % participant × coefficient × repetition
    averaged = nan(size(values, 1), 3);
    for i = 1:size(values, 1)
        for k = 1:3
            samples = reshape(values(i, k, :), 1, []);
            averaged(i, k) = mean(rmoutliers(samples));
        end
    end
    pred_beta.(field) = averaged;
end
save(fullfile(output_dir, 'all_GLM_beta_pred_2k.mat'), 'pred_beta', 'cohortID');

%% Generate conflict-trial predictions for each participant from fitted parameters
function predictions = generate_model_predictions(results, exp_data)
fixed = struct('b0', 300, 'delta', 0.01, 'plapse', 0.05 / 8);
predictions = struct([]);
for i = 1:numel(results)
    % Original model parameter order: v0, bs_ratio, bp_ratio, ks, kp, sigma, lambda, tau, t0.
    params = results(i).optParams([1:6, 8, 9, 7]);
    predictions = [predictions; predictSubject(exp_data(i), params, fixed)];
end
end

function prediction = predictSubject(data, params, fixed)
% Extract only conflict-trial variables used for prediction.
conflict = data.conflict == 1;
selfConf = data.conf1(conflict);
initialChoice = data.answer1(conflict);
actualChoice = data.answer2(conflict);
actualConf = data.conf2(conflict);
actualRT = data.rtime2(conflict) * 1000;
if isfield(data, 'inforate')
    otherConf = data.inforate(conflict);
else
    otherConf = data.o1conf1(conflict);
end

names = {'SelfConfidence', 'OtherConfidence', 'InitialChoice', ...
    'PredictedChoice', 'PredictedConfidence', 'PredictedRT', ...
    'ActualFinalChoice', 'ActualFinalConfidence', 'ActualRT', 'ConflictFlag'};
rows = nan(numel(selfConf), 10);
for t = 1:numel(selfConf)
    selfStim = selfConf(t) / 16 - 1/32;
    if otherConf(t) >= 1
        otherStim = otherConf(t) / 16 - 1/32;
    else
        otherStim = otherConf(t) - 0.5;
    end

    [choice, confidence, rt] = runPrediction(params, fixed, selfStim, otherStim);
    if rand() < fixed.plapse * 8
        choice = randperm(2, 1);
        confidence = (rand() + 1) / 2;
    end
    confidence = discretizeConfidence(confidence) / 8 + 0.5 - 1/16;

    % accumulator 1 means retaining the initial choice, accumulator 2 means switching choice.
    finalChoice = initialChoice(t);
    if choice == 2
        finalChoice = 3 - initialChoice(t);
    end
    rows(t, :) = [selfConf(t), otherConf(t), initialChoice(t), ...
        finalChoice, confidence, rt, actualChoice(t), actualConf(t), actualRT(t), 1];
end

prediction = struct();
for f = 1:numel(names)
    prediction.(names{f}) = rows(:, f);
end
end

function [choice, confidence, rt] = runPrediction(params, fixed, selfStim, otherStim)
v0 = params(1);
bs = fixed.b0 * max(0, params(2));
bp = fixed.b0 * max(0, params(3));
ks = params(4);
kp = params(5);
sigma = params(6);
tau = params(7);
t0 = params(8);
lambda = params(9);

vs = v0 + ks * selfStim;
vp = v0 + kp * otherStim;
xs = bs;
xp = bp;
rt = 1000;
choice = randi(2);

% Preserve the original simulation rules: at most 1000 steps, Check boundary crossing before updating evidence, RT=step count.
for step = 1:1000
    boundary = fixed.b0;
    if step >= t0
        boundary = fixed.b0 * (1 - 0.9 * (step - t0) / (step - t0 + tau));
    end
    if xs >= boundary
        choice = 1;
        rt = step;
        break;
    elseif xp >= boundary
        choice = 2;
        rt = step;
        break;
    end

    noiseSelf = sigma * randn;
    noiseOther = sigma * randn;
    dxs = -lambda * xp + (vs - fixed.delta * xs) + noiseSelf;
    dxp = -lambda * xs + (vp - fixed.delta * xp) + noiseOther;
    xs = max(0, xs + dxs);
    xp = max(0, xp + dxp);
end

% Compute confidence from the difference between accumulators at termination.
es = min(0.5 + xs / (2 * fixed.b0), 1 - 8 * fixed.plapse);
ep = min(0.5 + xp / (2 * fixed.b0), 1 - 8 * fixed.plapse);
es = max(min(es, 1 - eps), eps);
ep = max(min(ep, 1 - eps), eps);
difference = (log(es) - log(1 - es)) - (log(ep) - log(1 - ep));
if choice == 2
    difference = -difference;
end
confidence = 1 / (1 + exp(-difference));
end

function level = discretizeConfidence(confidence)
confidence = max(min(confidence, 1), 0);
if confidence < 0.625
    level = 1;
elseif confidence < 0.75
    level = 2;
elseif confidence < 0.875
    level = 3;
else
    level = 4;
end
end

%% Eight mixed-effects models: Four experimental conditions × choice / confidence
function beta = getPredictionGLMVar(os, ps, on, pn, cohortID)
groups = {os, ps, on, pn};
prefixes = {'os', 'ps', 'on', 'pn'};
n = numel(os);
cohortID = cohortID(:);
assert(isnumeric(cohortID) && numel(cohortID) == n && all(isfinite(cohortID)), ...
    'cohortID must match the number and order of predicted subjects.');
beta = struct();
for g = 1:4
    assert(numel(groups{g}) == n, 'Condition subject counts differ.');
    data = addPredictionVariables(groups{g}, g <= 2);
    beta.([prefixes{g}, '_isc1']) = fitPredictionGLM( ...
        data, 'isstay', {'lo_plus', 'lo_minus'}, true, cohortID);
    beta.([prefixes{g}, '_conf1']) = fitPredictionGLM( ...
        data, 'lo_ps2', {'lo_cplus', 'lo_cminus'}, false, cohortID);
end
end

function data = addPredictionVariables(data, isSocial)
logOdds = @(p) log(p ./ (1 - p));
for i = 1:numel(data)
    d = data(i);
    changed = d.InitialChoice ~= d.PredictedChoice;
    selfOdds = logOdds(d.SelfConfidence / 16 + 0.5 - 1/32);
    if isSocial
        otherOdds = logOdds(d.OtherConfidence / 16 + 0.5 - 1/32);
    else
        otherOdds = logOdds(d.OtherConfidence);
    end
    chosen = selfOdds;
    unchosen = otherOdds;
    chosen(changed) = otherOdds(changed);
    unchosen(changed) = selfOdds(changed);

    data(i).ischange = changed;
    data(i).isstay = double(~changed);
    data(i).lo_ps2 = logOdds(d.PredictedConfidence);
    data(i).lo_plus = selfOdds + otherOdds;
    data(i).lo_minus = selfOdds - otherOdds;
    data(i).lo_cplus = chosen + unchosen;
    data(i).lo_cminus = chosen - unchosen;
end
end

function beta = fitPredictionGLM(data, response, predictors, isChoice, cohortID)
% Process predictors within each participant, then combine into the fitting table.
blocks = cell(numel(data), 1);
for i = 1:numel(data)
    y = data(i).(response)(:, 1);
    x = [data(i).(predictors{1})(:, 1), data(i).(predictors{2})(:, 1)];
    if isChoice
        valid = ~isnan(y) & ~any(isnan(x), 2);
        y = y(valid);
        x = zscore(x(valid, :), 0, 1);
    else
        % As in the original linearmixedmodel_demean consistent: First exclude only y 's NaN.
        valid = ~isnan(y);
        y = y(valid);
        x = x(valid, :);
        x = x - mean(x, 1, 'omitnan');
    end
    blocks{i} = table(y, repmat(i, numel(y), 1), x(:, 1), x(:, 2), ...
        repmat(cohortID(i), numel(y), 1), ...
        'VariableNames', {response, 'id', predictors{1}, predictors{2}, 'exp_id'});
end
tbl = vertcat(blocks{:});

% Multiple cohort: numeric cohort fixed effect + cohort random intercept.
% Single cohort: Omit the constant cohort term; Always retain participant-specific random intercepts and slopes.
terms = strjoin(predictors, ' + ');
if numel(unique(cohortID)) > 1
    formula = [response, ' ~ ', terms, ' + exp_id + ( ', terms, ' | id ) + (1 | exp_id)'];
else
    formula = [response, ' ~ ', terms, ' + ( ', terms, ' | id )'];
end
if isChoice
    model = fitglme(tbl, formula, 'Distribution', 'Binomial', 'Link', 'logit');
else
    model = fitlme(tbl, formula);
end

% Participant by participant, Match coefficients by name, Extract only id random effects, Do not include cohort coefficient.
[re, ~, stats] = randomEffects(model);
fixed = fixedEffects(model);
beta = nan(numel(data), 3);
coefficientNames = {'(Intercept)', predictors{1}, predictors{2}};
subjectID = str2double(string(stats.Level));
for i = 1:numel(data)
    for j = 1:3
        ri = strcmp(stats.Group, 'id') & subjectID == i & ...
            strcmp(stats.Name, coefficientNames{j});
        fi = strcmp(model.CoefficientNames, coefficientNames{j});
        assert(sum(ri) == 1 && sum(fi) == 1, 'Cannot map subject coefficient.');
        beta(i, j) = fixed(fi) + re(ri);
    end
end
end
