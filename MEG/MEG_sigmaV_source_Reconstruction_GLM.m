%% trajectory_sourceEFR_GLM_beta_ROIgamma30_100hz_reorganized
% For each participant, at each time point perform ROI source signals GLM, Save only the final required variables.
%
% Output variables: 
%   sb_chosen0, cb_chosen0       stay/change, Including condition predictors
%   sb_chosen, cb_chosen          stay/change, self/other predictors
%   otsbs_chosen, otsbc_chosen    condition 1 's stay/change GLM
%   mb_chosen                     stay + change Balanced mixed sampling GLM mean
%   otsb_chosen1, plsb_chosen1,
%   otnb_chosen1, plnb_chosen1   Four condition 's chosen/unchosen
%                                  sum/difference GLM
%   vchan_time                    time axis
%
% Each GLM Coefficient order is: 
%   [intercept, predictor_1, predictor_2, ...]

clearvars -except cfg_user
clc;

%% configuration
cfg = struct();

% Gamma-band ROI Source-result directory
cfg.data_path = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
    'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\SigmaV_timeBetaFilter_vectorMagnitude_ROI'];

cfg.output_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT_all\MEG\figure_data\trajectory_sourceEFR_GLM_beta_sumq.mat';

cfg.roi_idx = 2;
cfg.analysis_window = [-1 0.5];

% Ordinary GLM minimum trial count.Set to 3 Ensures sufficient degrees of freedom with two predictors.
cfg.min_glm_trials = 3;

% Consistent with the original second code section: sb_chosen/cb_chosen Retain only trials with at least this within the analysis window 20 items trial.
cfg.min_stay_change_trials = 20;

% sb_chosen, cb_chosen, mb_chosen Whether to use only conflict==1 's trial.
% Consistent with the original second section's definition.
cfg.require_conflict = true;

% mb_chosen 's stay/change Balanced-sampling parameters
cfg.n_subsample = 5000;
cfg.mixed_fraction = 0.90;
cfg.min_mixed_trials = 20;
cfg.base_seed = 100000;
cfg.use_parallel = true;

if exist('cfg_user','var') && isstruct(cfg_user)
    cfg = merge_cfg_local(cfg,cfg_user);
end

assert(isfolder(cfg.data_path), 'data_path 不存在：%s', cfg.data_path);

if ~isfolder(fileparts(cfg.output_file))
    mkdir(fileparts(cfg.output_file));
end

files = dir(fullfile(cfg.data_path,'*.mat'));
assert(~isempty(files), 'data_path 中没有 MAT 文件：%s', cfg.data_path);
[~,order] = sort(lower(string({files.name})));
files = files(order);
n_sub = numel(files);

fprintf('输入被试数：%d\n', n_sub);
fprintf('ROI：%d；分析时间窗：[%.3f %.3f] s\n', ...
    cfg.roi_idx,cfg.analysis_window(1),cfg.analysis_window(2));

%% Read and validate time axes and data dimensions
first = load(fullfile(files(1).folder,files(1).name),'roi_amp');
assert(isfield(first,'roi_amp'),'文件缺少 roi_amp：%s',files(1).name);
check_roi_amp_local(first.roi_amp,files(1).name,cfg.roi_idx);

vchan_time = double(first.roi_amp.time(:)');
n_time = numel(vchan_time);
time_keep = vchan_time >= cfg.analysis_window(1) & ...
    vchan_time <= cfg.analysis_window(2);
sample_idx = find(time_keep);
assert(numel(sample_idx)>=3,'analysis_window 内时间点少于 3 个。');

% Preallocate arrays required for saving
sb_chosen0 = nan(n_sub,n_time,4);  % intercept + s1 + o1 + condition
cb_chosen0 = nan(n_sub,n_time,4);
sb_chosen  = nan(n_sub,n_time,3);  % intercept + s1 + o1
cb_chosen  = nan(n_sub,n_time,3);

otsbs_chosen = nan(n_sub,n_time,3);
otsbc_chosen = nan(n_sub,n_time,3);

otsb_chosen1 = nan(n_sub,n_time,3);
plsb_chosen1 = nan(n_sub,n_time,3);
otnb_chosen1 = nan(n_sub,n_time,3);
plnb_chosen1 = nan(n_sub,n_time,3);

mb_chosen = nan(n_sub,n_time,3);

%% Compute participant by participant
for i_sub = 1:n_sub
    file_now = fullfile(files(i_sub).folder,files(i_sub).name);
    S = load(file_now,'roi_amp');
    assert(isfield(S,'roi_amp'),'文件缺少 roi_amp：%s',files(i_sub).name);
    roi_amp = S.roi_amp;
    check_roi_amp_local(roi_amp,files(i_sub).name,cfg.roi_idx);

    this_time = double(roi_amp.time(:)');
    assert(numel(this_time)==n_time && ...
        max(abs(this_time-vchan_time),[],'omitnan')<1e-8, ...
        '被试 %d 的时间轴与第一个被试不一致。',i_sub);

    bhv = roi_amp.behavior;
    A = double(roi_amp.unitnoisegain_no_baseline);
    Y = reshape(A(:,cfg.roi_idx,:),size(A,1),n_time);

    % Convert behavioral variables to column vectors
    lo_chosen   = double(bhv.lo_chosen(:));
    lo_unchosen = double(bhv.lo_unchosen(:));
    lo_s1       = double(bhv.lo_s1(:));
    lo_o1       = double(bhv.lo_o1(:));
    condition   = double(bhv.condition(:));
    ischange    = double(bhv.ischange(:));

    base_valid = isfinite(lo_chosen) & isfinite(lo_unchosen) & ...
        isfinite(lo_s1) & isfinite(lo_o1) & ...
        isfinite(condition) & isfinite(ischange);

    % Original first-half analysis: No additional restriction on conflict
    % stay/change + condition：
    % [z(lo_s1), z(lo_o1), condition]
    idx_stay = base_valid & ischange==0;
    idx_change = base_valid & ischange==1;

    X_stay0 = [lo_s1(idx_stay),lo_o1(idx_stay),condition(idx_stay)];
    X_change0 = [lo_s1(idx_change),lo_o1(idx_change),condition(idx_change)];

    sb_chosen0(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_stay,:),X_stay0, ...
        [true true false],cfg.min_glm_trials,n_time),n_time);
    cb_chosen0(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_change,:),X_change0, ...
        [true true false],cfg.min_glm_trials,n_time),n_time);

    % condition 1 's stay/change chosen GLM
    idx_otsbs = base_valid & condition==1 & ischange==0;
    idx_otsbc = base_valid & condition==1 & ischange==1;
    otsbs_chosen(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_otsbs,:), ...
        [lo_s1(idx_otsbs),lo_o1(idx_otsbs)],true(1,2), ...
        cfg.min_glm_trials,n_time),n_time);
    otsbc_chosen(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_otsbc,:), ...
        [lo_s1(idx_otsbc),lo_o1(idx_otsbc)],true(1,2), ...
        cfg.min_glm_trials,n_time),n_time);

    % Four condition 's chosen/unchosen sum/difference GLM
    cond_id = [1 2 3 4];
    cond_output = {'otsb_chosen1','plsb_chosen1','otnb_chosen1','plnb_chosen1'};
    for k = 1:numel(cond_id)
        idx_cond = base_valid & condition==cond_id(k);
        x_sum = lo_chosen(idx_cond) + lo_unchosen(idx_cond);
        x_diff = lo_chosen(idx_cond) - lo_unchosen(idx_cond);
        beta_cond = fit_glm_timewise_local(Y(idx_cond,:), ...
            [x_sum,x_diff],true(1,2),cfg.min_glm_trials,n_time);
        beta_cond = place_time_window_local(beta_cond,n_time);
        switch cond_output{k}
            case 'otsb_chosen1', otsb_chosen1(i_sub,:,:) = beta_cond;
            case 'plsb_chosen1', plsb_chosen1(i_sub,:,:) = beta_cond;
            case 'otnb_chosen1', otnb_chosen1(i_sub,:,:) = beta_cond;
            case 'plnb_chosen1', plnb_chosen1(i_sub,:,:) = beta_cond;
        end
    end

    % Second analysis section: stay/change Use only conflict==1 's trial
    % The original second section requires a trial throughout the entire analysis_window to contain finite signals.
    signal_valid_window = all(isfinite(Y(:,sample_idx)),2);
    analysis_valid = base_valid & signal_valid_window;
    if cfg.require_conflict
        assert(ismember('conflict',bhv.Properties.VariableNames), ...
            '文件 %s 缺少 conflict 字段。',files(i_sub).name);
        conflict = double(bhv.conflict(:));
        analysis_valid = analysis_valid & isfinite(conflict) & conflict==1;
    end

    idx_stay_a = analysis_valid & ischange==0;
    idx_change_a = analysis_valid & ischange==1;

    sb_chosen(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_stay_a,sample_idx), ...
        [lo_s1(idx_stay_a),lo_o1(idx_stay_a)],true(1,2), ...
        cfg.min_stay_change_trials,numel(sample_idx)), ...
        sample_idx,n_time);
    cb_chosen(i_sub,:,:) = place_time_window_local( ...
        fit_glm_timewise_local(Y(idx_change_a,sample_idx), ...
        [lo_s1(idx_change_a),lo_o1(idx_change_a)],true(1,2), ...
        cfg.min_stay_change_trials,numel(sample_idx)), ...
        sample_idx,n_time);

    % mb_chosen: in conflict After filtering stay/change trial sample equal numbers, 
    % Fit a common model for each draw self/other GLM, Finally average over 5000 draws.
    n_stay = sum(idx_stay_a);
    n_change = sum(idx_change_a);
    n_each = floor(cfg.mixed_fraction*min(n_stay,n_change));
    if n_each < cfg.min_mixed_trials
        warning('被试 %d 的混合抽样 trial 不足，mb_chosen 保持 NaN。',i_sub);
        continue;
    end

    stay_rows = find(idx_stay_a);
    change_rows = find(idx_change_a);
    seed = cfg.base_seed + i_sub*10000;
    stream = RandStream('Threefry','Seed',seed);
    stay_draw = zeros(cfg.n_subsample,n_each,'uint32');
    change_draw = zeros(cfg.n_subsample,n_each,'uint32');
    for d = 1:cfg.n_subsample
        stay_draw(d,:) = uint32(randperm(stream,n_stay,n_each));
        change_draw(d,:) = uint32(randperm(stream,n_change,n_each));
    end

    Y_mix = Y(:,sample_idx);
    beta_draw = nan(cfg.n_subsample,numel(sample_idx),3,'single');
    min_for_mix = cfg.min_mixed_trials;

    if cfg.use_parallel
        parfor d = 1:cfg.n_subsample
            s = stay_rows(double(stay_draw(d,:)));
            c = change_rows(double(change_draw(d,:)));
            rows = [s;c];
            X_mix = [lo_s1(rows),lo_o1(rows)];
            beta_d = fit_glm_timewise_local(Y_mix(rows,:),X_mix, ...
                true(1,2),min_for_mix,numel(sample_idx));
            beta_draw(d,:,:) = single(beta_d);
        end
    else
        for d = 1:cfg.n_subsample
            s = stay_rows(double(stay_draw(d,:)));
            c = change_rows(double(change_draw(d,:)));
            rows = [s;c];
            X_mix = [lo_s1(rows),lo_o1(rows)];
            beta_d = fit_glm_timewise_local(Y_mix(rows,:),X_mix, ...
                true(1,2),min_for_mix,numel(sample_idx));
            beta_draw(d,:,:) = single(beta_d);
        end
    end

    beta_mean = squeeze(mean(beta_draw,1,'omitnan'));
    if isvector(beta_mean)
        beta_mean = reshape(beta_mean,numel(sample_idx),3);
    end
    mb_full = nan(n_time,3);
    mb_full(sample_idx,:) = double(beta_mean);
    mb_chosen(i_sub,:,:) = reshape(mb_full,1,n_time,3);

    fprintf('被试 %d/%d 完成：stay=%d，change=%d，mixed=%d，draw=%d。\n', ...
        i_sub,n_sub,n_stay,n_change,n_each,cfg.n_subsample);
end

%% Final checks and save specified variables only
assert(numel(vchan_time)==n_time,'vchan_time 长度异常。');
assert(size(sb_chosen0,1)==n_sub && size(sb_chosen0,2)==n_time, ...
    '输出变量尺寸异常。');

save(cfg.output_file, ...
    'sb_chosen0','sb_chosen','cb_chosen0','cb_chosen', ...
    'otsbs_chosen','otsbc_chosen','mb_chosen', ...
    'otsb_chosen1','plsb_chosen1','otnb_chosen1','plnb_chosen1', ...
    'vchan_time','-v7.3');

fprintf('\n已保存：%s\n',cfg.output_file);
fprintf('保存字段：sb_chosen0, sb_chosen, cb_chosen0, cb_chosen,\n');
fprintf('          otsbs_chosen, otsbc_chosen, mb_chosen,\n');
fprintf('          otsb_chosen1, plsb_chosen1, otnb_chosen1, plnb_chosen1, vchan_time\n');

%% Local functions
function out = merge_cfg_local(defaults,user)
out = defaults;
names = fieldnames(user);
for i = 1:numel(names)
    out.(names{i}) = user.(names{i});
end
end

function check_roi_amp_local(roi_amp,file_name,roi_idx)
required = {'time','behavior','unitnoisegain_no_baseline'};
assert(all(isfield(roi_amp,required)), ...
    '文件 %s 缺少 roi_amp 必需字段。',file_name);
assert(istable(roi_amp.behavior), ...
    '文件 %s 中 roi_amp.behavior 必须是 table。',file_name);
assert(ndims(roi_amp.unitnoisegain_no_baseline)==3, ...
    '文件 %s 的 unitnoisegain_no_baseline 必须是 trial×ROI×time。',file_name);
assert(roi_idx<=size(roi_amp.unitnoisegain_no_baseline,2), ...
    '文件 %s 不包含 ROI %d。',file_name,roi_idx);
assert(size(roi_amp.unitnoisegain_no_baseline,1)==height(roi_amp.behavior), ...
    '文件 %s 的 trial 数与 behavior 行数不一致。',file_name);
required_bhv = {'lo_chosen','lo_unchosen','lo_s1','lo_o1', ...
    'condition','ischange'};
assert(all(ismember(required_bhv,roi_amp.behavior.Properties.VariableNames)), ...
    '文件 %s 缺少一个或多个行为字段。',file_name);
end

function beta = fit_glm_timewise_local(Y,X,standardize_mask,min_trials,n_time)
% Fit separately at each time point GLM.
% Y：trial×time；X：trial×predictor。
% standardize_mask Specify which predictors in currently valid trial should undergo z-score.
Y = double(Y);
X = double(X);
if nargin<5 || isempty(n_time)
    n_time = size(Y,2);
end

n_predictor = size(X,2);
beta = nan(n_time,n_predictor+1);
assert(size(Y,1)==size(X,1),'Y 与 X 的 trial 数不一致。');
assert(numel(standardize_mask)==n_predictor, ...
    'standardize_mask 长度与预测变量数不一致。');

for t = 1:min(n_time,size(Y,2))
    valid = isfinite(Y(:,t)) & all(isfinite(X),2);
    if sum(valid)<min_trials
        continue;
    end

    Xt = X(valid,:);
    yt = Y(valid,t);

    % Preserve the original definition: Y At each time point z-score; 
    % Predictors in specified columns within currently valid trial observations z-score.
    for j = 1:n_predictor
        if standardize_mask(j)
            mu = mean(Xt(:,j),'omitnan');
            sd = std(Xt(:,j),0,'omitnan');
            if ~isfinite(sd) || sd<=eps
                Xt = [];
                break;
            end
            Xt(:,j) = (Xt(:,j)-mu)./sd;
        end
    end
    if isempty(Xt)
        continue;
    end

    y_mu = mean(yt,'omitnan');
    y_sd = std(yt,0,'omitnan');
    if ~isfinite(y_sd) || y_sd<=eps
        continue;
    end
    yt = (yt-y_mu)./y_sd;

    D = [ones(size(Xt,1),1),Xt];
    if rank(D)<size(D,2)
        continue;
    end
    beta(t,:) = (D\yt)';
end
end

function beta_full = place_time_window_local(beta,n_time_or_idx,varargin)
% Place full-time-axis or analysis-window results back into n_time time points.
if isempty(varargin)
    n_time = n_time_or_idx;
    assert(size(beta,1)==n_time,'GLM 输出时间维度不匹配。');
    beta_full = beta;
    return;
end

sample_idx = n_time_or_idx;
n_time = varargin{1};
assert(size(beta,1)==numel(sample_idx), ...
    '分析窗 GLM 输出长度与 sample_idx 不一致。');
beta_full = nan(n_time,size(beta,2));
beta_full(sample_idx,:) = beta;
end
