%% 
% cross test part1: conjunction map(gamma-band source ∩ oxytocin receptor expression map) → 
% coefficient(ΣV) source filter → single reconstraction → trajectory GLM
%% ============================== Configuration =========================
cfg = struct();
cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(cfg.ft_path);
ft_defaults;

% Significant source cluster produced by the SigmaV-beta source analysis.
cfg.stat_file = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\' ...
    'Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)\stat_cluster.mat'];
cfg.cluster_id = [];                 % [] = significant positive cluster with minimum p
cfg.cluster_alpha = 0.01;

% Which ROI time courses to save: cluster | conjunction | both
cfg.roi_output_mode = 'both';
cfg.target_roi_file = ['E:\xianliang\matlab_m\social_decision_m\data\' ...
    'OT-N2-MEG\MEG_data\MEG_mask\article\oxtr_mask\oxtr_50.nii'];
cfg.target_roi_parameter = 'anatomy';
cfg.target_roi_threshold = 0.5;
cfg.aal_file = fullfile(cfg.ft_path,'template','atlas','aal', ...
    'ROI_MNI_V4.nii');

% Raw MEG/behavior/leadfield inputs. The script rebuilds the vector filter
% from the same target/ITI time-resolved beta_SigmaV used by source analysis;
% it does not reuse a scalar, raw-covariance or ROI-SVD filter.
cfg.num_subjects = 19;              % [] = all subjects in behavior file
cfg.base_data_path = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MEG-Preproocess2\pressPutton'];
cfg.behavior_file = ['E:\xianliang\matlab_m\social_decision_m\data\' ...
    'OT-N2-MEG\behavior\analysis_data\behaviour_meg.mat'];
cfg.base_leadfield_path = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MRI_Preprocess\lead_field_pressPutton_3run(250hz)'];
cfg.head_movement_file = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MRI_Preprocess\out_run_4condition\checkHead_idx_2day.mat'];
cfg.template_grid_file = fullfile(cfg.ft_path,'template', ...
    'sourcemodel','standard_sourcemodel3d5mm.mat');
cfg.analysis_window = [-7.2 1];
cfg.trials_per_run = 80;
cfg.behavior_fields = {'conflict','ischange','lo_s1','lo_o1', ...
    'condition','run_id','run_uid','oiti','osti2','rtime2', ...
    'original_trial_id','lo_chosen','lo_unchosen','osti'};
cfg.condition_dirs = {
    fullfile(cfg.base_data_path,'social_OT_306(250hz)')
    fullfile(cfg.base_data_path,'social_PL_306(250hz)')
    fullfile(cfg.base_data_path,'nonsocial_OT_306(250hz)')
    fullfile(cfg.base_data_path,'nonsocial_PL_306(250hz)')};
cfg.leadfield_dirs = {
    fullfile(cfg.base_leadfield_path,'social-OT(250hz)')
    fullfile(cfg.base_leadfield_path,'social-PL(250hz)')
    fullfile(cfg.base_leadfield_path,'nonsocial-OT(250hz)')
    fullfile(cfg.base_leadfield_path,'nonsocial-PL(250hz)')};

% Filter training: must match the SigmaV source-localization configuration.
cfg.filter_target_window = [-0.5 -0.3];
cfg.filter_iti_delay = 0.0;
cfg.filter_iti_duration = 0.2;
cfg.filter_min_iti_duration = 0.19;

% GLM and beta covariance settings. Y is intentionally not z-scored.
% Predictor scaling is shared across OT/PL, so the two daily beta maps have
% the same physical sensor units per common predictor SD.
cfg.beta_filter_mode = 'sigmaV_target_iti_common';
cfg.beta_covariance_mode = 'secondmoment'; % secondmoment | covariance
cfg.beta_smooth_samples = 1;         % smooth beta(t) before covariance; 1=no
cfg.add_change_nuisance = false;
cfg.add_condition_nuisance = false;
cfg.add_run_nuisance = false;
cfg.regularization_mode = 'per_day'; % per_day recommended | joint_global sensitivity
cfg.lambda = 0.10;
cfg.min_rcond = 1e-10;
cfg.min_filter_trial_per_day = 20;
cfg.apply_day_gain_calibration = true;

% Projection/reconstruction: can differ from filter windows.
cfg.reconstruction_window = [-1.5 1.0];

% Trial-matched ITI amplitude baseline. This does not affect the no-baseline output.
cfg.baseline_iti_delay = cfg.filter_iti_delay;
cfg.baseline_iti_duration = cfg.filter_iti_duration;
cfg.baseline_min_iti_duration = cfg.filter_min_iti_duration;

% Match sensor sqrt(g1^2+g2^2). 'power' is only a sensitivity analysis.
cfg.amplitude_measure = 'magnitude'; % magnitude | power
cfg.log_transform = false;
cfg.log_offset_ratio = 1e-6;

% Downsample only after source projection and magnitude calculation.
cfg.fs_new = 100;
cfg.filter_storage_class = 'single';
cfg.overwrite = false;
cfg.output_dir = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\source_valid\cross_test\' ...
    'SigmaV_timeBetaFilter_vectorMagnitude_ROI_rank3_UG0_UNG1'];

cfg.beta_filter_mode = char(lower(string(cfg.beta_filter_mode)));
assert(strcmp(cfg.beta_filter_mode,'sigmav_target_iti_common'), ...
    'This reconstruction requires cfg.beta_filter_mode=sigmaV_target_iti_common.');
cfg.beta_covariance_mode = char(lower(string(cfg.beta_covariance_mode)));
assert(ismember(cfg.beta_covariance_mode,{'secondmoment','covariance'}));
assert(cfg.beta_smooth_samples>=1 && ...
    mod(cfg.beta_smooth_samples,1)==0);
assert(ismember(cfg.regularization_mode,{'joint_global','per_day'}));
assert(cfg.lambda>0 && cfg.min_rcond>0);
assert(cfg.filter_iti_duration>=cfg.filter_min_iti_duration);
assert(ismember(cfg.filter_storage_class,{'single','double'}));
assert(isempty(regexp(cfg.output_dir,'[<>|?*"]','once')), ...
    'cfg.output_dir contains a Windows-invalid filename character.');

%% ======================= Load cluster and define ROI(s) ================
assert(isfile(cfg.stat_file),'Cannot find stat file: %s',cfg.stat_file);
Sstat = load(cfg.stat_file,'stat_res');
assert(isfield(Sstat,'stat_res'),'%s lacks stat_res.',cfg.stat_file);
stat_res = Sstat.stat_res;
[cluster_idx,cluster_info] = select_significant_cluster_local( ...
    stat_res,cfg.cluster_id,cfg.cluster_alpha);

roi_defs = struct('label',{},'grid_idx',{});
switch lower(cfg.roi_output_mode)
    case 'cluster'
        conjunction_idx = [];
        roi_defs(1).label = 'SigmaV_cluster_vectorMagnitude';
        roi_defs(1).grid_idx = cluster_idx(:);
    case 'conjunction'
        [conjunction_idx,target_roi_info] = ...
            make_target_roi_conjunction_local(stat_res,cluster_idx,cfg);
        cluster_info.target_roi = target_roi_info;
        roi_defs(1).label = 'SigmaV_cluster_x_OXTR_vectorMagnitude';
        roi_defs(1).grid_idx = conjunction_idx(:);
    case 'both'
        [conjunction_idx,target_roi_info] = ...
            make_target_roi_conjunction_local(stat_res,cluster_idx,cfg);
        cluster_info.target_roi = target_roi_info;
        roi_defs(1).label = 'SigmaV_cluster_vectorMagnitude';
        roi_defs(1).grid_idx = cluster_idx(:);
        roi_defs(2).label = 'SigmaV_cluster_x_OXTR_vectorMagnitude';
        roi_defs(2).grid_idx = conjunction_idx(:);
    otherwise
        error('Unknown cfg.roi_output_mode: %s',cfg.roi_output_mode);
end
requested_grid_idx = unique(vertcat(roi_defs.grid_idx),'stable');
cluster_info.conjunction_idx = conjunction_idx;
cluster_info.roi_output_mode = lower(cfg.roi_output_mode);
cluster_info.roi_output_labels = {roi_defs.label};

fprintf('\nSelected cluster %d: p=%.6f, %d points.\n', ...
    cluster_info.cluster_id,cluster_info.cluster_p,numel(cluster_idx));
fprintf('Requested vector filters: %d grid points; ROI outputs: %d.\n', ...
    numel(requested_grid_idx),numel(roi_defs));

%% ======================== Load subjects/run metadata ===================
cfg_main = cfg;
B = load(cfg_main.behavior_file);
required_behavior = {'ots_meg','pls_meg','otn_meg','pln_meg'};
for i = 1:numel(required_behavior)
    assert(isfield(B,required_behavior{i}), ...
        'Behavior file lacks %s.',required_behavior{i});
end
data_bhv = {B.ots_meg;B.pls_meg;B.otn_meg;B.pln_meg};
subject_ids = get_and_verify_subject_ids_local(data_bhv);
if ~isempty(cfg.num_subjects)
    assert(cfg.num_subjects<=numel(subject_ids), ...
        'cfg.num_subjects exceeds the number in the behavior file.');
    subject_ids = subject_ids(1:cfg.num_subjects);
end

head_mov = load(cfg_main.head_movement_file, ...
    'optimal_idx','isolated_bad_run');
% Project-specific audited corrections retained from the current script.
% Project-specific audited correction retained from the latest script.
head_mov.isolated_bad_run.ot{18} = [3 6];
assert(isfield(head_mov,'optimal_idx') && ...
    isfield(head_mov,'isolated_bad_run'));

source_grid = load_source_grid_local(cfg_main.template_grid_file);
source_grid = ft_convert_units(source_grid,'mm');
verify_stat_source_geometry_local(stat_res,source_grid);
if ~exist(cfg.output_dir,'dir'),mkdir(cfg.output_dir);end

all_roi_file = cell(numel(subject_ids),1);
all_quality = cell(numel(subject_ids),1);

%% ================================ Subjects =============================
for i_sub = 1:numel(subject_ids)
    subject_id = subject_ids{i_sub};
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] Time-beta-filtered vector-amplitude ROI: %s\n', ...
        i_sub,numel(subject_ids),subject_id);

    [ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf,load_quality] = ...
        prepare_subject_local(i_sub,subject_id,data_bhv,head_mov, ...
        cfg_main,cfg_main.analysis_window);
    [ot_data,label_ot,lf_idx_ot] = align_data_lf_local( ...
        ot_data,ot_lf,'MEGGRAD','OT');
    [pl_data,label_pl,lf_idx_pl] = align_data_lf_local( ...
        pl_data,pl_lf,'MEGGRAD','PL');
    assert(numel(label_ot)==204 && numel(label_pl)==204, ...
        'Expected 204 gradiometers/day; found OT=%d, PL=%d.', ...
        numel(label_ot),numel(label_pl));

    [filter_bank,filter_quality] = build_vector_filter_bank_local( ...
        ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf, ...
        lf_idx_ot,lf_idx_pl,requested_grid_idx,cfg);

    [roi_amp,amp_quality] = reconstruct_vector_roi_local( ...
        ot_data,pl_data,bhv_ot,bhv_pl,filter_bank,roi_defs,cfg);

    % Convenience FieldTrip structures. vchan is kept as a backward-
    % compatible alias of the recommended no-baseline amplitude output.
    vchan_no_baseline = array_to_ft_raw_local( ...
        roi_amp.unitnoisegain_no_baseline, ...
        roi_amp.time,roi_amp.label,roi_amp.behavior,roi_amp.fsample);
    vchan_iti_corrected = array_to_ft_raw_local( ...
        roi_amp.unitnoisegain_iti_corrected, ...
        roi_amp.time,roi_amp.label,roi_amp.behavior,roi_amp.fsample);
    vchan = vchan_no_baseline;

    quality = struct('loading',load_quality,'filter',filter_quality, ...
        'amplitude',amp_quality);
    output_file = fullfile(cfg.output_dir,sprintf( ...
        '%s_SigmaV_TB_R2_%s.mat',subject_id, ...
        lower(cfg.roi_output_mode)));
    assert(numel(output_file)<240, ...
        'Output path is too long for reliable MATLAB v7.3/HDF5 writing.');

    % Recover from an interrupted save that left a zero-byte formal file.
    % A nonzero file remains protected unless cfg.overwrite=true.
    if isfile(output_file)
        old_info = dir(output_file);
        if old_info.bytes==0
            warning('Deleting zero-byte incomplete output: %s',output_file);
            delete(output_file);
        end
    end
    if isfile(output_file) && ~cfg.overwrite
        error('Output exists and cfg.overwrite=false: %s',output_file);
    end

    % Save atomically: complete the HDF5 write under a temporary filename
    % and expose the formal output only after the write succeeds.
    partial_file = fullfile(cfg.output_dir,sprintf( ...
        '.%s_SigmaV_TB_R2_%s.partial.mat', ...
        subject_id,lower(cfg.roi_output_mode)));
    assert(numel(partial_file)<240, ...
        'Temporary output path is too long for reliable v7.3 writing.');
    if isfile(partial_file),delete(partial_file);end
    try
        save(partial_file,'roi_amp','vchan','vchan_no_baseline', ...
            'vchan_iti_corrected','quality','filter_bank', ...
            'cluster_info','cfg','subject_id','-v7.3');
    catch ME
        if isfile(partial_file),delete(partial_file);end
        error('SubjectSave:WriteFailed', ...
            'Failed to write temporary MAT file for %s:\n%s', ...
            subject_id,ME.message);
    end
    [move_ok,move_msg] = movefile(partial_file,output_file,'f');
    assert(move_ok,'Could not finalize %s: %s',output_file,move_msg);

    all_roi_file{i_sub} = output_file;
    all_quality{i_sub} = quality;
    fprintf('Saved: %s\n',output_file);
end

group_file = fullfile(cfg.output_dir,sprintf( ...
    'GROUP_SigmaV_TB_R2_%s.mat', ...
    lower(cfg.roi_output_mode)));
group_partial_file = fullfile(cfg.output_dir,sprintf( ...
    '.GROUP_SigmaV_TB_R2_%s.partial.mat', ...
    lower(cfg.roi_output_mode)));
if isfile(group_partial_file),delete(group_partial_file);end
try
    save(group_partial_file,'all_roi_file','all_quality','subject_ids', ...
        'cluster_info','cfg','-v7.3');
catch ME
    if isfile(group_partial_file),delete(group_partial_file);end
    error('GroupSave:WriteFailed', ...
        'Failed to write temporary group MAT file:\n%s',ME.message);
end
[move_ok,move_msg] = movefile(group_partial_file,group_file,'f');
assert(move_ok,'Could not finalize %s: %s',group_file,move_msg);
fprintf('\nAll subjects completed: %s\n',group_file);

%% ============================= Local functions =========================
function [bank,info] = build_vector_filter_bank_local( ...
    ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf, ...
    lf_idx_ot,lf_idx_pl,grid_idx,cfg)
fs_ot = get_fsample_local(ot_data);
fs_pl = get_fsample_local(pl_data);
assert(abs(fs_ot-fs_pl)<0.5,'OT/PL sampling rates differ.');

% ---------------------------------------------------------------------
% 1. Extract raw target/ITI cubes only for fitting time-resolved GLMs.
%    The raw sensor cubes themselves are NOT used to estimate C.
% ---------------------------------------------------------------------
[Yt_ot,time_t_ot,t_ot] = extract_fixed_cube_local( ...
    ot_data,cfg.filter_target_window,'OT filter target');
[Yt_pl,time_t_pl,t_pl] = extract_fixed_cube_local( ...
    pl_data,cfg.filter_target_window,'PL filter target');
assert_same_time_local(time_t_ot,time_t_pl,'OT/PL filter target');
[Yi_ot,time_i_ot,i_ot] = extract_iti_cube_local( ...
    ot_data,bhv_ot,cfg.filter_iti_duration,cfg.filter_iti_delay, ...
    cfg.filter_min_iti_duration,fs_ot,'OT filter ITI');
[Yi_pl,time_i_pl,i_pl] = extract_iti_cube_local( ...
    pl_data,bhv_pl,cfg.filter_iti_duration,cfg.filter_iti_delay, ...
    cfg.filter_min_iti_duration,fs_pl,'PL filter ITI');
assert_same_time_local(time_i_ot,time_i_pl,'OT/PL filter ITI');

% A fixed trial set and complete ITI windows are required at every beta
% time point. Partial ITIs remain allowed later for amplitude baselining,
% but are not allowed to construct the beta covariance.
finite_predictor_ot = isfinite(bhv_ot.lo_chosen) & ...
    isfinite(bhv_ot.lo_unchosen);
finite_predictor_pl = isfinite(bhv_pl.lo_chosen) & ...
    isfinite(bhv_pl.lo_unchosen);
valid_ot = t_ot.valid & i_ot.valid & i_ot.full_length & ...
    finite_predictor_ot;
valid_pl = t_pl.valid & i_pl.valid & i_pl.full_length & ...
    finite_predictor_pl;
assert(sum(valid_ot)>=cfg.min_filter_trial_per_day && ...
    sum(valid_pl)>=cfg.min_filter_trial_per_day, ...
    ['Too few common target/ITI filter trials: OT=%d, PL=%d; ' ...
     'minimum/day=%d.'],sum(valid_ot),sum(valid_pl), ...
    cfg.min_filter_trial_per_day);

% ---------------------------------------------------------------------
% 2. Build SigmaV/DeltaV and share each predictor's scaling across OT/PL.
%    Y is intentionally not z-scored. Thus beta_SigmaV has physical sensor
%    units per one pooled OT/PL SD of SigmaV for this subject.
% ---------------------------------------------------------------------
sigma_ot = double(bhv_ot.lo_chosen) + double(bhv_ot.lo_unchosen);
delta_ot = double(bhv_ot.lo_chosen) - double(bhv_ot.lo_unchosen);
sigma_pl = double(bhv_pl.lo_chosen) + double(bhv_pl.lo_unchosen);
delta_pl = double(bhv_pl.lo_chosen) - double(bhv_pl.lo_unchosen);

sigma_pool = [sigma_ot(valid_ot);sigma_pl(valid_pl)];
delta_pool = [delta_ot(valid_ot);delta_pl(valid_pl)];
scale_info = struct();
scale_info.sigma_mean = mean(sigma_pool,'omitnan');
scale_info.sigma_sd = std(sigma_pool,0,'omitnan');
scale_info.delta_mean = mean(delta_pool,'omitnan');
scale_info.delta_sd = std(delta_pool,0,'omitnan');
scale_info.definition = [ ...
    'SigmaV=chosen+unchosen; DeltaV=chosen-unchosen; ' ...
    'each predictor scaled with its own pooled OT/PL mean and SD'];
assert(all(isfinite([scale_info.sigma_sd,scale_info.delta_sd])) && ...
    scale_info.sigma_sd>eps && scale_info.delta_sd>eps, ...
    'SigmaV/DeltaV predictor scaling is invalid.');

sigma_z_ot = (sigma_ot-scale_info.sigma_mean)/scale_info.sigma_sd;
delta_z_ot = (delta_ot-scale_info.delta_mean)/scale_info.delta_sd;
sigma_z_pl = (sigma_pl-scale_info.sigma_mean)/scale_info.sigma_sd;
delta_z_pl = (delta_pl-scale_info.delta_mean)/scale_info.delta_sd;

% ---------------------------------------------------------------------
% 3. Fit Y(t) ~ beta_SigmaV(t)*SigmaV + beta_DeltaV(t)*DeltaV.
% ---------------------------------------------------------------------
[bt_sigma_ot,bt_delta_ot,glm_t_ot] = fit_time_glm_local( ...
    Yt_ot,sigma_z_ot,delta_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT filter target');
[bt_sigma_pl,bt_delta_pl,glm_t_pl] = fit_time_glm_local( ...
    Yt_pl,sigma_z_pl,delta_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL filter target');
[bi_sigma_ot,bi_delta_ot,glm_i_ot] = fit_time_glm_local( ...
    Yi_ot,sigma_z_ot,delta_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT filter ITI');
[bi_sigma_pl,bi_delta_pl,glm_i_pl] = fit_time_glm_local( ...
    Yi_pl,sigma_z_pl,delta_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL filter ITI');

if cfg.beta_smooth_samples>1
    beta_cell = {bt_sigma_ot,bt_delta_ot,bt_sigma_pl,bt_delta_pl, ...
        bi_sigma_ot,bi_delta_ot,bi_sigma_pl,bi_delta_pl};
    for i_beta = 1:numel(beta_cell)
        beta_cell{i_beta} = smoothdata(beta_cell{i_beta},2,'movmean', ...
            cfg.beta_smooth_samples,'omitnan');
    end
    [bt_sigma_ot,bt_delta_ot,bt_sigma_pl,bt_delta_pl, ...
        bi_sigma_ot,bi_delta_ot,bi_sigma_pl,bi_delta_pl] = beta_cell{:};
end

% Only beta_SigmaV constructs the spatial filter. beta_DeltaV is retained
% as a diagnostic/nuisance-model coefficient and is not mixed into C_beta.
Bt_ot = bt_sigma_ot;
Bt_pl = bt_sigma_pl;
Bi_ot = bi_sigma_ot;
Bi_pl = bi_sigma_pl;

% ---------------------------------------------------------------------
% 4. Equal target/ITI beta_SigmaV second-order matrices per day.
% ---------------------------------------------------------------------
[Ct_ot,cov_t_ot] = beta_time_covariance_local( ...
    Bt_ot,cfg.beta_covariance_mode,'OT target beta_SigmaV');
[Ci_ot,cov_i_ot] = beta_time_covariance_local( ...
    Bi_ot,cfg.beta_covariance_mode,'OT ITI beta_SigmaV');
[Ct_pl,cov_t_pl] = beta_time_covariance_local( ...
    Bt_pl,cfg.beta_covariance_mode,'PL target beta_SigmaV');
[Ci_pl,cov_i_pl] = beta_time_covariance_local( ...
    Bi_pl,cfg.beta_covariance_mode,'PL ITI beta_SigmaV');
C_ot = real((Ct_ot+Ci_ot)/2);
C_pl = real((Ct_pl+Ci_pl)/2);
[C_joint,reg_info] = regularize_joint_local(C_ot,C_pl,cfg);
invC = pinv(C_joint);

% ---------------------------------------------------------------------
% 5. Rank-2 joint LCMV at the requested source grids.
% ---------------------------------------------------------------------
n_ot = numel(ot_data.label);
n_pl = numel(pl_data.label);
n_request = numel(grid_idx);
W_ot_raw = nan(n_request,2,n_ot,cfg.filter_storage_class);
W_pl_raw = nan(n_request,2,n_pl,cfg.filter_storage_class);
W_ot_ung = nan(n_request,2,n_ot,cfg.filter_storage_class);
W_pl_ung = nan(n_request,2,n_pl,cfg.filter_storage_class);
Vr_all = nan(n_request,3,2,'single');
valid = false(n_request,1);
gain_rcond_ot = nan(n_request,1);
gain_rcond_pl = nan(n_request,1);

for k = 1:n_request
    g = grid_idx(k);
    Lot0 = ot_lf.leadfield{g};
    Lpl0 = pl_lf.leadfield{g};
    if isempty(Lot0) || isempty(Lpl0),continue;end
    Lot = double(Lot0(lf_idx_ot,:));
    Lpl = double(Lpl0(lf_idx_pl,:));
    Ljoint = [Lot;Lpl];
    if size(Ljoint,2)~=3 || ~all(isfinite(Ljoint(:))) || ...
            rank(Ljoint)<2,continue;end
    [~,~,V] = svd(Ljoint,'econ');
    Vr = V(:,1:2);
    Lrot = Lot*Vr;
    Lrpl = Lpl*Vr;
    Lr = [Lrot;Lrpl];
    G = real((Lr'*invC*Lr+(Lr'*invC*Lr)')/2);
    if rcond(G)<cfg.min_rcond,continue;end
    Wjoint = pinv(G)*(Lr'*invC);
    Wot = real(Wjoint(:,1:n_ot));
    Wpl = real(Wjoint(:,n_ot+1:end));

    Got = Wot*Lrot;
    Gpl = Wpl*Lrpl;
    gain_rcond_ot(k) = rcond(Got);
    gain_rcond_pl(k) = rcond(Gpl);
    if cfg.apply_day_gain_calibration
        if gain_rcond_ot(k)<cfg.min_rcond || ...
                gain_rcond_pl(k)<cfg.min_rcond,continue;end
        Wot = pinv(Got)*Wot;
        Wpl = pinv(Gpl)*Wpl;
    end

    % True two-component unit white-noise-gain normalization.  The beta
    % domain has no independent sensor noise covariance, so do not relabel
    % empirical ITI whitening as unit-noise-gain.
    Not = real((Wot*Wot'+(Wot*Wot')')/2);
    Npl = real((Wpl*Wpl'+(Wpl*Wpl')')/2);
    [Not_invhalf,ok_ot] = invsqrt_psd_local(Not,cfg.min_rcond);
    [Npl_invhalf,ok_pl] = invsqrt_psd_local(Npl,cfg.min_rcond);
    if ~ok_ot || ~ok_pl,continue;end
    Wot_ung = Not_invhalf*Wot;
    Wpl_ung = Npl_invhalf*Wpl;

    W_ot_raw(k,:,:) = reshape(cast(Wot,cfg.filter_storage_class),1,2,n_ot);
    W_pl_raw(k,:,:) = reshape(cast(Wpl,cfg.filter_storage_class),1,2,n_pl);
    W_ot_ung(k,:,:) = reshape(cast(Wot_ung,cfg.filter_storage_class),1,2,n_ot);
    W_pl_ung(k,:,:) = reshape(cast(Wpl_ung,cfg.filter_storage_class),1,2,n_pl);
    Vr_all(k,:,:) = reshape(single(Vr),1,3,2);
    valid(k) = true;
end
assert(any(valid),'No valid rank-2 ROI filter was constructed.');

bank = struct();
bank.grid_index = grid_idx(valid);
bank.W_ot_raw = W_ot_raw(valid,:,:);
bank.W_pl_raw = W_pl_raw(valid,:,:);
bank.W_ot_unitnoisegain = W_ot_ung(valid,:,:);
bank.W_pl_unitnoisegain = W_pl_ung(valid,:,:);
bank.orientation_basis_xyz = Vr_all(valid,:,:);
bank.label_ot = ot_data.label;
bank.label_pl = pl_data.label;
bank.normalization = 'two-component unit white-noise gain: (W*W'')^(-1/2)W';

info = struct('n_requested',n_request,'n_valid',sum(valid), ...
    'valid_requested_mask',valid,'gain_rcond_ot',gain_rcond_ot, ...
    'gain_rcond_pl',gain_rcond_pl,'valid_trial_ot',valid_ot, ...
    'valid_trial_pl',valid_pl,'predictor_scale',scale_info, ...
    'filter_target_time',time_t_ot,'filter_iti_time',time_i_ot, ...
    'glm_target_ot',glm_t_ot,'glm_target_pl',glm_t_pl, ...
    'glm_iti_ot',glm_i_ot,'glm_iti_pl',glm_i_pl, ...
    'sensor_beta',struct( ...
        'target_sigma_ot',bt_sigma_ot,'target_delta_ot',bt_delta_ot, ...
        'target_sigma_pl',bt_sigma_pl,'target_delta_pl',bt_delta_pl, ...
        'iti_sigma_ot',bi_sigma_ot,'iti_delta_ot',bi_delta_ot, ...
        'iti_sigma_pl',bi_sigma_pl,'iti_delta_pl',bi_delta_pl), ...
    'beta_filter_mode',cfg.beta_filter_mode, ...
    'beta_covariance_mode',cfg.beta_covariance_mode, ...
    'beta_covariance_target_sigma_ot',cov_t_ot, ...
    'beta_covariance_iti_sigma_ot',cov_i_ot, ...
    'beta_covariance_target_sigma_pl',cov_t_pl, ...
    'beta_covariance_iti_sigma_pl',cov_i_pl, ...
    'C_beta_target_ot',Ct_ot,'C_beta_target_pl',Ct_pl, ...
    'C_beta_iti_ot',Ci_ot,'C_beta_iti_pl',Ci_pl, ...
    'C_beta_ot',C_ot,'C_beta_pl',C_pl, ...
    'regularization',reg_info, ...
    'definition',['Y(t)~SigmaV+DeltaV; target/ITI time-resolved ' ...
    'beta_SigmaV -> equal target/ITI second-order matrix/day -> ' ...
    'rank-2 joint LCMV; separate day matrix unit gain; true ' ...
    'unit-white-noise-gain output retained']);
end


function [roi_amp,info] = reconstruct_vector_roi_local( ...
    ot_data,pl_data,bhv_ot,bhv_pl,bank,roi_defs,cfg)
[ot_crop,crop_ot] = crop_raw_local(ot_data,cfg.reconstruction_window);
[pl_crop,crop_pl] = crop_raw_local(pl_data,cfg.reconstruction_window);
assert(numel(ot_crop.trial)==height(bhv_ot), ...
    'OT reconstruction crop changed the MEG/behavior trial count.');
assert(numel(pl_crop.trial)==height(bhv_pl), ...
    'PL reconstruction crop changed the MEG/behavior trial count.');

fs_ot = get_fsample_local(ot_data);
fs_pl = get_fsample_local(pl_data);
[Yi_ot,~,iti_ot] = extract_iti_cube_local( ...
    ot_data,bhv_ot,cfg.baseline_iti_duration,cfg.baseline_iti_delay, ...
    cfg.baseline_min_iti_duration,fs_ot,'OT amplitude baseline');
[Yi_pl,~,iti_pl] = extract_iti_cube_local( ...
    pl_data,bhv_pl,cfg.baseline_iti_duration,cfg.baseline_iti_delay, ...
    cfg.baseline_min_iti_duration,fs_pl,'PL amplitude baseline');

[raw_ot,raw_base_ot,roi_count,time_ot,valid_time_ot] = ...
    project_day_roi_local( ...
    ot_crop,Yi_ot,iti_ot.valid,bank.W_ot_raw,bank.grid_index, ...
    roi_defs,cfg);
[raw_pl,raw_base_pl,~,time_pl,valid_time_pl] = ...
    project_day_roi_local( ...
    pl_crop,Yi_pl,iti_pl.valid,bank.W_pl_raw,bank.grid_index, ...
    roi_defs,cfg);
[ung_ot,ung_base_ot,~,time_ot_ung] = project_day_roi_local( ...
    ot_crop,Yi_ot,iti_ot.valid,bank.W_ot_unitnoisegain, ...
    bank.grid_index,roi_defs,cfg);
[ung_pl,ung_base_pl,~,time_pl_ung] = project_day_roi_local( ...
    pl_crop,Yi_pl,iti_pl.valid,bank.W_pl_unitnoisegain, ...
    bank.grid_index,roi_defs,cfg);
assert_same_time_local(time_ot,time_pl,'OT/PL reconstruction');
assert_same_time_local(time_ot,time_ot_ung,'OT raw/UNG reconstruction');
assert_same_time_local(time_pl,time_pl_ung,'PL raw/UNG reconstruction');
assert(numel(raw_ot)==height(bhv_ot) && ...
    numel(ung_ot)==height(bhv_ot), ...
    'OT projected trial count differs from behavior.');
assert(numel(raw_pl)==height(bhv_pl) && ...
    numel(ung_pl)==height(bhv_pl), ...
    'PL projected trial count differs from behavior.');

bhv_ot.day = ones(height(bhv_ot),1);
bhv_pl.day = 2*ones(height(bhv_pl),1);
behavior = [bhv_ot;bhv_pl];
raw_cells = [raw_ot;raw_pl];
ung_cells = [ung_ot;ung_pl];
baseline_raw = [raw_base_ot;raw_base_pl];
baseline_ung = [ung_base_ot;ung_base_pl];
labels = {roi_defs.label}';

[raw_no_base,time_new,fs_new] = resample_roi_cells_local( ...
    raw_cells,time_ot,labels,fs_ot,cfg.fs_new);
[ung_no_base,time_new_ung,fs_new_ung] = resample_roi_cells_local( ...
    ung_cells,time_ot,labels,fs_ot,cfg.fs_new);
assert_same_time_local(time_new,time_new_ung,'raw/UNG resampling');
assert(abs(fs_new-fs_new_ung)<1e-8);
assert(size(raw_no_base,1)==height(behavior) && ...
    size(ung_no_base,1)==height(behavior), ...
    'Resampled ROI data/behavior trial counts differ.');

% Baseline correction is performed AFTER source projection and the
% nonlinear vector-amplitude operation.  Explicit reshape is required:
% task amplitude is trial x ROI x time, whereas the trial-matched ITI
% amplitude mean is trial x ROI.
raw_corrected = raw_no_base - ...
    reshape(baseline_raw,size(baseline_raw,1),size(baseline_raw,2),1);
ung_corrected = ung_no_base - ...
    reshape(baseline_ung,size(baseline_ung,1),size(baseline_ung,2),1);

roi_amp = struct();
roi_amp.label = labels;
roi_amp.time = time_new;
roi_amp.fsample = fs_new;
roi_amp.behavior = behavior;
roi_amp.raw_no_baseline = raw_no_base;
roi_amp.unitnoisegain_no_baseline = ung_no_base;
roi_amp.raw_iti_baseline = baseline_raw;
roi_amp.unitnoisegain_iti_baseline = baseline_ung;
roi_amp.raw_iti_corrected = raw_corrected;
roi_amp.unitnoisegain_iti_corrected = ung_corrected;
roi_amp.iti_valid = [iti_ot.valid;iti_pl.valid];
roi_amp.dimord = 'rpt_chan_time';
roi_amp.filter_training_method = ...
    ['Y(t)~SigmaV+DeltaV; equal target/ITI time-resolved ' ...
     'beta_SigmaV second-order matrices'];
roi_amp.filter_beta_filter_mode = cfg.beta_filter_mode;
roi_amp.filter_beta_covariance_mode = cfg.beta_covariance_mode;
roi_amp.definition = ['SigmaV-beta-informed rank-2 filter applied to raw ' ...
    'trials; voxel-wise vector magnitude, then ROI arithmetic average; ' ...
    'no-baseline and trial-matched ITI-corrected outputs'];

info = struct('n_roi',numel(roi_defs),'n_voxel_by_roi',roi_count, ...
    'n_trial_ot',numel(raw_ot),'n_trial_pl',numel(raw_pl), ...
    'n_valid_iti_ot',sum(iti_ot.valid), ...
    'n_valid_iti_pl',sum(iti_pl.valid), ...
    'n_truncated_reconstruction_trial_ot',sum(~all(valid_time_ot,2)), ...
    'n_truncated_reconstruction_trial_pl',sum(~all(valid_time_pl,2)), ...
    'n_valid_reconstruction_trial_by_time_ot',sum(valid_time_ot,1), ...
    'n_valid_reconstruction_trial_by_time_pl',sum(valid_time_pl,1), ...
    'crop_ot',crop_ot,'crop_pl',crop_pl, ...
    'reconstruction_time_original',time_ot, ...
    'fs_original',fs_ot,'fs_saved',fs_new);
end


function [trial_amp,baseline,roi_count,time_ref,valid_time] = ...
    project_day_roi_local( ...
    data_crop,Yiti,iti_valid,Wbank,bank_grid,roi_defs,cfg)
n_trial = numel(data_crop.trial);
n_roi = numel(roi_defs);
fs = get_fsample_local(data_crop);
sample_ref = round(cfg.reconstruction_window(1)*fs): ...
    round(cfg.reconstruction_window(2)*fs);
time_ref = double(sample_ref/fs);
n_time = numel(time_ref);
n_bank = numel(bank_grid);
membership = false(n_bank,n_roi);
for r = 1:n_roi
    membership(:,r) = ismember(bank_grid,roi_defs(r).grid_idx);
end
roi_count = sum(membership,1);
assert(all(roi_count>0),'At least one requested ROI has no valid filter.');
trial_amp = cell(n_trial,1);
baseline = nan(n_trial,n_roi);
valid_time = false(n_trial,n_time);

for i_trial = 1:n_trial
    X0 = double(data_crop.trial{i_trial});
    time0 = double(data_crop.time{i_trial}(:)');
    assert(size(X0,2)==numel(time0), ...
        'Trial %d data/time lengths differ.',i_trial);
    assert(size(X0,1)==size(Wbank,3), ...
        'Trial %d channel count does not match the filter.',i_trial);

    % Align every trial to the requested common time axis. A truncated
    % trial is retained, but unavailable samples remain NaN. Alignment is
    % based on time stamps rather than assuming that missing samples are
    % always located at the end.
    ref_idx = round((time0-time_ref(1))*fs)+1;
    in_ref = ref_idx>=1 & ref_idx<=n_time;
    X = nan(size(X0,1),n_time);
    if any(in_ref)
        alignment_error = abs(time0(in_ref)-time_ref(ref_idx(in_ref)));
        assert(all(alignment_error<=max(1e-8,0.25/fs)), ...
            'Trial %d time samples are off the common sampling grid.', ...
            i_trial);
        assert(numel(unique(ref_idx(in_ref)))==sum(in_ref), ...
            'Trial %d has duplicate samples on the common time grid.', ...
            i_trial);
        X(:,ref_idx(in_ref)) = X0(:,in_ref);
    end
    valid_time(i_trial,:) = all(isfinite(X),1);

    has_iti = logical(iti_valid(i_trial));
    if has_iti
        Xi = reshape(double(Yiti(i_trial,:,:)),size(Yiti,2),[]);
        Xi = Xi(:,all(isfinite(Xi),1));
        has_iti = ~isempty(Xi);
    else
        Xi = [];
    end
    sum_task = zeros(n_roi,n_time);
    n_task = zeros(n_roi,n_time);
    sum_base = zeros(1,n_roi);
    n_base = zeros(1,n_roi);
    for k = 1:n_bank
        W = reshape(double(Wbank(k,:,:)),2,[]);
        Qt = W*X;
        if has_iti,Qi=W*Xi;else,Qi=[];end
        switch lower(cfg.amplitude_measure)
            case 'magnitude'
                at = sqrt(sum(Qt.^2,1));
                if has_iti,ai=sqrt(sum(Qi.^2,1));else,ai=[];end
            case 'power'
                at = sum(Qt.^2,1);
                if has_iti,ai=sum(Qi.^2,1);else,ai=[];end
            otherwise
                error('Unknown amplitude_measure: %s',cfg.amplitude_measure);
        end
        if cfg.log_transform
            pooled = [at(isfinite(at)),ai(isfinite(ai))];
            if isempty(pooled)
                offset = realmin('double');
            else
                offset = max( ...
                    median(pooled,'omitnan')*cfg.log_offset_ratio, ...
                    realmin('double'));
            end
            at = log(at+offset);
            if has_iti,ai=log(ai+offset);end
        end
        for r = find(membership(k,:))
            valid_at = isfinite(at);
            sum_task(r,valid_at) = sum_task(r,valid_at)+at(valid_at);
            n_task(r,valid_at) = n_task(r,valid_at)+1;
            if has_iti
                base_k = mean(ai,'omitnan');
                if isfinite(base_k)
                    sum_base(r) = sum_base(r)+base_k;
                    n_base(r) = n_base(r)+1;
                end
            end
        end
    end
    task = sum_task./n_task;
    task(n_task==0) = NaN;
    trial_amp{i_trial} = task;
    valid_base_roi = n_base>0;
    baseline(i_trial,valid_base_roi) = ...
        sum_base(valid_base_roi)./n_base(valid_base_roi);
end
end


function [array,time_new,fs_new] = resample_roi_cells_local( ...
    trial_cells,time,labels,fs,requested_fs)
time = double(time(:)');
n_trial = numel(trial_cells);
n_roi = numel(labels);

if isempty(requested_fs) || abs(requested_fs-fs)<1e-8
    fs_new = fs;
    time_new = time;
    array = nan(n_trial,n_roi,numel(time_new));
    for i = 1:n_trial
        Xi = double(trial_cells{i});
        assert(isequal(size(Xi),[n_roi,numel(time_new)]));
        array(i,:,:) = Xi;
    end
    return;
end

% Do not send NaN-padded trials directly to ft_resampledata: filtering can
% spread an unavailable tail into otherwise valid samples. Build one
% variable-length FieldTrip dataset from the finite contiguous segment of
% every trial, resample all retained trials together, and only then place
% them on the exact common output grid with NaN padding.
fs_new = double(requested_fs);
sample_new = round(time(1)*fs_new):round(time(end)*fs_new);
time_new = double(sample_new/fs_new);
array = nan(n_trial,n_roi,numel(time_new));

v0 = struct();
v0.label = labels(:);
v0.fsample = fs;
v0.time = {};
v0.trial = {};
v0.dimord = 'rpt_chan_time';
original_trial_index = zeros(0,1);
valid_start = zeros(0,1);
valid_end = zeros(0,1);

for i = 1:n_trial
    Xi = double(trial_cells{i});
    assert(isequal(size(Xi),[n_roi,numel(time)]));
    valid_col = all(isfinite(Xi),1);
    if ~any(valid_col),continue;end
    first_valid = find(valid_col,1,'first');
    last_valid = find(valid_col,1,'last');
    assert(all(valid_col(first_valid:last_valid)), ...
        'Trial %d has an internal NaN gap; cannot safely resample.',i);

    j = numel(original_trial_index)+1;
    original_trial_index(j,1) = i;
    valid_start(j,1) = time(first_valid);
    valid_end(j,1) = time(last_valid);
    v0.time{j} = time(first_valid:last_valid);
    v0.trial{j} = Xi(:,first_valid:last_valid);
end

if isempty(original_trial_index)
    return;
end

cfgr = [];
cfgr.resamplefs = fs_new;
cfgr.method = 'resample';
cfgr.detrend = 'no';
cfgr.demean = 'no';
cfgr.feedback = 'no';
v = ft_resampledata(cfgr,v0);

for j = 1:numel(original_trial_index)
    i = original_trial_index(j);
    ti = double(v.time{j}(:)');
    Yi = double(v.trial{j});
    assert(size(Yi,1)==n_roi && size(Yi,2)==numel(ti), ...
        'Resampled trial %d data/time dimensions differ.',i);
    assert(numel(ti)>=2 && all(diff(ti)>0), ...
        'Resampled trial %d has an invalid time axis.',i);

    % ft_resampledata shifts its output by a sub-sample amount so that the
    % centres of the old and new time axes agree. The shift can approach
    % half a 100-Hz sample (5 ms), and can depend on trial length. Rather
    % than merely relabelling the shifted samples, interpolate the already
    % anti-aliased signal onto the exact common event-locked grid.
    tol_old = max(1e-8,0.25/fs);
    out_idx = find(time_new>=valid_start(j)-tol_old & ...
        time_new<=valid_end(j)+tol_old);
    if isempty(out_idx),continue;end
    target_time = time_new(out_idx);
    edge_gap = max([0,ti(1)-target_time(1), ...
        target_time(end)-ti(end)]);
    assert(edge_gap<=0.51/fs_new+tol_old, ...
        ['Resampled trial %d is more than half a new sample away ' ...
        'from the common output grid.'],i);
    Yi_exact = interp1(ti(:),Yi.',target_time(:), ...
        'linear','extrap').';
    array(i,:,out_idx) = reshape( ...
        Yi_exact,1,n_roi,numel(out_idx));
end
end


function v = array_to_ft_raw_local(array,time,labels,bhv,fs)
n_trial = size(array,1);
assert(n_trial==height(bhv), ...
    'ROI array has %d trials but behavior has %d rows.', ...
    n_trial,height(bhv));
v = struct();
v.label = labels(:);
v.fsample = fs;
v.time = repmat({double(time(:)')},1,n_trial);
v.trial = cell(1,n_trial);
for i = 1:n_trial
    v.trial{i} = reshape(double(array(i,:,:)),size(array,2),size(array,3));
end
v.dimord = 'rpt_chan_time';
v.trialinfo_label = bhv.Properties.VariableNames;
v.trialinfo = table2array(bhv);
end


function [sigma_beta,delta_beta,info] = fit_time_glm_local( ...
    Y,sigma_z,delta_z,bhv,valid,cfg,name)
% Fit all channels and all time points in one matrix solve. Y remains in
% physical sensor units. SigmaV and DeltaV use pooled OT/PL scaling.
valid = logical(valid(:));
n = sum(valid);
D = [ones(n,1),sigma_z(valid),delta_z(valid)];
column_names = {'intercept','SigmaV','DeltaV'};
assert(rank(D)==size(D,2), ...
    '%s primary GLM design is rank deficient.',name);

candidate = zeros(n,0);
candidate_names = {};
if cfg.add_change_nuisance
    candidate = [candidate,double(bhv.ischange(valid))];
    candidate_names{end+1} = 'ischange';
end
if cfg.add_condition_nuisance
    [D0,N0] = make_dummies_local( ...
        double(bhv.condition(valid)),'condition');
    candidate = [candidate,D0];
    candidate_names = [candidate_names,N0];
end
if cfg.add_run_nuisance
    assert(ismember('run_uid',bhv.Properties.VariableNames), ...
        '%s behavior lacks run_uid.',name);
    [D0,N0] = make_dummies_local(double(bhv.run_uid(valid)),'run_uid');
    candidate = [candidate,D0];
    candidate_names = [candidate_names,N0];
end

dropped = {};
for i = 1:size(candidate,2)
    if rank([D,candidate(:,i)])>rank(D)
        D = [D,candidate(:,i)]; %#ok<AGROW>
        column_names{end+1} = candidate_names{i}; %#ok<AGROW>
    else
        dropped{end+1} = candidate_names{i}; %#ok<AGROW>
    end
end

n_chan = size(Y,2);
n_time = size(Y,3);
Y2 = reshape(Y(valid,:,:),n,n_chan*n_time);
assert(all(isfinite(Y2(:))),'%s contains non-finite GLM Y.',name);
B = D\Y2;
sigma_beta = reshape(B(2,:),n_chan,n_time);
delta_beta = reshape(B(3,:),n_chan,n_time);
info = struct('name',name,'n_trial',n, ...
    'design_columns',{column_names}, ...
    'dropped_redundant_nuisance',{dropped}, ...
    'rank',rank(D),'n_column',size(D,2), ...
    'X_standardized',true,'Y_standardized',false, ...
    'primary_beta_for_filter','SigmaV');
end


function [D,names] = make_dummies_local(code,prefix)
levels = unique(code(:),'stable');
D = zeros(numel(code),max(numel(levels)-1,0));
names = cell(1,size(D,2));
for i = 2:numel(levels)
    D(:,i-1) = code==levels(i);
    names{i-1} = sprintf('%s_%g',prefix,levels(i));
end
end


function [C,info] = beta_time_covariance_local(B,mode,name)
% B is channel x beta-time. In secondmoment mode, a sustained WTA beta is
% retained. Covariance mode removes each channel's temporal beta mean.
assert(size(B,2)>=2,'%s has fewer than two beta time points.',name);
assert(all(isfinite(B(:))),'%s contains non-finite beta values.',name);
switch lower(mode)
    case 'covariance'
        B_used = B-mean(B,2);
        denominator = size(B,2)-1;
    case 'secondmoment'
        B_used = B;
        denominator = size(B,2);
    otherwise
        error('Unknown cfg.beta_covariance_mode: %s',mode);
end
C = (B_used*B_used')/denominator;
C = real((C+C')/2);
assert(all(isfinite(C(:))) && trace(C)>0 && norm(C,'fro')>0, ...
    '%s beta covariance/second moment is invalid.',name);
info = struct('name',name,'mode',lower(mode), ...
    'n_beta_time',size(B,2),'rank',rank(C), ...
    'trace',trace(C),'rcond',rcond(C));
end


function [C,info] = regularize_joint_local(Cot,Cpl,cfg)
switch lower(cfg.regularization_mode)
    case 'joint_global'
        C0 = blkdiag(Cot,Cpl);
        scale = trace(C0)/size(C0,1);
        assert(isfinite(scale) && scale>0);
        C = C0+cfg.lambda*scale*eye(size(C0));
        info = struct('mode','joint_global','lambda',cfg.lambda, ...
            'scale_joint',scale);
    case 'per_day'
        sot = trace(Cot)/size(Cot,1);
        spl = trace(Cpl)/size(Cpl,1);
        assert(all(isfinite([sot,spl])) && sot>0 && spl>0);
        C = blkdiag(Cot+cfg.lambda*sot*eye(size(Cot)), ...
            Cpl+cfg.lambda*spl*eye(size(Cpl)));
        info = struct('mode','per_day','lambda',cfg.lambda, ...
            'scale_ot',sot,'scale_pl',spl);
    otherwise
        error('Unknown regularization_mode: %s',cfg.regularization_mode);
end
C = real((C+C')/2);
info.rcond = rcond(C);
end


function [A,ok] = invsqrt_psd_local(M,min_rcond)
M = real((M+M')/2);
[V,D] = eig(M,'vector');
D = real(D);
dmax = max(D);
ok = all(isfinite(D)) && dmax>0 && min(D)>dmax*min_rcond;
if ok
    A = V*diag(1./sqrt(D))*V';
    A = real((A+A')/2);
else
    A = nan(size(M));
end
end


function [Y,time_ref,info] = extract_fixed_cube_local(data,window,name)
n_trial = numel(data.trial);
n_chan = numel(data.label);
t0 = double(data.time{1}(:)');
idx0 = find(t0>=window(1) & t0<window(2));
assert(numel(idx0)>=2,'%s contains fewer than two samples.',name);
time_ref = t0(idx0);
Y = nan(n_trial,n_chan,numel(idx0));
valid = false(n_trial,1);
for i_trial = 1:n_trial
    t = double(data.time{i_trial}(:)');
    idx = find(t>=window(1) & t<window(2));
    if numel(idx)~=numel(idx0) || ...
            max(abs(t(idx)-time_ref))>1e-8,continue;end
    X = double(data.trial{i_trial}(:,idx));
    if ~all(isfinite(X(:))),continue;end
    Y(i_trial,:,:) = X;
    valid(i_trial) = true;
end
info = struct('name',name,'window',window,'valid',valid, ...
    'n_valid',sum(valid),'n_invalid',sum(~valid));
end


function [Y,relative_time,info] = extract_iti_cube_local( ...
    data,bhv,duration,delay,min_duration,fs,name)
n_trial = numel(data.trial);
n_chan = numel(data.label);
n_requested = round(duration*fs);
n_min = ceil(min_duration*fs-1e-8);
assert(n_requested>=2 && n_min>=2 && n_min<=n_requested);
relative_time = (0:n_requested-1)/fs;
Y = nan(n_trial,n_chan,n_requested);
valid = false(n_trial,1);
n_available = zeros(n_trial,1);

start_sec = double(bhv.oiti)-double(bhv.osti2)- ...
    double(bhv.rtime2)+delay;

for i_trial = 1:n_trial
    t = double(data.time{i_trial}(:)');
    if ~isfinite(start_sec(i_trial)) || isempty(t),continue;end
    idx_start = round((start_sec(i_trial)-t(1))*fs)+1;
    if idx_start<1 || idx_start>numel(t),continue;end
    n_take = min(n_requested,numel(t)-idx_start+1);
    if n_take<n_min,continue;end
    idx = idx_start+(0:n_take-1);
    X = double(data.trial{i_trial}(:,idx));
    if ~all(isfinite(X(:))),continue;end
    Y(i_trial,:,1:n_take) = X;
    valid(i_trial) = true;
    n_available(i_trial) = n_take;
end
info = struct('name',name,'duration',duration,'delay',delay, ...
    'min_duration',min_duration,'start_sec',start_sec, ...
    'n_available_sample',n_available,'valid',valid, ...
    'full_length',valid & n_available==n_requested, ...
    'n_valid',sum(valid),'n_partial',sum(valid & n_available<n_requested));
end


function [data,info] = crop_raw_local(data,window)
% Crop each trial without allowing FieldTrip to remove a non-overlapping
% trial. Keeping an empty trial placeholder is essential because behavior,
% trial-matched ITI and task data must retain exactly the same row order.
assert(numel(window)==2 && all(isfinite(window)) && window(1)<window(2));
n_trial = numel(data.trial);
fs = get_fsample_local(data);
n_requested = round((window(2)-window(1))*fs)+1;
n_available = zeros(n_trial,1);

for i_trial = 1:n_trial
    X = data.trial{i_trial};
    t = double(data.time{i_trial}(:)');
    assert(size(X,2)==numel(t), ...
        'Trial %d data/time lengths differ before cropping.',i_trial);
    tol = max(1e-8,0.25/fs);
    select = t>=window(1)-tol & t<=window(2)+tol;
    data.trial{i_trial} = X(:,select);
    data.time{i_trial} = t(select);
    n_available(i_trial) = sum(select);
end

% sampleinfo no longer describes a single contiguous equal-length epoch
% after the trial-wise crop and is not used in reconstruction.
data = rmfield_if_exists_local(data,{'sampleinfo'});
no_overlap = n_available==0;
partial = n_available>0 & n_available<n_requested;
info = struct( ...
    'window',double(window(:)'), ...
    'n_requested_sample',n_requested, ...
    'n_available_sample',n_available, ...
    'no_overlap_trial_index',find(no_overlap), ...
    'partial_trial_index',find(partial), ...
    'n_no_overlap',sum(no_overlap), ...
    'n_partial',sum(partial));
if any(no_overlap)
    warning(['%d trial(s) do not overlap reconstruction_window and are ' ...
        'retained as all-NaN rows to preserve behavior alignment.'], ...
        sum(no_overlap));
end
end


function [idx,info] = select_significant_cluster_local( ...
    stat_res,requested_id,alpha)
assert(isfield(stat_res,'posclusterslabelmat') && ...
    isfield(stat_res,'posclusters') && isfield(stat_res,'stat'));
prob = [stat_res.posclusters.prob];
if isempty(requested_id)
    candidate = find(isfinite(prob) & prob<=alpha);
    assert(~isempty(candidate),'No positive cluster survives p<=%.4f.',alpha);
    [~,j] = min(prob(candidate));
    cluster_id = candidate(j);
else
    cluster_id = double(requested_id);
    assert(isscalar(cluster_id) && cluster_id>=1 && ...
        cluster_id<=numel(prob) && cluster_id==round(cluster_id));
end
idx = find(double(stat_res.posclusterslabelmat(:))==cluster_id);
assert(~isempty(idx),'Cluster %d is empty.',cluster_id);
[peak,j] = max(double(stat_res.stat(idx)));
info = struct('cluster_id',cluster_id,'cluster_p',prob(cluster_id), ...
    'alpha',alpha,'n_grid',numel(idx),'grid_idx',idx, ...
    'peak_grid_idx',idx(j),'peak_stat',peak);
end


function [idx,info] = make_target_roi_conjunction_local( ...
    stat_res,cluster_idx,cfg)
assert(isfile(cfg.target_roi_file),'Cannot find %s',cfg.target_roi_file);
target = ft_read_mri(cfg.target_roi_file);
cfgi = [];
cfgi.parameter = cfg.target_roi_parameter;
cfgi.interpmethod = 'nearest';
target_grid = ft_sourceinterpolate(cfgi,target,stat_res);
value = double(target_grid.(cfg.target_roi_parameter)(:));
idx_target = find(isfinite(value) & value>cfg.target_roi_threshold);
idx = intersect(cluster_idx(:),idx_target(:),'stable');
assert(~isempty(idx),'Cluster has no intersection with target ROI.');
[peak,j] = max(double(stat_res.stat(idx)));
peak_idx = idx(j);
peak_label = 'Unlabeled';
if isfile(cfg.aal_file)
    aal = ft_read_atlas(cfg.aal_file);
    cfga = [];
    cfga.parameter = 'tissue';
    cfga.interpmethod = 'nearest';
    aal_grid = ft_sourceinterpolate(cfga,aal,stat_res);
    ti = double(aal_grid.tissue(peak_idx));
    if isfinite(ti) && ti>=1 && ti<=numel(aal.tissuelabel)
        peak_label = aal.tissuelabel{ti};
    end
end
info = struct('target_roi_file',cfg.target_roi_file, ...
    'target_roi_parameter',cfg.target_roi_parameter, ...
    'target_roi_threshold',cfg.target_roi_threshold, ...
    'n_target_grid',numel(idx_target), ...
    'n_conjunction_grid',numel(idx), ...
    'grid_idx',idx,'peak_grid_idx',peak_idx, ...
    'peak_stat',peak,'peak_aal_label',peak_label);
end


function [ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf,quality] = ...
    prepare_subject_local(i_sub,subject_id,data_bhv,head_mov, ...
    cfg_main,processing_window)
all_run_data = cell(12,1);
all_run_bhv = cell(12,1);
all_run_meta = cell(12,1);
run_counter = 0;
ot_slots = [1:3 7:9];
pl_slots = [4:6 10:12];
bad_ot = get_subject_index_vector_local( ...
    head_mov.isolated_bad_run.ot,i_sub);
bad_pl = get_subject_index_vector_local( ...
    head_mov.isolated_bad_run.pl,i_sub);
optimal_ot = get_subject_scalar_local(head_mov.optimal_idx.ot,i_sub);
optimal_pl = get_subject_scalar_local(head_mov.optimal_idx.pl,i_sub);
validate_day_run_indices_local(bad_ot,'OT isolated bad run');
validate_day_run_indices_local(bad_pl,'PL isolated bad run');
validate_day_run_indices_local(optimal_ot,'OT optimal run');
validate_day_run_indices_local(optimal_pl,'PL optimal run');
assert(~ismember(optimal_ot,bad_ot) && ~ismember(optimal_pl,bad_pl), ...
    'An optimal leadfield run is also marked as a bad run.');
optimal_ot_slot = ot_slots(optimal_ot);
optimal_pl_slot = pl_slots(optimal_pl);

for i_condition = 1:4
    raw_file = find_unique_subject_file_local( ...
        cfg_main.condition_dirs{i_condition},subject_id);
    data_process = load_raw_with_rejection_local(raw_file);
    bhv0 = data_bhv{i_condition}(i_sub);
    assert(strcmp(normalize_subject_id_local(bhv0.subsName),subject_id));
    [bhv0,bad_global,bad_by_run] = ...
        remove_rejected_behavior_trials_local( ...
        bhv0,data_process.rejection_info,cfg_main.trials_per_run);
    assert(numel(data_process.trial)==get_behavior_trial_count_local(bhv0));
    bhv0.run_id = make_run_id_after_rejection_local( ...
        cfg_main.trials_per_run,bad_global);
    n_bhv = get_behavior_trial_count_local(bhv0);
    bhv0.condition = repmat(i_condition,n_bhv,1);
    bhv0.original_trial_id = setdiff( ...
        (1:3*cfg_main.trials_per_run)',bad_global,'stable');
    bhv0.run_uid = (i_condition-1)*3+bhv0.run_id;
    data_process = rmfield_if_exists_local(data_process,{'rejection_info'});

    % Keep the full long epoch here because the matched pre-trial ITI can
    % lie outside reconstruction_window. The task segment is cropped only
    % inside reconstruct_vector_roi_local.
    % cfgc = [];
    % cfgc.toilim = processing_window;
    % data_process = ft_redefinetrial(cfgc,data_process);
    for i_run = 1:3
        run_counter = run_counter+1;
        run_idx = find(bhv0.run_id==i_run & double(bhv0.conflict)==1);
        all_run_meta{run_counter} = struct('condition',i_condition, ...
            'run',i_run,'raw_file',raw_file, ...
            'bad_trials_local',bad_by_run{i_run}, ...
            'n_conflict_trial',numel(run_idx));
        if isempty(run_idx),continue;end
        cfgs = [];
        cfgs.trials = run_idx;
        cfgs.channel = 'MEGGRAD';
        run_data = ft_selectdata(cfgs,data_process);
        run_data = rmfield_if_exists_local(run_data, ...
            {'grad','elec','hdr','sampleinfo','rejection_info'});
        run_bhv = behavior_struct_to_table_local( ...
            bhv0,run_idx,cfg_main.behavior_fields);
        assert(numel(run_data.trial)==height(run_bhv));
        all_run_data{run_counter} = run_data;
        all_run_bhv{run_counter} = run_bhv;
    end
end
assert(run_counter==12);

ot_keep = setdiff(1:6,bad_ot,'stable');
pl_keep = setdiff(1:6,bad_pl,'stable');
ot_used_slots = ot_slots(ot_keep);
pl_used_slots = pl_slots(pl_keep);
[ot_data,bhv_ot] = append_selected_slots_local( ...
    all_run_data,all_run_bhv,ot_used_slots,'OT');
[pl_data,bhv_pl] = append_selected_slots_local( ...
    all_run_data,all_run_bhv,pl_used_slots,'PL');
bhv_ot = ensure_wta_predictors_local(bhv_ot,'OT');
bhv_pl = ensure_wta_predictors_local(bhv_pl,'PL');

optimal_ot_condition = ceil(optimal_ot_slot/3);
optimal_pl_condition = ceil(optimal_pl_slot/3);
optimal_ot_run = mod(optimal_ot_slot-1,3)+1;
optimal_pl_run = mod(optimal_pl_slot-1,3)+1;
ot_lf_file = find_unique_subject_file_local( ...
    cfg_main.leadfield_dirs{optimal_ot_condition},subject_id);
pl_lf_file = find_unique_subject_file_local( ...
    cfg_main.leadfield_dirs{optimal_pl_condition},subject_id);
Sot = load(ot_lf_file,'leadfield_individual');
Spl = load(pl_lf_file,'leadfield_individual');
assert(iscell(Sot.leadfield_individual) && ...
    numel(Sot.leadfield_individual)==3);
assert(iscell(Spl.leadfield_individual) && ...
    numel(Spl.leadfield_individual)==3);
ot_lf = ft_convert_units( ...
    Sot.leadfield_individual{optimal_ot_run},'mm');
pl_lf = ft_convert_units( ...
    Spl.leadfield_individual{optimal_pl_run},'mm');

quality = struct('subject_id',subject_id, ...
    'removed_ot_run_indices',bad_ot,'removed_pl_run_indices',bad_pl, ...
    'ot_used_slots',ot_used_slots,'pl_used_slots',pl_used_slots, ...
    'optimal_ot_index_within_day',optimal_ot, ...
    'optimal_pl_index_within_day',optimal_pl, ...
    'optimal_ot_slot',optimal_ot_slot,'optimal_pl_slot',optimal_pl_slot, ...
    'optimal_ot_leadfield_file',ot_lf_file, ...
    'optimal_pl_leadfield_file',pl_lf_file, ...
    'optimal_ot_run_within_condition',optimal_ot_run, ...
    'optimal_pl_run_within_condition',optimal_pl_run, ...
    'run_metadata',{all_run_meta},'n_trial_ot',height(bhv_ot), ...
    'n_trial_pl',height(bhv_pl), ...
    'configured_analysis_window',processing_window);
end


function [data,labels,lf_idx] = align_data_lf_local( ...
    data,lf,selection,name)
cfgs = [];
cfgs.channel = selection;
data = ft_selectdata(cfgs,data);
shared = ismember(data.label,lf.label);
labels = data.label(shared);
assert(~isempty(labels),'%s has no channels shared with leadfield.',name);
cfgs = [];
cfgs.channel = labels;
data = ft_selectdata(cfgs,data);
[found,lf_idx] = ismember(data.label,lf.label);
assert(all(found),'%s data/leadfield alignment failed.',name);
labels = data.label(:);
end


function [data,bhv] = append_selected_slots_local( ...
    all_data,all_bhv,slots,name)
slots = double(slots(:)');
slots = slots(~cellfun(@isempty,all_data(slots)));
assert(~isempty(slots),'%s has no retained runs.',name);
data_cells = all_data(slots);
bhv_cells = all_bhv(slots);
reference_label = data_cells{1}.label(:);
for i = 1:numel(data_cells)
    this_label = data_cells{i}.label(:);
    assert(numel(this_label)==numel(reference_label) && ...
        all(ismember(reference_label,this_label)) && ...
        all(ismember(this_label,reference_label)), ...
        '%s retained runs have different channel sets.',name);
    if ~isequal(this_label,reference_label)
        cfgs = [];
        cfgs.channel = reference_label;
        data_cells{i} = ft_selectdata(cfgs,data_cells{i});
    end
end
cfga = [];
cfga.keepsampleinfo = 'no';
data = ft_appenddata(cfga,data_cells{:});
bhv = vertcat(bhv_cells{:});
assert(numel(data.trial)==height(bhv));
end


function file = find_unique_subject_file_local(folder,subject_id)
items = dir(fullfile(folder,[subject_id '*']));
items = items(~[items.isdir]);
assert(isscalar(items),'Expected one file for %s in %s; found %d.', ...
    subject_id,folder,numel(items));
file = fullfile(items(1).folder,items(1).name);
end


function data = load_raw_with_rejection_local(file)
S = load(file);
names = fieldnames(S);
data = [];
for i = 1:numel(names)
    x = S.(names{i});
    if isstruct(x) && isfield(x,'trial') && isfield(x,'time') && ...
            isfield(x,'label')
        assert(isempty(data),'Multiple raw structures in %s.',file);
        data = x;
    end
end
assert(~isempty(data),'No FieldTrip raw structure found in %s.',file);
assert(isfield(data,'rejection_info') && ...
    iscell(data.rejection_info) && numel(data.rejection_info)==3, ...
    'Missing three-cell rejection_info in %s.',file);
end


function [bhv,bad_global,bad_by_run] = ...
    remove_rejected_behavior_trials_local( ...
    bhv,rejection_info,trials_per_run)
n_original = get_behavior_trial_count_local(bhv);
assert(n_original==3*trials_per_run);
bad_global = [];
bad_by_run = cell(3,1);
for i_run = 1:3
    b = unique(double(rejection_info{i_run}.bad_trials(:)),'stable');
    assert(all(b>=1 & b<=trials_per_run));
    bad_by_run{i_run} = b;
    bad_global = [bad_global;b+(i_run-1)*trials_per_run]; %#ok<AGROW>
end
bad_global = unique(bad_global,'stable');
fields = fieldnames(bhv);
for i = 1:numel(fields)
    value = bhv.(fields{i});
    if (isnumeric(value) || islogical(value) || iscell(value) || ...
            isstring(value)) && isvector(value) && numel(value)==n_original
        value(bad_global) = [];
        bhv.(fields{i}) = value;
    end
end
end


function n = get_behavior_trial_count_local(bhv)
candidates = {'answer1','conflict','ischange','lo_s1'};
n = [];
for i = 1:numel(candidates)
    if isfield(bhv,candidates{i})
        n = numel(bhv.(candidates{i}));
        break;
    end
end
assert(~isempty(n),'Cannot determine behavior trial count.');
end


function run_id = make_run_id_after_rejection_local( ...
    trials_per_run,bad_global)
run_id = repelem((1:3)',trials_per_run);
run_id(bad_global) = [];
end


function tbl = behavior_struct_to_table_local(bhv,idx,fields)
tbl = table();
for i = 1:numel(fields)
    f = fields{i};
    if ~isfield(bhv,f)
        if ismember(f,{'lo_chosen','lo_unchosen'}),continue;end
        error('Behavior lacks required field %s.',f);
    end
    value = bhv.(f);
    assert(numel(value)==get_behavior_trial_count_local(bhv));
    tbl.(f) = double(value(idx));
    tbl.(f) = tbl.(f)(:);
end
end


function tbl = ensure_wta_predictors_local(tbl,name)
required = {'conflict','ischange','lo_s1','lo_o1'};
assert(all(ismember(required,tbl.Properties.VariableNames)), ...
    '%s behavior lacks fields required to derive chosen/unchosen.',name);
assert(all(tbl.conflict==1), ...
    '%s reconstruction data unexpectedly contain congruent trials.',name);
stay = tbl.ischange==0;
change = tbl.ischange==1;
assert(all(stay | change),'%s contains invalid ischange values.',name);
chosen = nan(height(tbl),1);
unchosen = nan(height(tbl),1);
chosen(stay) = tbl.lo_s1(stay);
unchosen(stay) = tbl.lo_o1(stay);
chosen(change) = tbl.lo_o1(change);
unchosen(change) = tbl.lo_s1(change);
if ismember('lo_chosen',tbl.Properties.VariableNames)
    delta = abs(tbl.lo_chosen-chosen);
    assert(all(delta<1e-10 | ~isfinite(delta)), ...
        '%s existing lo_chosen is inconsistent with ischange.',name);
else
    tbl.lo_chosen = chosen;
end
if ismember('lo_unchosen',tbl.Properties.VariableNames)
    delta = abs(tbl.lo_unchosen-unchosen);
    assert(all(delta<1e-10 | ~isfinite(delta)), ...
        '%s existing lo_unchosen is inconsistent with ischange.',name);
else
    tbl.lo_unchosen = unchosen;
end
end


function data = rmfield_if_exists_local(data,names)
for i = 1:numel(names)
    if isfield(data,names{i}),data=rmfield(data,names{i});end
end
end


function subject_ids = get_and_verify_subject_ids_local(data_bhv)
n_set = numel(data_bhv);
assert(n_set==4,'Expected four behavior condition sets.');
n_subject = numel(data_bhv{1});
subject_ids = cell(n_subject,1);
for i_sub = 1:n_subject
    subject_ids{i_sub} = normalize_subject_id_local( ...
        data_bhv{1}(i_sub).subsName);
end
assert(numel(unique(subject_ids))==n_subject, ...
    'Behavior subject IDs are not unique.');
for i_set = 2:n_set
    assert(numel(data_bhv{i_set})==n_subject, ...
        'Behavior condition sets contain different subject counts.');
    for i_sub = 1:n_subject
        this_id = normalize_subject_id_local( ...
            data_bhv{i_set}(i_sub).subsName);
        assert(strcmp(this_id,subject_ids{i_sub}), ...
            'Behavior subject order differs in condition set %d.',i_set);
    end
end
end


function id = normalize_subject_id_local(value)
if iscell(value),value=value{1};end
id = strtrim(char(string(value)));
end


function value = get_subject_index_vector_local(container,i_sub)
if iscell(container),value=container{i_sub};else,value=container(i_sub,:);end
value = double(value(:)');
value = value(isfinite(value) & value~=0);
end


function value = get_subject_scalar_local(container,i_sub)
if iscell(container),value=container{i_sub};else,value=container(i_sub);end
value = double(value);
assert(isscalar(value) && isfinite(value));
end


function validate_day_run_indices_local(idx,name)
assert(all(idx==round(idx) & idx>=1 & idx<=6), ...
    '%s must contain integers from 1 to 6.',name);
end


function grid = load_source_grid_local(file)
S = load(file);
names = fieldnames(S);
grid = [];
for i = 1:numel(names)
    x = S.(names{i});
    if isstruct(x) && isfield(x,'pos')
        grid = x;
        break;
    end
end
assert(~isempty(grid),'No source grid found in %s.',file);
end


function verify_stat_source_geometry_local(stat_res,source)
assert(numel(stat_res.stat)==size(source.pos,1), ...
    'stat_res and source grid have different point counts.');
if isfield(stat_res,'pos')
    p1 = double(stat_res.pos);
    p2 = double(source.pos);
    d0 = max(abs(p1(:)-p2(:)),[],'omitnan');
    d10 = max(abs(10*p1(:)-p2(:)),[],'omitnan');
    d01 = max(abs(p1(:)-10*p2(:)),[],'omitnan');
    assert(min([d0,d10,d01])<1e-4, ...
        'stat_res and source-grid positions/order differ.');
end
end


function fs = get_fsample_local(data)
if isfield(data,'fsample') && isfinite(data.fsample)
    fs = double(data.fsample);
else
    fs = 1/median(diff(double(data.time{1})));
end
end


function assert_same_time_local(a,b,name)
a = double(a(:)');
b = double(b(:)');
assert(numel(a)==numel(b) && max(abs(a-b))<1e-8, ...
    '%s time axes differ.',name);
end






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross test part2: conjunction map(coefficient(ΣV) source ∩ oxytocin receptor expression map) → 
% gamma-band source filter → single reconstraction → time-frequency
clear;
%% Configuration
cfg = struct();
cfg.n_subject = 19;
cfg.filter_time = [-0.5 -0.3];
cfg.tfr_time_high = [-1.2 0.8];
cfg.tfr_time_low  = [-1.8 1.3];
% Final time range retained in the saved TFR structures.
cfg.output_tfr_time = [-1 0.5];
cfg.high_fields_contain_merged_low_high = true;
cfg.tfr_time = [min(cfg.tfr_time_high(1),cfg.tfr_time_low(1)), ...
                max(cfg.tfr_time_high(2),cfg.tfr_time_low(2))];
% Assigned below after the frequency-dependent low-frequency windows are
% defined. The raw trials must cover the full convolution support.
cfg.long_time = [];
cfg.expected_fsample = 250;
cfg.target_foi = [30 100];
cfg.lambda = 0.10;
cfg.gain_pinv_rtol = 1e-10;

% Match the uploaded whole-brain DICS source code.
cfg.apply_day_unit_gain = false;
cfg.grid_ung_for_source_svd = true;
cfg.final_scalar_ung = true;
cfg.tfr_high_freq = 30:2:120;
cfg.tfr_high_win = 0.25*ones(size(cfg.tfr_high_freq));
% The uploaded workflow uses a frequency-dependent window for 2--30 Hz.
cfg.tfr_low_freq = 2:1:30;
cfg.tfr_low_cycles = linspace(3,6,numel(cfg.tfr_low_freq));
cfg.tfr_low_win = cfg.tfr_low_cycles ./ cfg.tfr_low_freq;
cfg.tfr_step = 0.05;
cfg.tfr_trim_percent = 4;
cfg.tfr_edge_padding = 0.5*max([cfg.tfr_high_win cfg.tfr_low_win]) + 0.02;
cfg.long_time = [cfg.tfr_time(1)-cfg.tfr_edge_padding, ...
                 cfg.tfr_time(2)+cfg.tfr_edge_padding];
cfg.tfr_high_smooth = min(10*ones(size(cfg.tfr_high_freq)), ...
    cfg.expected_fsample/2-1-cfg.tfr_high_freq);
assert(all(cfg.tfr_high_smooth>0),'High-frequency TFR smoothing reaches Nyquist.');

cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(cfg.ft_path);
ft_defaults;

cfg.template_file = fullfile(cfg.ft_path,'template','sourcemodel', ...
    'standard_sourcemodel3d5mm.mat');
T = load(cfg.template_file,'sourcemodel');
template_grid = ft_convert_units(T.sourcemodel,'mm');

cfg.anatomy_dir = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\head_model';
cfg.data_root = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
cfg.lf_root = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\lead_field_pressPutton_3run(250hz)';
cfg.behavior_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\pressPutton102(250hz)\data_vdx.mat';
cfg.head_file = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\out_run_4condition\checkHead_idx_2day.mat';

cfg.condition_dirs = {
    fullfile(cfg.data_root,'social_OT_306(250hz)')
    fullfile(cfg.data_root,'social_PL_306(250hz)')
    fullfile(cfg.data_root,'nonsocial_OT_306(250hz)')
    fullfile(cfg.data_root,'nonsocial_PL_306(250hz)')};

cfg.lf_dirs = {
    fullfile(cfg.lf_root,'social-OT(250hz)')
    fullfile(cfg.lf_root,'social-PL(250hz)')
    fullfile(cfg.lf_root,'nonsocial-OT(250hz)')
    fullfile(cfg.lf_root,'nonsocial-PL(250hz)')};

cfg.stat_file = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
    'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\stat_cluster.mat'];
cfg.oxtr_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\MEG_mask\article\oxtr_mask\oxtr_50.nii';
cfg.cluster_id = 1;
cfg.oxtr_threshold = 0.5;

cfg.output_dir = fullfile( ...
    ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\source_valid\cross_test\TFRdata_ROI_sigmaV'],'gridwise_raw_sourceSignalSVD');
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

%% Group metadata and ROI
B = load(cfg.behavior_file,'data_vdx');
data_vdx = B.data_vdx;
H = load(cfg.head_file,'optimal_idx','isolated_bad_run');
mri_dir = dir(fullfile(cfg.anatomy_dir,'*_anatomy.mat'));
assert(numel(mri_dir)>=cfg.n_subject,'Not enough anatomy files.');

S = load(cfg.stat_file,'stat_res');
stat_res = S.stat_res;
oxtr = ft_read_mri(cfg.oxtr_file);
ci = [];
ci.parameter = 'anatomy';
ci.interpmethod = 'nearest';
oxtr_grid = ft_sourceinterpolate(ci,oxtr,template_grid);

roi_idx = find(stat_res.posclusterslabelmat(:)==cfg.cluster_id & ...
    oxtr_grid.anatomy(:)>cfg.oxtr_threshold & template_grid.inside(:));
assert(~isempty(roi_idx),'Functional cluster and OXTR mask have no overlap.');
fprintf('ROI grids: %d\n',numel(roi_idx));

roi_conf = cell(cfg.n_subject,4);
roi_cong = cell(cfg.n_subject,4);
% roi_conf/roi_cong remain the high-frequency trial-level TFR for backward
% compatibility.  The following fields contain both frequency branches.
roi_conf_low = cell(cfg.n_subject,4);
roi_cong_low = cell(cfg.n_subject,4);
high_mean_conf = cell(cfg.n_subject,4);
high_mean_cong = cell(cfg.n_subject,4);
low_mean_conf = cell(cfg.n_subject,4);
low_mean_cong = cell(cfg.n_subject,4);
% All-condition results: concatenate all trials from the four conditions
% first, then apply the same robust trimming and averaging procedure.
high_all_conf = cell(cfg.n_subject,1);
high_all_cong = cell(cfg.n_subject,1);
low_all_conf = cell(cfg.n_subject,1);
low_all_cong = cell(cfg.n_subject,1);
high_all_mean_conf = cell(cfg.n_subject,1);
high_all_mean_cong = cell(cfg.n_subject,1);
low_all_mean_conf = cell(cfg.n_subject,1);
low_all_mean_cong = cell(cfg.n_subject,1);
qc_all = cell(cfg.n_subject,1);

%% Subject loop
for s = 1:cfg.n_subject
    sub_name = regexprep(mri_dir(s).name,'_anatomy.mat','');
    fprintf('\nSubject %d/%d: %s\n',s,cfg.n_subject,sub_name);
    sub_idx = find(strcmp(data_vdx(:,1),sub_name));
    assert(isscalar(sub_idx),'Subject ID is not unique.');

    [D,L] = load_runs_local(sub_name,sub_idx,data_vdx,cfg);
    D_ot = D([1:3 7:9]);
    D_pl = D([4:6 10:12]);
    LOT = L([1:3 7:9]);
    LPL = L([4:6 10:12]);
    run_ot_original = (1:6)';
    run_pl_original = (1:6)';

    bad_ot = finite_integer_vector(H.isolated_bad_run.ot(s));
    bad_pl = finite_integer_vector(H.isolated_bad_run.pl(s));
    D_ot(bad_ot) = [];
    D_pl(bad_pl) = [];
    LOT(bad_ot) = [];
    LPL(bad_pl) = [];

    run_ot_original(bad_ot) = [];
    run_pl_original(bad_pl) = [];
    keep_ot = ~cellfun(@isempty,D_ot);
    keep_pl = ~cellfun(@isempty,D_pl);
    D_ot = D_ot(keep_ot);
    D_pl = D_pl(keep_pl);
    LOT = LOT(keep_ot);
    LPL = LPL(keep_pl);
    run_ot_original = run_ot_original(keep_ot);
    run_pl_original = run_pl_original(keep_pl);

    OT = append_runs_local(D_ot);
    PL = append_runs_local(D_pl);

    iot = find(run_ot_original==double(H.optimal_idx.ot(s)),1);
    ipl = find(run_pl_original==double(H.optimal_idx.pl(s)),1);
    assert(~isempty(iot)&&~isempty(ipl), ...
        'Optimal leadfield points to a removed or empty run.');
    lf_ot = ft_convert_units(LOT{iot},'mm');
    lf_pl = ft_convert_units(LPL{ipl},'mm');

    cells_long = split_cells_local(OT,PL);
    cells_short = cell(4,2);
    for c = 1:4
        for k = 1:2
            cells_short{c,k} = ft_redefinetrial( ...
                struct('toilim',cfg.filter_time),cells_long{c,k});
        end

    end

    [Cot_conf,Fot_conf] = pooled_csd_local( ...
        cells_short([1 3],:),true,cfg);
    [Cot_cong,~] = pooled_csd_local( ...
        cells_short([1 3],:),false,cfg);
    [Cpl_conf,Fpl_conf] = pooled_csd_local( ...
        cells_short([2 4],:),true,cfg);
    [Cpl_cong,~] = pooled_csd_local( ...
        cells_short([2 4],:),false,cfg);

    % The filter covariance pools both task classes and both contexts within
    % each acquisition day.  It is therefore common to all four conditions.
    Cjoint = blkdiag((Cot_conf+Cot_cong)/2, ...
                     (Cpl_conf+Cpl_cong)/2);
    [Wot,Wpl,qot,qpl,Qf] = build_grid_filter_local( ...
        Cjoint,lf_ot,lf_pl,OT.label,PL.label,roi_idx,cfg);

    [u,Qsvd] = train_source_svd_local(cells_short,Wot,Wpl,qot,qpl,cfg);

    for c = 1:4
        if ismember(c,[1 3]), W=Wot; q=qot;
        else, W=Wpl; q=qpl; end
        roi_conf{s,c} = reconstruct_tfr_local( ...
            cells_long{c,1},W,q,u,cfg,'conflict','high');
        roi_cong{s,c} = reconstruct_tfr_local( ...
            cells_long{c,2},W,q,u,cfg,'congruent','high');
        roi_conf_low{s,c} = reconstruct_tfr_local( ...
            cells_long{c,1},W,q,u,cfg,'conflict','low');
        roi_cong_low{s,c} = reconstruct_tfr_local( ...
            cells_long{c,2},W,q,u,cfg,'congruent','low');

        high_mean_conf{s,c} = trim_mean_tfr_local(roi_conf{s,c},cfg.tfr_trim_percent);
        high_mean_cong{s,c} = trim_mean_tfr_local(roi_cong{s,c},cfg.tfr_trim_percent);
        low_mean_conf{s,c} = trim_mean_tfr_local(roi_conf_low{s,c},cfg.tfr_trim_percent);
        low_mean_cong{s,c} = trim_mean_tfr_local(roi_cong_low{s,c},cfg.tfr_trim_percent);
    end

    % Combine all conflict and all congruent trials across the four
    % conditions before trimming. This is a trial-weighted grand condition
    % result, rather than an unweighted mean of four condition means.
    high_all_conf{s} = concat_tfr_trials_local(roi_conf(s,:));
    high_all_cong{s} = concat_tfr_trials_local(roi_cong(s,:));
    low_all_conf{s} = concat_tfr_trials_local(roi_conf_low(s,:));
    low_all_cong{s} = concat_tfr_trials_local(roi_cong_low(s,:));

    high_all_mean_conf{s} = trim_mean_tfr_local( ...
        high_all_conf{s},cfg.tfr_trim_percent);
    high_all_mean_cong{s} = trim_mean_tfr_local( ...
        high_all_cong{s},cfg.tfr_trim_percent);
    low_all_mean_conf{s} = trim_mean_tfr_local( ...
        low_all_conf{s},cfg.tfr_trim_percent);
    low_all_mean_cong{s} = trim_mean_tfr_local( ...
        low_all_cong{s},cfg.tfr_trim_percent);

    Qsvd.subject = sub_name;
    Qsvd.filter = Qf;
    Qsvd.filter_frequency = Fot_conf;
    Qsvd.filter_frequency_pl = Fpl_conf;
    % Qf.grid_idx contains only grids for which both OT and PL filters were
    % successfully constructed and passed the gain checks.
    Qsvd.grid_idx = Qf.grid_idx;
    Qsvd.apply_day_unit_gain = cfg.apply_day_unit_gain;
    Qsvd.grid_ung_for_source_svd = cfg.grid_ung_for_source_svd;
    Qsvd.final_scalar_ung = cfg.final_scalar_ung;
    qc_all{s} = Qsvd;
    fprintf('PC1 explained variance: %.2f%%\n',100*Qsvd.explained);
end

cfg_saved = cfg;
% Condition aliases follow the uploaded naming convention:
% 1=OT-social, 2=PL-social, 3=OT-nonsocial, 4=PL-nonsocial.
% These aliases intentionally point to the trimmed condition means, matching
% the uploaded script's high_* / low_* save convention.  The unaveraged
% trial-level data remain available in roi_conf/roi_cong and their low fields.
high_os_tfr_conf = high_mean_conf(:,1); high_os_tfr_cong = high_mean_cong(:,1);
high_ps_tfr_conf = high_mean_conf(:,2); high_ps_tfr_cong = high_mean_cong(:,2);
high_on_tfr_conf = high_mean_conf(:,3); high_on_tfr_cong = high_mean_cong(:,3);
high_pn_tfr_conf = high_mean_conf(:,4); high_pn_tfr_cong = high_mean_cong(:,4);
low_os_tfr_conf = low_mean_conf(:,1); low_os_tfr_cong = low_mean_cong(:,1);
low_ps_tfr_conf = low_mean_conf(:,2); low_ps_tfr_cong = low_mean_cong(:,2);
low_on_tfr_conf = low_mean_conf(:,3); low_on_tfr_cong = low_mean_cong(:,3);
low_pn_tfr_conf = low_mean_conf(:,4); low_pn_tfr_cong = low_mean_cong(:,4);
high_all_tfr_conf = high_all_mean_conf;
high_all_tfr_cong = high_all_mean_cong;
low_all_tfr_conf = low_all_mean_conf;
low_all_tfr_cong = low_all_mean_cong;

% Merge only the trial-averaged TFRs. Preserve high_mean_* and low_mean_*
% as the original band-specific results; high_* aliases become broadband.
high_alias_conf={high_os_tfr_conf,high_ps_tfr_conf,high_on_tfr_conf,high_pn_tfr_conf};
high_alias_cong={high_os_tfr_cong,high_ps_tfr_cong,high_on_tfr_cong,high_pn_tfr_cong};
low_alias_conf={low_os_tfr_conf,low_ps_tfr_conf,low_on_tfr_conf,low_pn_tfr_conf};
low_alias_cong={low_os_tfr_cong,low_ps_tfr_cong,low_on_tfr_cong,low_pn_tfr_cong};
for c=1:4
    for j=1:numel(high_alias_conf{c})
        high_alias_conf{c}{j}=merge_low_high_tfr_crop_local( ...
            low_alias_conf{c}{j},high_alias_conf{c}{j},cfg.output_tfr_time);
        high_alias_cong{c}{j}=merge_low_high_tfr_crop_local( ...
            low_alias_cong{c}{j},high_alias_cong{c}{j},cfg.output_tfr_time);
    end
end
os_tfr_conf=high_alias_conf{1}; os_tfr_cong=high_alias_cong{1};
ps_tfr_conf=high_alias_conf{2}; ps_tfr_cong=high_alias_cong{2};
on_tfr_conf=high_alias_conf{3}; on_tfr_cong=high_alias_cong{3};
pn_tfr_conf=high_alias_conf{4}; pn_tfr_cong=high_alias_cong{4};
all_tfr_conf = cell(cfg.n_subject,1);
all_tfr_cong = cell(cfg.n_subject,1);

for j = 1:cfg.n_subject
    all_tfr_conf{j} = merge_low_high_tfr_crop_local( ...
        low_all_tfr_conf{j},high_all_tfr_conf{j},cfg.output_tfr_time);

    all_tfr_cong{j} = merge_low_high_tfr_crop_local( ...
        low_all_tfr_cong{j},high_all_tfr_cong{j},cfg.output_tfr_time);
end

subject_ids = cellfun(@(x)x.subject,qc_all,'UniformOutput',false);
save(fullfile(cfg.output_dir,'gridwise_raw_sourceSignalSVD_TFR.mat'), ...
    'os_tfr_conf','os_tfr_cong', ...
    'ps_tfr_conf','ps_tfr_cong', ...
    'on_tfr_conf','on_tfr_cong', ...
    'pn_tfr_conf','pn_tfr_cong', ...
    'all_tfr_conf','all_tfr_cong', ...
    'subject_ids','cfg_saved','roi_idx','-v7.3');

%% ----------------------------- Functions -------------------------------
function [D,L] = load_runs_local(name,si,b,cfg)
D=cell(1,12); L=cell(1,12); n=1;
for c=1:4
    bh=b{si,c+1};
    M=load_one_main_local(dir(fullfile(cfg.condition_dirs{c},[name '*'])));
    assert(numel(M.trial)==numel(bh.is_valid),'MEG/behavior trial mismatch.');
    LF=load(fullfile(dir(fullfile(cfg.lf_dirs{c},[name '*'])).folder, ...
        dir(fullfile(cfg.lf_dirs{c},[name '*'])).name));
    for r=1:3
        ii=find(double(bh.is_valid(:))==1 & double(bh.run_idx(:))==r);
        if isempty(ii),n=n+1;continue;end
        X=ft_selectdata(struct('trials',ii),M);
        X=ft_redefinetrial(struct('toilim',cfg.long_time),X);
        X.trialinfo=[double(bh.is_conflict(ii)),repmat(c,numel(ii),1)];
        D{n}=X; L{n}=LF.leadfield_individual{r}; n=n+1;
    end
end
end

function M=load_one_main_local(d)
assert(isscalar(d),'Expected one MEG file.');
S=load(fullfile(d.folder,d.name)); f=fieldnames(S);
assert(isscalar(f),'MEG MAT must contain one main variable.');
M=S.(f{1});
end

function out=append_runs_local(x)
x=x(~cellfun(@isempty,x)); assert(~isempty(x),'No valid runs.');
% Raw files contain 306 MEG channels plus three auxiliary channels.  The
% leadfields and DICS filters in this pipeline use the 204 planar/gradiometer
% channels only, so remove non-MEG and magnetometer channels before comparing
% labels or appending runs.
for k=1:numel(x)
    cfg_grad=[];
    cfg_grad.channel='MEGGRAD';
    x{k}=ft_selectdata(cfg_grad,x{k});
end

ref=cellstr(x{1}.label(:));
assert(numel(unique(ref))==numel(ref), ...
    'Reference run contains duplicate channel labels.');
for k=1:numel(x)
    lab=cellstr(x{k}.label(:));
    assert(numel(unique(lab))==numel(lab), ...
        'Run %d contains duplicate channel labels.',k);
    [tf,ord]=ismember(ref,lab);
    assert(all(tf) && numel(lab)==numel(ref), ...
        'Run %d has missing or extra channels relative to the reference run.',k);
    if ~isequal(ord(:),(1:numel(ref))')
        % Reorder the numeric trial matrices explicitly.  This avoids a
        % FieldTrip-version-dependent behavior in ft_selectdata, which may
        % select the requested channels but retain their original order.
        for it=1:numel(x{k}.trial)
            X=x{k}.trial{it};
            assert(size(X,1)==numel(lab), ...
                'Run %d trial %d has an invalid channel dimension.',k,it);
            x{k}.trial{it}=X(ord,:);
        end
        x{k}.label=x{k}.label(ord);

        % Keep sensor metadata in the same order when it is present.
        if isfield(x{k},'grad') && isfield(x{k}.grad,'label')
            glab=cellstr(x{k}.grad.label(:));
            [gtf,gord]=ismember(ref,glab);
            if all(gtf) && numel(glab)==numel(ref)
                x{k}.grad.label=x{k}.grad.label(gord);
                x{k}.grad=reorder_grad_fields_local(x{k}.grad,gord);
            end
        end

        % Header channel metadata, when present, follows the data-channel
        % order. Trial/sample metadata are deliberately left unchanged.
        if isfield(x{k},'hdr') && isstruct(x{k}.hdr) && ...
                isfield(x{k}.hdr,'label')
            hlabel=cellstr(x{k}.hdr.label(:));
            [htf,hord]=ismember(ref,hlabel);
            if all(htf) && numel(hlabel)==numel(ref)
                x{k}.hdr.label=x{k}.hdr.label(hord);
                for hf={'chantype','chanunit'}
                    fn=hf{1};
                    if isfield(x{k}.hdr,fn) && ...
                            numel(x{k}.hdr.(fn))==numel(hord)
                        x{k}.hdr.(fn)=x{k}.hdr.(fn)(hord);
                    end
                end
            end
        end
    end
    assert(isequal(cellstr(x{k}.label(:)),ref), ...
        'Explicit channel reordering failed for run %d.',k);
end
out=ft_appenddata(struct('keepsampleinfo','no'),x{:});
end

function grad=reorder_grad_fields_local(grad,ord)
% Reorder row-wise sensor metadata without changing coil geometry itself.
% chanpos/chanori/chantype/chanunit are channel-level fields.
% coilpos and coilori are coil-level fields and must not be reordered by
% channel index merely because their row count happens to match.
row_fields={'chanpos','chanori','chantype','chanunit'};
for i=1:numel(row_fields)
    f=row_fields{i};
    if isfield(grad,f)
        v=grad.(f);
        if size(v,1)==numel(ord)
            grad.(f)=v(ord,:);
        elseif isvector(v) && numel(v)==numel(ord)
            grad.(f)=v(ord);
        end
    end
end
if isfield(grad,'tra') && size(grad.tra,1)==numel(ord)
    grad.tra=grad.tra(ord,:);
end
end

function cells=split_cells_local(OT,PL)
cells=cell(4,2);
for c=1:4
    if ismember(c,[1 3]),D=OT;else,D=PL;end
    % trialinfo(:,1): 1=conflict, 0=congruent
    % trialinfo(:,2): condition code, 1=OT-social, 2=PL-social,
    %                 3=OT-nonsocial, 4=PL-nonsocial
    idx_condition = D.trialinfo(:,2)==c;

    idx_conflict = find(idx_condition & D.trialinfo(:,1)==1);
    idx_congruent = find(idx_condition & D.trialinfo(:,1)==0);

    assert(~isempty(idx_conflict), ...
        'Condition %d has no valid conflict trials.',c);
    assert(~isempty(idx_congruent), ...
        'Condition %d has no valid congruent trials.',c);

    cells{c,1}=ft_selectdata( ...
        struct('trials',idx_conflict),D);
    cells{c,2}=ft_selectdata( ...
        struct('trials',idx_congruent),D);

    assert(all(cells{c,1}.trialinfo(:,2)==c) && ...
           all(cells{c,2}.trialinfo(:,2)==c), ...
        'Condition split failed for condition %d.',c);
end
end

function [C,F]=pooled_csd_local(cells,is_conflict,cfg)
% cells are [social; nonsocial] for one day. Each row contains
% {conflict, congruent}.  For the selected trial type, concatenate all
% social and nonsocial trials first, then estimate one pooled CSD.  This
% matches the whole-brain DICS implementation with pooled-trial weighting.
if is_conflict
    col=1;
else
    col=2;
end
D_social = cells{1,col};
D_nonsocial = cells{2,col};
assert(~isempty(D_social.trial) && ~isempty(D_nonsocial.trial), ...
    'Both social and nonsocial trial sets are required for pooled CSD.');

% append_runs_local also verifies/reorders the 204 MEGGRAD labels.
D_pool = append_runs_local({D_social,D_nonsocial});

fc=[]; fc.method='mtmfft'; fc.output='fourier'; fc.taper='dpss';
fc.channel='MEGGRAD'; fc.foi=mean(cfg.target_foi);
fc.tapsmofrq=diff(cfg.target_foi)/2; fc.keeptrials='yes';
F=ft_freqanalysis(fc,D_pool);
C=csd_local(F);
end

function C=csd_local(F)
X=F.fourierspctrm(:,:,1);
C=(X'*X)/size(X,1);
C=(C+C')/2;
end

function [Wot,Wpl,qot,qpl,Q]=build_grid_filter_local( ...
    C,L1,L2,lab1,lab2,roi,cfg)
[tf1,i1]=ismember(lab1,L1.label(:));
[tf2,i2]=ismember(lab2,L2.label(:));
assert(all(tf1)&&all(tf2),'CSD/leadfield labels do not match.');
scale=real(trace(C))/size(C,1);
invC=pinv((C+C')/2+cfg.lambda*scale*eye(size(C)));
ng=numel(roi); Wot=cell(ng,1);Wpl=cell(ng,1);
qot=nan(ng,1);qpl=nan(ng,1); Q=struct(); Q.n_valid=0;
Q.gain=cell(ng,2);
for k=1:ng
    g=roi(k); A=double(L1.leadfield{g}(i1,:));
    B=double(L2.leadfield{g}(i2,:));
    if isempty(A)||isempty(B),continue;end
    Lj=[A;B]; D=Lj'*invC*Lj;
    W=pinv(D)*(Lj'*invC);
    W1=W(:,1:numel(lab1)); W2=W(:,numel(lab1)+1:end);
    [W1,Q1]=gain_local(W1,A,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);
    [W2,Q2]=gain_local(W2,B,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);
    if any(~isfinite([W1(:);W2(:)])),continue;end
    Wot{k}=W1;Wpl{k}=W2;qot(k)=real(trace(W1*W1'));
    qpl(k)=real(trace(W2*W2'));Q.gain{k,1}=Q1;Q.gain{k,2}=Q2;
    Q.n_valid=Q.n_valid+1;
end
ok=~cellfun(@isempty,Wot)&~cellfun(@isempty,Wpl)&qot>0&qpl>0;
assert(any(ok),'No valid grid filter.');
Wot=Wot(ok);Wpl=Wpl(ok);qot=qot(ok);qpl=qpl(ok);
Q.grid_idx=roi(ok);Q.n_valid=sum(ok);Q.rank=3;
end

function [u,Q]=train_source_svd_local(cells,Wot,Wpl,qot,qpl,cfg)
ng=numel(Wot); r=size(Wot{1},1); nf=ng*r;
R=zeros(nf,nf); trc=zeros(4,2);
for c=1:4
    if ismember(c,[1 3]),W=Wot;q=qot;else,W=Wpl;q=qpl;end
    for k=1:2
        D=cells{c,k}; RR=zeros(nf,nf); nobs=0;
        for t=1:numel(D.trial)
            X=double(D.trial{t}); if any(~isfinite(X(:))),continue;end
            Z=grid_project_local(X,W,q,cfg.grid_ung_for_source_svd);
            RR=RR+Z*Z'; nobs=nobs+size(Z,2);
        end
        assert(nobs>0,'Empty source-SVD cell.');
        RR=RR/nobs; trc(c,k)=real(trace(RR)); R=R+RR/8;
    end
end
R=(R+R')/2; [U,d]=eig(R,'vector'); [d,ord]=sort(real(d),'descend');
assert(d(1)>0,'Invalid source second moment.');
u=real(U(:,ord(1))); [~,a]=max(abs(u)); if u(a)<0,u=-u;end
Q.explained=d(1)/sum(max(d,0));Q.eigenvalues=d;
Q.n_grid=ng;Q.n_orientation=r;Q.cell_trace=trc;
Q.method='gridwise raw source projection followed by source-signal SVD';
end

function Z=grid_project_local(X,W,q,useung)
ng=numel(W); r=size(W{1},1); Z=zeros(ng*r,size(X,2));
for k=1:ng
    w=W{k}; if useung,w=w/sqrt(q(k));end
    Z((k-1)*r+(1:r),:)=w*X;
end
end

function T=reconstruct_tfr_local(D,W,q,u,cfg,tag,band)
V=D;V.label={'ROI_sourceSignalSVD'};V.trial=cell(numel(D.trial),1);
for t=1:numel(D.trial)
    Z=grid_project_local(double(D.trial{t}),W,q,cfg.grid_ung_for_source_svd);
    y=real(u'*Z);
    V.trial{t}=y;
end
if cfg.final_scalar_ung
    w_eff = effective_filter_local(W,q,u,cfg.grid_ung_for_source_svd);
    q_eff = real(w_eff*w_eff');
    assert(isfinite(q_eff)&&q_eff>0,'Invalid final scalar white-noise gain.');
    for t=1:numel(V.trial),V.trial{t}=V.trial{t}/sqrt(q_eff);end
end
fc=[];fc.method='mtmconvol';fc.output='pow';fc.channel='all';
switch lower(band)
    case 'high'
        fc.taper='dpss';
        fc.foi=cfg.tfr_high_freq;
        fc.t_ftimwin=cfg.tfr_high_win;
        fc.tapsmofrq=cfg.tfr_high_smooth;
        fc.toi=cfg.tfr_time_high(1):cfg.tfr_step:cfg.tfr_time_high(2);
    case 'low'
        fc.taper='hanning';
        fc.foi=cfg.tfr_low_freq;
        fc.t_ftimwin=cfg.tfr_low_win;
        fc.toi=cfg.tfr_time_low(1):cfg.tfr_step:cfg.tfr_time_low(2);
    otherwise
        error('Unknown TFR band: %s',band);
end
fc.keeptrials='yes';fc.pad='nextpow2';
T=ft_freqanalysis(fc,V);
% FieldTrip cfg.previous can retain a large processing history. Analysis
% settings are saved once in cfg_saved, rather than in every trial TFR.
if isfield(T,'cfg'), T=rmfield(T,'cfg'); end
if ~isfield(T,'trialinfo'), T.trialinfo=D.trialinfo; end
T=crop_tfr_time_local(T,cfg.output_tfr_time);
T.source_svd_method=tag;
T.tfr_band=band;
end

function T=crop_tfr_time_local(T,timewin)
% Crop the already computed TFR without recomputing the convolution.
assert(numel(timewin)==2 && timewin(2)>timewin(1), ...
    'Invalid output TFR time window.');
idx=T.time>=timewin(1)-1e-10 & T.time<=timewin(2)+1e-10;
assert(any(idx),'Output TFR time window has no sampled time points.');
T.time=T.time(idx);
if isfield(T,'powspctrm')
    if ndims(T.powspctrm)==4
        T.powspctrm=T.powspctrm(:,:,:,idx);
    else
        T.powspctrm=T.powspctrm(:,:,idx);
    end
end
end

function Tm=trim_mean_tfr_local(T,pct)
% Robust point-wise MAD/RMS trial rejection, followed by a trial mean.
assert(isfield(T,'powspctrm'),'TFR has no powspctrm.');
assert(strcmp(T.dimord,'rpt_chan_freq_time'), ...
    'Expected trial-level TFR with dimord rpt_chan_freq_time.');
P=double(T.powspctrm);
P=reshape(P,size(P,1),numel(T.label),numel(T.freq),numel(T.time));
n=size(P,1);
[Ptrim,removed_idx]=trim_by_percent_strict(P,1,pct);
kept_idx=setdiff((1:n)',removed_idx(:));
Tm=T;
if isfield(Tm,'cfg'), Tm=rmfield(Tm,'cfg'); end
Tm.dimord='chan_freq_time';
Pm=mean(Ptrim,1,'omitnan');
Tm.powspctrm=reshape(Pm,[numel(T.label) numel(T.freq) numel(T.time)]);
Tm.trim_percent=pct;
Tm.out_idx=removed_idx;
Tm.removed_trial_idx=removed_idx;
Tm.kept_trial_idx=kept_idx;
Tm.n_trials_before=n;
Tm.n_trials_used=numel(kept_idx);
Tm.actual_trim_percent=100*numel(removed_idx)/n;
Tm.trim_method='pointwise median/MAD robust Z-score aggregated by RMS';
trial_fields={'trialinfo','cumtapcnt','cumsumcnt'};
for j=1:numel(trial_fields)
    if isfield(Tm,trial_fields{j}), Tm=rmfield(Tm,trial_fields{j}); end
end
end

function [trimmed_data,removed_idx]=trim_by_percent_strict(data,dim,percent)
% Same scoring/count rule as the uploaded function for percent>0 and n>1.
% Kept local so the reconstruction remains a single .m file.
validateattributes(percent,{'numeric'},{'scalar','real','finite','>=',0,'<=',100});
validateattributes(dim,{'numeric'},{'scalar','integer','>=',1,'<=',ndims(data)});
sz=size(data); n=sz(dim);
assert(n>=1,'Cannot trim an empty trial dimension.');
if percent==0 || n==1
    trimmed_data=data; removed_idx=[]; return;
end
n_remove=max(1,min(round(n*percent/100),n-1));
perm_order=1:ndims(data);
perm_order(perm_order==dim)=[];
perm_order=[dim perm_order];
data_perm=permute(data,perm_order);
matrix_data=reshape(data_perm,n,[]);
assert(~any(isinf(matrix_data(:))),'TFR contains Inf; inspect the input data.');
voxel_med=median(matrix_data,1,'omitnan');
voxel_mad=median(abs(matrix_data-voxel_med),1,'omitnan');
zero_mad_mask=voxel_mad==0;
if any(zero_mad_mask)
    voxel_std=std(matrix_data(:,zero_mad_mask),1,'omitnan');
    voxel_mad(zero_mad_mask)=voxel_std+eps;
end
Z_matrix=abs(matrix_data-voxel_med)./(1.4826*voxel_mad);
slice_scores=sqrt(mean(Z_matrix.^2,2,'omitnan'));
% Explicitly rank completely missing trials as most abnormal.
slice_scores(~isfinite(slice_scores))=Inf;
[~,sorted_idx]=sort(slice_scores,'descend');
removed_idx=sort(sorted_idx(1:n_remove)');
kept_idx=setdiff(1:n,removed_idx);
new_sz=sz; new_sz(dim)=numel(kept_idx);
trimmed_perm=reshape(matrix_data(kept_idx,:),new_sz(perm_order));
[~,inv_perm]=sort(perm_order);
trimmed_data=permute(trimmed_perm,inv_perm);
fprintf('Robust trial trimming: removed %d/%d trials (%.2f%%).\n', ...
    n_remove,n,100*n_remove/n);
end

function T=concat_tfr_trials_local(Tcell)
% Concatenate keeptrials='yes' TFR structures along the trial dimension.
% The input is one subject's four condition-specific TFRs.
Tcell=Tcell(~cellfun(@isempty,Tcell));
assert(~isempty(Tcell),'No condition-level TFR is available.');
T=Tcell{1};
if isfield(T,'cfg'), T=rmfield(T,'cfg'); end
assert(strcmp(T.dimord,'rpt_chan_freq_time'), ...
    'Expected rpt_chan_freq_time TFR for trial concatenation.');
P=T.powspctrm;
for k=2:numel(Tcell)
    Tk=Tcell{k};
    assert(strcmp(Tk.dimord,'rpt_chan_freq_time'), ...
        'All TFR inputs must have rpt_chan_freq_time dimord.');
    assert(isequal(T.label,Tk.label) && isequal(T.freq,Tk.freq) && ...
        isequal(T.time,Tk.time), ...
        'TFR grids differ across conditions.');
    P=cat(1,P,Tk.powspctrm);
end
T.powspctrm=P;
T.dimord='rpt_chan_freq_time';
T.trialinfo=[];
T.n_trials_total=size(P,1);
end

function [Wc,Q]=gain_local(W,L,apply,rtol)
G=W*L;[U,S,V]=svd(G,'econ');s=real(diag(S));
keep=s>rtol*max(s);
if isempty(s)||~any(keep),Wc=nan(size(W));Q=[];return;end
P=V(:,keep)*diag(1./s(keep))*U(:,keep)';
if apply,Wc=P*W;else,Wc=W;end
Q.applied=apply;Q.supported_rank=sum(keep);
Q.error_before=norm(G-P*G,'fro')/max(norm(P*G,'fro'),realmin);
Q.error_after=norm(Wc*L-P*G,'fro')/max(norm(P*G,'fro'),realmin);
end

function v=finite_integer_vector(x)
if iscell(x), x=cell2mat(x(:)); end
v=double(x(:));v=v(isfinite(v)&v==round(v)&v>=1);
end

function j=adjusted_index(x,bad,n)
j=double(x);j=j-sum(bad<j);assert(j>=1&&j<=n,'Bad adjusted leadfield index.');
end

function tfr_merge = merge_low_high_tfr_crop_local(tfr_low, tfr_high, timewin)
% =========================================================================
% Combine low and high frequencies TFR, and crop to the specified time window
%
% New rule: 
%   freq < 31 Hz  -> Use low-frequency TFR
%   freq > 31 Hz  -> Use high-frequency TFR
%   freq == 31 Hz -> Discard
%
% Input requirements: 
%   tfr_low.powspctrm  : chan × freq × time
%   tfr_high.powspctrm : chan × freq × time
%
% Output: 
%   freq Sorted in ascending order
%   powspctrm frequency dimension and freq order must match exactly
% =========================================================================

assert(isfield(tfr_low, 'powspctrm') && isfield(tfr_high, 'powspctrm'), ...
    '输入必须包含 powspctrm。');

assert(isfield(tfr_low, 'freq') && isfield(tfr_high, 'freq'), ...
    '输入必须包含 freq。');

assert(isfield(tfr_low, 'time') && isfield(tfr_high, 'time'), ...
    '输入必须包含 time。');

assert(numel(timewin) == 2 && timewin(2) > timewin(1), ...
    'timewin 必须是 [start end]，例如 [-1.2 0.8]。');

if isfield(tfr_low, 'dimord')
    assert(strcmp(tfr_low.dimord, 'chan_freq_time'), ...
        '当前函数只支持 tfr_low.dimord = chan_freq_time。');
end

if isfield(tfr_high, 'dimord')
    assert(strcmp(tfr_high.dimord, 'chan_freq_time'), ...
        '当前函数只支持 tfr_high.dimord = chan_freq_time。');
end

assert(isequal(tfr_low.label(:),tfr_high.label(:)), ...
    'Low/high TFR channel labels or order differ.');

tol = 1e-6;
cut_freq = 31;

low_freq_raw  = tfr_low.freq(:)';
high_freq_raw = tfr_high.freq(:)';

% ---------------------------------------------------------
% 0. Select frequencies by the rule: 
%    Retain only low-frequency <30, Retain only high-frequency >30
% ---------------------------------------------------------
low_keep_raw  = low_freq_raw  < cut_freq - tol;
high_keep_raw = high_freq_raw > cut_freq + tol;

if ~any(low_keep_raw)
    error('低频 TFR 中没有 freq < %.1f Hz 的频点。当前 low freq 为：%s', ...
        cut_freq, mat2str(low_freq_raw));
end

if ~any(high_keep_raw)
    error('高频 TFR 中没有 freq > %.1f Hz 的频点。当前 high freq 为：%s', ...
        cut_freq, mat2str(high_freq_raw));
end

low_freq_keep_raw  = low_freq_raw(low_keep_raw);
high_freq_keep_raw = high_freq_raw(high_keep_raw);

% ---------------------------------------------------------
% 1. Check for duplicates within retained frequencies
% ---------------------------------------------------------
low_freq_check = sort(low_freq_keep_raw);
high_freq_check = sort(high_freq_keep_raw);

if any(diff(low_freq_check) <= tol)
    error('保留后的 low freq 内部存在重复频点：%s', ...
        mat2str(low_freq_keep_raw));
end

if any(diff(high_freq_check) <= tol)
    error('保留后的 high freq 内部存在重复频点：%s', ...
        mat2str(high_freq_keep_raw));
end

% ---------------------------------------------------------
% 2. Retain low-frequency <30 Hz, and sort ascending; Reorder concurrently power
% ---------------------------------------------------------
[low_freq, low_sort_idx] = sort(low_freq_keep_raw, 'ascend');
low_keep_idx = find(low_keep_raw);
low_freq_idx_sorted = low_keep_idx(low_sort_idx);

low_pow_sorted = double(tfr_low.powspctrm(:, low_freq_idx_sorted, :));

% ---------------------------------------------------------
% 3. Retain high-frequency >30 Hz, and sort ascending; Reorder concurrently power
% ---------------------------------------------------------
[high_freq, high_sort_idx] = sort(high_freq_keep_raw, 'ascend');
high_keep_idx = find(high_keep_raw);
high_freq_idx_sorted = high_keep_idx(high_sort_idx);

high_time_idx = find(tfr_high.time >= timewin(1)-tol & ...
                     tfr_high.time <= timewin(2)+tol);

if isempty(high_time_idx)
    error('高频 TFR 中没有落在目标时间窗 [%.3f %.3f] 内的时间点。', ...
        timewin(1), timewin(2));
end

high_pow_crop = double(tfr_high.powspctrm(:, high_freq_idx_sorted, high_time_idx));

% ---------------------------------------------------------
% 4. Use the high-frequency time axis within timewin the time interval
% ---------------------------------------------------------
merge_time = tfr_high.time(high_time_idx);
merge_time = merge_time(:)';

% ---------------------------------------------------------
% 5. Low-frequency power Interpolate onto the high-frequency time axis
% ---------------------------------------------------------
low_time = tfr_low.time(:)';

if min(merge_time) < min(low_time)-tol || max(merge_time) > max(low_time)+tol
    error(['低频 TFR 时间范围不能覆盖最终合并时间窗。' ...
           '低频范围为 [%.3f %.3f]，目标范围为 [%.3f %.3f]。'], ...
        min(low_time), max(low_time), min(merge_time), max(merge_time));
end

n_chan = size(low_pow_sorted, 1);
n_low_freq = size(low_pow_sorted, 2);
n_low_time = size(low_pow_sorted, 3);

low_pow_2d = reshape(low_pow_sorted, n_chan * n_low_freq, n_low_time)';
low_pow_interp_2d = interp1(low_time, low_pow_2d, merge_time, 'linear');

if any(isnan(low_pow_interp_2d(:)))
    warning('低频插值后出现 NaN，请检查 low_time 是否完整覆盖 merge_time。');
end

low_pow_crop = reshape(low_pow_interp_2d', ...
    n_chan, n_low_freq, numel(merge_time));

% ---------------------------------------------------------
% 6. Concatenate the frequency dimension
% ---------------------------------------------------------
merge_freq_pre = [low_freq, high_freq];
merge_pow_pre = cat(2, low_pow_crop, high_pow_crop);

% ---------------------------------------------------------
% 7. Sort globally again after concatenation, and reorder concurrently power
% ---------------------------------------------------------
[merge_freq, merge_sort_idx] = sort(merge_freq_pre, 'ascend');
merge_pow = merge_pow_pre(:, merge_sort_idx, :);

% ---------------------------------------------------------
% 8. Final checks: Frequency order, power Dimension, 30 Hz Whether it was removed
% ---------------------------------------------------------
if any(diff(merge_freq) <= tol)
    error('合并后频率轴仍存在重复或倒序。当前 merge_freq 为：%s', ...
        mat2str(merge_freq));
end

if any(abs(merge_freq - cut_freq) <= tol)
    error('合并后仍包含 %.1f Hz，请检查频率筛选逻辑。当前 merge_freq 为：%s', ...
        cut_freq, mat2str(merge_freq));
end

if any(merge_freq < cut_freq - tol) && any(merge_freq > cut_freq + tol)
    % Normal: Both frequency bands contain
else
    warning('合并后似乎只包含低频或只包含高频，请检查输入数据。');
end

if size(merge_pow, 2) ~= numel(merge_freq)
    error('powspctrm 的频率维度数量与 freq 数量不一致。');
end

if any(isnan(merge_pow(:)))
    warning('合并后的 powspctrm 中存在 NaN。');
end

% ---------------------------------------------------------
% 9. Organize output structure
% ---------------------------------------------------------
tfr_merge = [];
tfr_merge.label = tfr_low.label;
tfr_merge.dimord = 'chan_freq_time';
tfr_merge.freq = merge_freq;
tfr_merge.time = merge_time;
tfr_merge.powspctrm = merge_pow;

% Auxiliary information
tfr_merge.low_freq = low_freq;
tfr_merge.high_freq = high_freq;
tfr_merge.cut_freq = cut_freq;
tfr_merge.freq_rule = ...
    'freq < 31 Hz uses low-frequency TFR; freq > 31 Hz uses high-frequency TFR; freq == 31 Hz is removed.';

tfr_merge.original_low_freq = low_freq_raw;
tfr_merge.original_high_freq = high_freq_raw;
tfr_merge.original_low_time = tfr_low.time;
tfr_merge.original_high_time = tfr_high.time;
tfr_merge.merge_timewin = timewin;

tfr_merge.merge_note = ...
    'Merged low+high TFR. Low keeps freq < 31 Hz, high keeps freq > 31 Hz. Frequencies are sorted ascending and powspctrm is reordered with the same index.';

if isfield(tfr_low, 'out_idx')
    tfr_merge.low_out_idx = tfr_low.out_idx;
end

if isfield(tfr_high, 'out_idx')
    tfr_merge.high_out_idx = tfr_high.out_idx;
end

fprintf('✅ 合并完成：低频 %.1f–%.1f Hz，高频 %.1f–%.1f Hz；最终 freq %.1f–%.1f Hz，共 %d 个频点；time %.3f–%.3f s，共 %d 个时间点。\n', ...
    min(low_freq), max(low_freq), ...
    min(high_freq), max(high_freq), ...
    min(tfr_merge.freq), max(tfr_merge.freq), numel(tfr_merge.freq), ...
    min(tfr_merge.time), max(tfr_merge.time), numel(tfr_merge.time));

end

function w=effective_filter_local(W,q,u,useung)
ng=numel(W); r=size(W{1},1); nc=size(W{1},2);
Wflat=zeros(ng*r,nc);
for k=1:ng
    wk=W{k};
    if useung, wk=wk/sqrt(q(k)); end
    Wflat((k-1)*r+(1:r),:)=wk;
end
w=u(:)'*Wflat;
end
















