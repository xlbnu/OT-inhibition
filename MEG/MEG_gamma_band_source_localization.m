%% Whole-brain DICS matched to the filter-SVD reconstruction pipeline
% Pipeline
%   1) Compute OT/PL conflict/congruent CSDs in the reconstruction filter
%      window and with exactly the same frequency definition.
%   2) Average conflict/congruent CSDs equally within day and build one
%      408 x 408 block-diagonal common CSD.
%   3) Estimate the same 3-row joint DICS filter used before ROI filter-SVD
%      at every whole-brain source grid.
%   4) Split Wjoint into OT and PL channel blocks.
%   5) Optionally calibrate each block separately to unit gain.
%   6) Project conflict/congruent CSDs separately for OT and PL.
%   7) Average the two selected day-filter estimates with equal day weights.
%   8) Optionally use unit-white-noise-gain normalized power as avg.pow.
%
% IMPORTANT
% - The same selected OT filter is used for every OT condition, and the
%   same selected PL filter is used for every PL condition.
% - Do not estimate a different filter or gain matrix for each condition.
% - ROI filter-SVD is intentionally NOT performed here: before a cluster is
%   known, each grid is one candidate source. ROI filter-SVD is a subsequent
%   reconstruction step after significant grids have been selected.

% clear;


%% ----------------------------- Parameters -----------------------------
num_subjects = 19;
cfg_target = struct();
cfg_target.target_foi  = [30 100];
cfg_target.target_time = [-0.5 -0.3];
assert(isnumeric(cfg_target.target_foi) && ...
       numel(cfg_target.target_foi)==2 && ...
       all(isfinite(cfg_target.target_foi)) && ...
       cfg_target.target_foi(1)>0 && ...
       cfg_target.target_foi(2)>cfg_target.target_foi(1), ...
    'target_foi must be a finite increasing positive [low high] pair.');
assert(isnumeric(cfg_target.target_time) && ...
       numel(cfg_target.target_time)==2 && ...
       all(isfinite(cfg_target.target_time)) && ...
       cfg_target.target_time(2)>cfg_target.target_time(1), ...
    'target_time must be a finite increasing [start end] pair.');

cfg_dics = struct();
cfg_dics.lambda = 0.10;
cfg_dics.min_rcond = 1e-10;
cfg_dics.gain_pinv_rtol = 1e-10;

% Independently control the two gain operations.
% true : calibrate the split OT/PL joint-filter blocks separately in the
%        observable (normally rank-2) part of the 3-D source subspace.
% false: use the OT/PL blocks from the joint filter without day calibration.
cfg_dics.apply_day_unit_gain = false;

% false: avg.pow contains raw power from the selected day filters.
% true : avg.pow contains selected-filter unit-white-noise-gain power.
% Both versions are retained in generic pow_selectedfilter_* fields.
cfg_dics.use_unit_white_noise_gain = true;

validateattributes(cfg_dics.apply_day_unit_gain, ...
    {'logical','numeric'},{'scalar'});
validateattributes(cfg_dics.use_unit_white_noise_gain, ...
    {'logical','numeric'},{'scalar'});
cfg_dics.apply_day_unit_gain = logical(cfg_dics.apply_day_unit_gain);
cfg_dics.use_unit_white_noise_gain = ...
    logical(cfg_dics.use_unit_white_noise_gain);

% 'pooled_trials': preserves the original code's context weighting.
% 'equal_context': computes social/nonsocial CSDs separately and averages
%                  them equally within each day and trial type.
cfg_dics.context_weighting = 'pooled_trials';

assert(ismember(cfg_dics.context_weighting,{'pooled_trials','equal_context'}), ...
    'Unknown context weighting mode.');

% ------------------------------- Paths --------------------------------
ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(ft_path);
ft_defaults;

template_grid_file = fullfile(ft_path,'template','sourcemodel', ...
    'standard_sourcemodel3d5mm.mat');
load(template_grid_file,'sourcemodel');
sourcemodel = ft_convert_units(sourcemodel,'mm');

base_anatomy_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\head_model';
mri_dir = dir(fullfile(base_anatomy_path,'*_anatomy.mat'));
assert(numel(mri_dir)>=num_subjects,'Fewer anatomy files than subjects.');

base_leadfield_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\lead_field_pressPutton_3run(250hz)';
lf_path = {
    fullfile(base_leadfield_path,'social-OT(250hz)');
    fullfile(base_leadfield_path,'social-PL(250hz)');
    fullfile(base_leadfield_path,'nonsocial-OT(250hz)');
    fullfile(base_leadfield_path,'nonsocial-PL(250hz)')};

base_data_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
cond_dirs = {
    fullfile(base_data_path,'social_OT_306(250hz)'), ...
    fullfile(base_data_path,'social_PL_306(250hz)'), ...
    fullfile(base_data_path,'nonsocial_OT_306(250hz)'), ...
    fullfile(base_data_path,'nonsocial_PL_306(250hz)')};

load(['E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\' ...
    'MEG_data\pressPutton102(250hz)\data_vdx.mat'],'data_vdx');

head_mov = load(['F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\' ...
    'out_run_4condition\checkHead_idx_2day.mat'], ...
    'optimal_idx','isolated_bad_run');
assert(isfield(head_mov,'optimal_idx') && ...
       isfield(head_mov,'isolated_bad_run') && ...
       isfield(head_mov.optimal_idx,'ot') && ...
       isfield(head_mov.optimal_idx,'pl') && ...
       isfield(head_mov.isolated_bad_run,'ot') && ...
       isfield(head_mov.isolated_bad_run,'pl'), ...
    'Head-movement file lacks required OT/PL fields.');
assert(numel(head_mov.optimal_idx.ot)>=num_subjects && ...
       numel(head_mov.optimal_idx.pl)>=num_subjects && ...
       numel(head_mov.isolated_bad_run.ot)>=num_subjects && ...
       numel(head_mov.isolated_bad_run.pl)>=num_subjects, ...
    'Head-movement arrays contain fewer entries than num_subjects.');

if cfg_dics.apply_day_unit_gain
    gain_tag = 'dayUnitGain';
else
    gain_tag = 'noDayUnitGain';
end
if cfg_dics.use_unit_white_noise_gain
    noise_tag = 'unitWhiteNoiseGain';
else
    noise_tag = 'raw';
end
norm_tag = sprintf('%s_%s',gain_tag,noise_tag);
save_source_path = fullfile( ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain', ...
    sprintf('Source_408DICS_rank3_%s_%g-%gHz_(%g_%g)', ...
    norm_tag,cfg_target.target_foi(1),cfg_target.target_foi(2), ...
    cfg_target.target_time(1),cfg_target.target_time(2)));
if ~exist(save_source_path,'dir'), mkdir(save_source_path); end

% ------------------------- Frequency parameters -----------------------
cfg_freq = [];
cfg_freq.method      = 'mtmfft';
cfg_freq.output      = 'fourier';
cfg_freq.taper       = 'dpss';
cfg_freq.channel     = 'MEGGRAD';
cfg_freq.foi         = mean(cfg_target.target_foi);
cfg_freq.tapsmofrq   = diff(cfg_target.target_foi)/2;
% cfg_freq.pad         = 2;% beta_band Dedicated
cfg_freq.keeptrials  = 'yes';
assert(cfg_freq.foi-cfg_freq.tapsmofrq>0, ...
    'The requested DPSS smoothing band extends to zero/negative frequency.');

fprintf(['DICS frequency definition: centre %.1f Hz, smoothing +/-%.1f Hz, ' ...
    'time window [%.3f %.3f] s.\n'], ...
    cfg_freq.foi,cfg_freq.tapsmofrq, ...
    cfg_target.target_time(1),cfg_target.target_time(2));

%% ----------------------------- Outputs --------------------------------
all_source_conf = cell(num_subjects,1);
all_source_cong = cell(num_subjects,1);
all_source_quality = cell(num_subjects,1);
subject_ids = cell(num_subjects,1);

ot_freq_conf = cell(num_subjects,1);
pl_freq_conf = cell(num_subjects,1);
ot_freq_cong = cell(num_subjects,1);
pl_freq_cong = cell(num_subjects,1);

% --------------------------- Subject loop ------------------------------
for sub = 1:num_subjects
    sub_name = regexprep(mri_dir(sub).name,'_anatomy.mat','');
    subject_ids{sub} = sub_name;
    fprintf('\n=========================================\n');
    fprintf('Subject %d/%d: %s\n',sub,num_subjects,sub_name);

    sub_idx = find(strcmp(data_vdx(:,1),sub_name));
    assert(isscalar(sub_idx),'Cannot uniquely match subject %s in data_vdx.',sub_name);

    all_time_data = cell(1,12);
    all_leadfields = cell(1,12);
    run_counter = 1;

    cfg_win = [];
    cfg_win.toilim = cfg_target.target_time;

    for i_cond = 1:numel(cond_dirs)
        bhv_cond = data_vdx{sub_idx,i_cond+1};
        condition_vdx = [bhv_cond.is_valid,bhv_cond.is_conflict,bhv_cond.run_idx];

        meg_dir = dir(fullfile(cond_dirs{i_cond},[sub_name '*']));
        assert(isscalar(meg_dir), ...
            'Expected one MEG file for %s condition %d, found %d.', ...
            sub_name,i_cond,numel(meg_dir));
        tmpdata = load(fullfile(meg_dir.folder,meg_dir.name));
        var_names = fieldnames(tmpdata);
        assert(isscalar(var_names),'MEG MAT file must contain one main variable.');
        data_process0 = tmpdata.(var_names{1});
        clear tmpdata;

        assert(numel(data_process0.trial)==size(condition_vdx,1), ...
            ['MEG/behaviour trial count mismatch before run selection for %s, ' ...
             'condition %d: MEG=%d, behaviour=%d.'], ...
            sub_name,i_cond,numel(data_process0.trial),size(condition_vdx,1));

        lf_dir = dir(fullfile(lf_path{i_cond},[sub_name '*']));
        assert(isscalar(lf_dir), ...
            'Expected one leadfield file for %s condition %d.',sub_name,i_cond);
        tmp_lf = load(fullfile(lf_dir.folder,lf_dir.name), ...
            'leadfield_individual','sensor_grad0');
        leadfield_3run = tmp_lf.leadfield_individual;
        grad_3run = tmp_lf.sensor_grad0;

        time_data_win = ft_redefinetrial(cfg_win,data_process0);

        % Context coding: 1=social, 2=nonsocial.
        context_id = 1 + (i_cond>=3);

        for i_run = 1:3
            trial_idx = find(condition_vdx(:,3)==i_run & condition_vdx(:,1)==1);

            % Leadfield availability is independent of whether this run has
            % any surviving trials in the present analysis window.
            assert(i_run<=numel(leadfield_3run) && ...
                   ~isempty(leadfield_3run{i_run}), ...
                'Missing leadfield for %s condition %d run %d.', ...
                sub_name,i_cond,i_run);
            all_leadfields{run_counter} = leadfield_3run{i_run};

            if isempty(trial_idx)
                run_counter = run_counter+1;
                continue
            end

            cfg_run = [];
            cfg_run.trials = trial_idx;
            time_data_run = ft_selectdata(cfg_run,time_data_win);
            time_data_run.grad = ft_convert_units(grad_3run{i_run},'mm');
            assert(numel(time_data_run.trial)==numel(trial_idx), ...
                'Trial selection failed for %s condition %d run %d.', ...
                sub_name,i_cond,i_run);

            % Columns: conflict, run, context, original four-condition ID.
            time_data_run.trialinfo = [ ...
                condition_vdx(trial_idx,2), ...
                condition_vdx(trial_idx,3), ...
                repmat(context_id,numel(trial_idx),1), ...
                repmat(i_cond,numel(trial_idx),1)];

            all_time_data{run_counter} = time_data_run;
            run_counter = run_counter+1;
        end
    end

    % Original ordering: social OT(1:3), social PL(4:6),
    %                    nonsocial OT(7:9), nonsocial PL(10:12).
    ot_time_data = all_time_data([1:3 7:9]);
    pl_time_data = all_time_data([4:6 10:12]);
    ot_leadfields = all_leadfields([1:3 7:9]);
    pl_leadfields = all_leadfields([4:6 10:12]);

    bad_ot = unique(cell2mat(head_mov.isolated_bad_run.ot(sub)));
    bad_pl = unique(cell2mat(head_mov.isolated_bad_run.pl(sub)));
    assert(all(isfinite(bad_ot)) && all(bad_ot==round(bad_ot)) && ...
           all(bad_ot>=1 & bad_ot<=6), ...
        'Invalid OT bad-run index for %s.',sub_name);
    assert(all(isfinite(bad_pl)) && all(bad_pl==round(bad_pl)) && ...
           all(bad_pl>=1 & bad_pl<=6), ...
        'Invalid PL bad-run index for %s.',sub_name);
    if ~isempty(bad_ot), ot_time_data(bad_ot)=[]; end
    if ~isempty(bad_pl), pl_time_data(bad_pl)=[]; end
    ot_time_data = ot_time_data(~cellfun(@isempty,ot_time_data));
    pl_time_data = pl_time_data(~cellfun(@isempty,pl_time_data));
    assert(~isempty(ot_time_data) && ~isempty(pl_time_data), ...
        'No valid OT or PL runs remain for %s.',sub_name);

    idx_lf_ot = double(head_mov.optimal_idx.ot(sub));
    idx_lf_pl = double(head_mov.optimal_idx.pl(sub));
    assert(isscalar(idx_lf_ot) && isfinite(idx_lf_ot) && ...
           idx_lf_ot==round(idx_lf_ot) && ...
           isscalar(idx_lf_pl) && isfinite(idx_lf_pl) && ...
           idx_lf_pl==round(idx_lf_pl) && ...
           idx_lf_ot>=1 && idx_lf_ot<=numel(ot_leadfields) && ...
           idx_lf_pl>=1 && idx_lf_pl<=numel(pl_leadfields), ...
        'Optimal leadfield index is outside 1..6.');
    assert(~isempty(ot_leadfields{idx_lf_ot}) && ...
           ~isempty(pl_leadfields{idx_lf_pl}), ...
        'Selected optimal leadfield is empty.');
    assert(~ismember(idx_lf_ot,bad_ot), ...
        'Selected OT leadfield belongs to an excluded bad run for %s.',sub_name);
    assert(~ismember(idx_lf_pl,bad_pl), ...
        'Selected PL leadfield belongs to an excluded bad run for %s.',sub_name);

    ot_lf = ft_convert_units(ot_leadfields{idx_lf_ot},'mm');
    pl_lf = ft_convert_units(pl_leadfields{idx_lf_pl},'mm');
    assert_source_grid_geometry_local(ot_lf,sourcemodel, ...
        sprintf('%s OT',sub_name));
    assert_source_grid_geometry_local(pl_lf,sourcemodel, ...
        sprintf('%s PL',sub_name));

    % ft_appenddata should not be allowed to silently intersect or reorder
    % channels. All retained runs within a day must already have identical
    % channel labels and ordering.
    cfg_app = [];
    cfg_app.keepsampleinfo = 'no';
    cfg_app.appenddim = 'rpt';

    ot_all = ft_appenddata(cfg_app,ot_time_data{:});
    pl_all = ft_appenddata(cfg_app,pl_time_data{:});

    % Compute pooled Fourier structures for channel-level output and the
    % CSDs used by source analysis. With equal_context, the CSD is replaced
    % by the equal average of the two context-specific CSDs.
    [freq_ot_conf,C_ot_conf] = condition_csd_local( ...
        ot_all,1,cfg_freq,cfg_dics.context_weighting);
    [freq_ot_cong,C_ot_cong] = condition_csd_local( ...
        ot_all,0,cfg_freq,cfg_dics.context_weighting);
    [freq_pl_conf,C_pl_conf] = condition_csd_local( ...
        pl_all,1,cfg_freq,cfg_dics.context_weighting);
    [freq_pl_cong,C_pl_cong] = condition_csd_local( ...
        pl_all,0,cfg_freq,cfg_dics.context_weighting);

    assert(isequal(freq_ot_conf.label,freq_ot_cong.label), ...
        'OT conflict/congruent channel order differs.');
    assert(isequal(freq_pl_conf.label,freq_pl_cong.label), ...
        'PL conflict/congruent channel order differs.');
    assert(numel(freq_ot_conf.label)==204 && numel(freq_pl_conf.label)==204, ...
        'Frequency data must contain exactly 204 MEGGRAD channels per day.');

    ot_freq_conf{sub} = getPowspctrmCombine_local(freq_ot_conf,C_ot_conf);
    ot_freq_cong{sub} = getPowspctrmCombine_local(freq_ot_cong,C_ot_cong);
    pl_freq_conf{sub} = getPowspctrmCombine_local(freq_pl_conf,C_pl_conf);
    pl_freq_cong{sub} = getPowspctrmCombine_local(freq_pl_cong,C_pl_cong);

    % Equal conflict/congruent weights in the common filter within each day.
    C_ot_common = hermitian_local((C_ot_conf+C_ot_cong)/2);
    C_pl_common = hermitian_local((C_pl_conf+C_pl_cong)/2);
    C_global = blkdiag(C_ot_common,C_pl_common);

    result = dics_reconstruction_matched_power_local( ...
        C_global,C_ot_conf,C_pl_conf,C_ot_cong,C_pl_cong, ...
        ot_lf,pl_lf,freq_ot_conf.label,freq_pl_conf.label,cfg_dics);
    assert(numel(result.daycal_raw.conf)==size(sourcemodel.pos,1), ...
        'Leadfield/source-model grid count mismatch for %s.',sub_name);
    assert(result.quality.n_valid>0, ...
        'No valid source grid survived DICS estimation for %s.',sub_name);

    source_conf = sourcemodel;
    source_cong = sourcemodel;

    % Always retain both estimators obtained from the filters selected by
    % cfg_dics.apply_day_unit_gain. These generic names remain correct when
    % day-specific unit-gain calibration is disabled.
    source_conf.avg.pow_selectedfilter_raw = result.daycal_raw.conf;
    source_cong.avg.pow_selectedfilter_raw = result.daycal_raw.cong;
    source_conf.avg.pow_selectedfilter_unitnoisegain = result.daycal_ung.conf;
    source_cong.avg.pow_selectedfilter_unitnoisegain = result.daycal_ung.cong;
    source_conf.avg.pow_ot_selectedfilter_raw = result.day.ot_raw.conf;
    source_conf.avg.pow_pl_selectedfilter_raw = result.day.pl_raw.conf;
    source_cong.avg.pow_ot_selectedfilter_raw = result.day.ot_raw.cong;
    source_cong.avg.pow_pl_selectedfilter_raw = result.day.pl_raw.cong;
    source_conf.avg.pow_ot_selectedfilter_unitnoisegain = result.day.ot_ung.conf;
    source_conf.avg.pow_pl_selectedfilter_unitnoisegain = result.day.pl_ung.conf;
    source_cong.avg.pow_ot_selectedfilter_unitnoisegain = result.day.ot_ung.cong;
    source_cong.avg.pow_pl_selectedfilter_unitnoisegain = result.day.pl_ung.cong;

    % Legacy aliases retained for compatibility with earlier result files.
    % Consult source_*.cfg_dics.apply_day_unit_gain before interpreting the
    % word "daycal" when calibration is disabled.
    source_conf.avg.pow_daycal_raw = result.daycal_raw.conf;
    source_cong.avg.pow_daycal_raw = result.daycal_raw.cong;
    source_conf.avg.pow_daycal_unitnoisegain = result.daycal_ung.conf;
    source_cong.avg.pow_daycal_unitnoisegain = result.daycal_ung.cong;

    % Also retain day-specific powers from the selected OT/PL filters.
    source_conf.avg.pow_ot_daycal_raw = result.day.ot_raw.conf;
    source_conf.avg.pow_pl_daycal_raw = result.day.pl_raw.conf;
    source_cong.avg.pow_ot_daycal_raw = result.day.ot_raw.cong;
    source_cong.avg.pow_pl_daycal_raw = result.day.pl_raw.cong;

    if cfg_dics.use_unit_white_noise_gain
        source_conf.avg.pow = result.daycal_ung.conf;
        source_cong.avg.pow = result.daycal_ung.cong;
        if cfg_dics.apply_day_unit_gain
            normalization_label = ...
                'day-unit-gain + unit-white-noise-gain';
        else
            normalization_label = ...
                'no-day-unit-gain + unit-white-noise-gain';
        end
    else
        source_conf.avg.pow = result.daycal_raw.conf;
        source_cong.avg.pow = result.daycal_raw.cong;
        if cfg_dics.apply_day_unit_gain
            normalization_label = 'day-unit-gain raw power';
        else
            normalization_label = 'no-day-unit-gain raw power';
        end
    end
    source_conf.power_normalization = normalization_label;
    source_cong.power_normalization = normalization_label;

    source_conf.unit = 'mm';
    source_cong.unit = 'mm';
    source_conf.cfg_dics = cfg_dics;
    source_cong.cfg_dics = cfg_dics;

    all_source_conf{sub} = source_conf;
    all_source_cong{sub} = source_cong;
    all_source_quality{sub} = result.quality;

    fprintf('Valid source grids (%s): %d/%d\n', ...
        normalization_label, ...
        result.quality.n_valid,numel(result.quality.valid));
end

% ------------------------------- Save ---------------------------------
save(fullfile(save_source_path,'all_frequency.mat'), ...
    'ot_freq_conf','pl_freq_conf','ot_freq_cong','pl_freq_cong','-v7.3');
save(fullfile(save_source_path,'all_source.mat'), ...
    'all_source_conf','all_source_cong','all_source_quality', ...
    'subject_ids','cfg_target','cfg_dics','cfg_freq','-v7.3');

fprintf('\nSaved results to:\n%s\n',save_source_path);
%% ============================= Functions ===============================
function [freq_pooled,C] = condition_csd_local( ...
    data_all,conflict_code,cfg_freq,weighting)

idx = find(data_all.trialinfo(:,1)==conflict_code);
assert(~isempty(idx),'Requested conflict/congruent cell contains no trials.');
cfg_sel = [];
cfg_sel.trials = idx;
data_pool = ft_selectdata(cfg_sel,data_all);
freq_pooled = ft_freqanalysis(cfg_freq,data_pool);

switch weighting
    case 'pooled_trials'
        C = calculate_csd_local(freq_pooled);

    case 'equal_context'
        context_values = unique(data_pool.trialinfo(:,3));
        assert(numel(context_values)==2, ...
            'equal_context requires both social and nonsocial trials.');
        C_context = cell(numel(context_values),1);
        for k = 1:numel(context_values)
            cfg_ctx = [];
            cfg_ctx.trials = find(data_pool.trialinfo(:,3)==context_values(k));
            assert(~isempty(cfg_ctx.trials),'An expected context cell is empty.');
            data_ctx = ft_selectdata(cfg_ctx,data_pool);
            freq_ctx = ft_freqanalysis(cfg_freq,data_ctx);
            assert(isequal(freq_ctx.label,freq_pooled.label), ...
                'Context-specific channel order differs from pooled order.');
            C_context{k} = calculate_csd_local(freq_ctx);
        end
        C = mean(cat(3,C_context{:}),3);

    otherwise
        error('Unknown context weighting: %s',weighting);
end
C = hermitian_local(C);
end

function C = calculate_csd_local(freq_fourier)
X = freq_fourier.fourierspctrm(:,:,1);
X = reshape(X,size(X,1),size(X,2));
assert(all(isfinite(X(:))),'Fourier coefficients contain NaN/Inf.');
C = (X'*X)/size(X,1);
C = hermitian_local(C);
end

function result = dics_reconstruction_matched_power_local( ...
    C_global,C_ot_conf,C_pl_conf,C_ot_cong,C_pl_cong, ...
    lf_ot,lf_pl,label_ot,label_pl,cfg)
% Match the grid-filter definition used before ROI filter-SVD reconstruction.
% Each grid keeps the original three Cartesian leadfield columns. Because
% MEG usually observes a rank-2 subspace, pinv is used for the 3x3 filter
% denominator and for optional day-specific gain calibration.

label_ot = label_ot(:);
label_pl = label_pl(:);
n_ot = numel(label_ot);
n_pl = numel(label_pl);
assert(n_ot==204 && n_pl==204,'Expected 204 channels per day.');

[tf_ot,idx_ot] = ismember(label_ot,lf_ot.label(:));
[tf_pl,idx_pl] = ismember(label_pl,lf_pl.label(:));
assert(all(tf_ot) && all(tf_pl),'CSD label absent from leadfield labels.');
assert(numel(unique(idx_ot))==n_ot && numel(unique(idx_pl))==n_pl, ...
    'Channel-to-leadfield mapping is not one-to-one.');
assert(numel(lf_ot.leadfield)==numel(lf_pl.leadfield), ...
    'OT/PL leadfields have different grid counts.');
if isfield(lf_ot,'pos') && isfield(lf_pl,'pos')
    assert(isequal(size(lf_ot.pos),size(lf_pl.pos)), ...
        'OT/PL source-grid arrays differ in size.');
    pos_error = max(vecnorm(double(lf_ot.pos)-double(lf_pl.pos),2,2), ...
        [],'omitnan');
    assert(pos_error<1e-4, ...
        'OT/PL leadfields do not use the same source-grid geometry.');
end

C_global = hermitian_local(C_global);
C_ot_conf = hermitian_local(C_ot_conf);
C_pl_conf = hermitian_local(C_pl_conf);
C_ot_cong = hermitian_local(C_ot_cong);
C_pl_cong = hermitian_local(C_pl_cong);

scale = real(trace(C_global))/(n_ot+n_pl);
assert(isfinite(scale) && scale>0,'Invalid global CSD scale.');
C_reg = hermitian_local(C_global + ...
    cfg.lambda*scale*eye(n_ot+n_pl));

[~,chol_flag] = chol(C_reg);
if chol_flag==0
    C_solver = decomposition(C_reg,'chol');
else
    warning('Regularized CSD is not positive definite; using LU.');
    C_solver = decomposition(C_reg,'lu');
end

n_grid = numel(lf_ot.leadfield);
names = {'conf','cong'};
for k = 1:numel(names)
    result.daycal_raw.(names{k}) = nan(n_grid,1);
    result.daycal_ung.(names{k}) = nan(n_grid,1);
    result.day.ot_raw.(names{k}) = nan(n_grid,1);
    result.day.pl_raw.(names{k}) = nan(n_grid,1);
    result.day.ot_ung.(names{k}) = nan(n_grid,1);
    result.day.pl_ung.(names{k}) = nan(n_grid,1);
end

quality.valid = false(n_grid,1);
quality.effective_filter_rank = nan(n_grid,1);
quality.rcond_filter_denominator = nan(n_grid,1);
quality.joint_gain_error = nan(n_grid,1);
quality.ot_gain_error_before = nan(n_grid,1);
quality.pl_gain_error_before = nan(n_grid,1);
quality.ot_gain_error_after = nan(n_grid,1);
quality.pl_gain_error_after = nan(n_grid,1);
quality.ot_unit_gain_error = nan(n_grid,1);
quality.pl_unit_gain_error = nan(n_grid,1);
quality.ot_white_noise_gain = nan(n_grid,1);
quality.pl_white_noise_gain = nan(n_grid,1);
quality.apply_day_unit_gain = cfg.apply_day_unit_gain;
quality.use_unit_white_noise_gain = cfg.use_unit_white_noise_gain;

for g = 1:n_grid
    L1 = lf_ot.leadfield{g};
    L2 = lf_pl.leadfield{g};
    if isempty(L1) || isempty(L2), continue; end

    L1 = double(L1(idx_ot,:));
    L2 = double(L2(idx_pl,:));
    if size(L1,2)~=size(L2,2) || any(~isfinite([L1(:);L2(:)]))
        continue
    end

    Ljoint = [L1;L2];
    Z = C_solver\Ljoint;
    D = hermitian_local(Ljoint'*Z);
    singular_values = svd(D);
    if isempty(singular_values) || singular_values(1)<=0 || ...
            any(~isfinite(singular_values))
        continue
    end
    rank_tol = cfg.gain_pinv_rtol*singular_values(1);
    effective_rank = sum(singular_values>rank_tol);
    if effective_rank<1, continue; end

    % This is the same three-row joint filter used in the reconstruction
    % code before filters from all selected ROI grids are stacked for SVD.
    Wjoint = pinv(D)*Z';
    if any(~isfinite(Wjoint(:))), continue; end
    W1_raw = Wjoint(:,1:n_ot);
    W2_raw = Wjoint(:,n_ot+1:end);

    [~,Qjoint] = calibrate_day_unit_gain_local( ...
        Wjoint,Ljoint,false,cfg.gain_pinv_rtol);
    [W1_used,Q1] = calibrate_day_unit_gain_local( ...
        W1_raw,L1,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);
    [W2_used,Q2] = calibrate_day_unit_gain_local( ...
        W2_raw,L2,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);

    q1 = real(trace(W1_used*W1_used'));
    q2 = real(trace(W2_used*W2_used'));
    if ~isfinite(q1) || q1<=0 || ~isfinite(q2) || q2<=0
        continue
    end

    [p1_conf,a1] = source_power_local(W1_used,C_ot_conf);
    [p2_conf,a2] = source_power_local(W2_used,C_pl_conf);
    [p1_cong,a3] = source_power_local(W1_used,C_ot_cong);
    [p2_cong,a4] = source_power_local(W2_used,C_pl_cong);
    if ~(a1 && a2 && a3 && a4), continue; end

    result.day.ot_raw.conf(g) = p1_conf;
    result.day.pl_raw.conf(g) = p2_conf;
    result.day.ot_raw.cong(g) = p1_cong;
    result.day.pl_raw.cong(g) = p2_cong;
    result.day.ot_ung.conf(g) = p1_conf/q1;
    result.day.pl_ung.conf(g) = p2_conf/q2;
    result.day.ot_ung.cong(g) = p1_cong/q1;
    result.day.pl_ung.cong(g) = p2_cong/q2;

    % Equal day weights reproduce the intended two-session estimator
    % without allowing the day with larger absolute gain to dominate.
    result.daycal_raw.conf(g) = mean([p1_conf,p2_conf]);
    result.daycal_raw.cong(g) = mean([p1_cong,p2_cong]);
    result.daycal_ung.conf(g) = mean([p1_conf/q1,p2_conf/q2]);
    result.daycal_ung.cong(g) = mean([p1_cong/q1,p2_cong/q2]);

    quality.valid(g) = true;
    quality.effective_filter_rank(g) = effective_rank;
    quality.rcond_filter_denominator(g) = rcond(D);
    quality.joint_gain_error(g) = Qjoint.gain_error_before;
    quality.ot_gain_error_before(g) = Q1.gain_error_before;
    quality.pl_gain_error_before(g) = Q2.gain_error_before;
    quality.ot_gain_error_after(g) = Q1.gain_error_after;
    quality.pl_gain_error_after(g) = Q2.gain_error_after;
    quality.ot_unit_gain_error(g) = Q1.gain_error_after;
    quality.pl_unit_gain_error(g) = Q2.gain_error_after;
    quality.ot_white_noise_gain(g) = q1;
    quality.pl_white_noise_gain(g) = q2;
end

quality.n_valid = sum(quality.valid);
quality.n_grid = n_grid;
quality.regularization_scale = scale;
quality.filter_definition = ...
    'three-row reconstruction-matched joint DICS filter before ROI filter-SVD';
if cfg.apply_day_unit_gain
    result.filter_mode = 'day-unit-gain';
else
    result.filter_mode = 'no-day-unit-gain';
end
result.quality = quality;
end


function [Wcal,Q] = calibrate_day_unit_gain_local(W,L,do_apply,rtol)
G = W*L;
if any(~isfinite(G(:)))
    Wcal = nan(size(W));
    Q = struct('gain_error_before',nan,'gain_error_after',nan);
    return
end
[U,S,V] = svd(G,'econ');
s = real(diag(S));
if isempty(s) || ~isfinite(s(1)) || s(1)<=0
    Wcal = nan(size(W));
    Q = struct('gain_error_before',nan,'gain_error_after',nan);
    return
end
keep = s>rtol*s(1);
if ~any(keep)
    Wcal = nan(size(W));
    Q = struct('gain_error_before',nan,'gain_error_after',nan);
    return
end
Gpinv = V(:,keep)*diag(1./s(keep))*U(:,keep)';
Pobs = Gpinv*G;
if do_apply
    Wcal = Gpinv*W;
else
    Wcal = W;
end
Q = struct();
Q.applied = logical(do_apply);
Q.supported_rank = sum(keep);
Q.singular_values = s;
Q.relative_tolerance = rtol;
Q.gain_error_before = norm(G-Pobs,'fro')/ ...
    max(norm(Pobs,'fro'),realmin('double'));
Q.gain_error_after = norm(Wcal*L-Pobs,'fro')/ ...
    max(norm(Pobs,'fro'),realmin('double'));
Q.white_noise_gain_before = real(trace(W*W'));
Q.white_noise_gain_after = real(trace(Wcal*Wcal'));
end

function assert_source_grid_geometry_local(leadfield,template_grid,name)
% Check output-array length only, Do not compare individual coordinates with MNI template coordinates

n_lf   = numel(leadfield.leadfield);
n_grid = size(template_grid.pos,1);

assert(n_lf == n_grid, ...
    ['%s leadfield contains %d grid points, but the output template ' ...
     'contains %d grid points.'], ...
    name,n_lf,n_grid);
end

function [p,is_valid] = source_power_local(W,C)
S = hermitian_local(W*C*W');
p = real(trace(S));
tol = 1e-10*max(norm(S,'fro'),1);
is_valid = isfinite(p) && p>0 && min(real(eig(S)))>=-tol;
if ~is_valid, p=nan; end
end

function C = hermitian_local(C)
C = (C+C')/2;
end

function out = getPowspctrmCombine_local(freq_fourier,CSD_matrix)
tmp = [];
tmp.label = freq_fourier.label;
tmp.freq = freq_fourier.freq;
tmp.dimord = 'chan_freq';
tmp.powspctrm = real(diag(CSD_matrix));
if isfield(freq_fourier,'grad'), tmp.grad=freq_fourier.grad; end
cfg_cmb = [];
cfg_cmb.method = 'sum';
out = ft_combineplanar(cfg_cmb,tmp);
end
