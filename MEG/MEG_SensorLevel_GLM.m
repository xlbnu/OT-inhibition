%% Prepare sensor-level planar-combined, ITI-baseline-corrected data
% ============================== Config ===================================
cfg_main = struct();
cfg_main.num_subjects = 19;

cfg_main.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(cfg_main.ft_path);
ft_defaults;

cfg_main.analysis_window = [-7.2 1];
cfg_main.target_window = [-2 1];

cfg_main.original_fsample = 250;
cfg_main.resamplefs = 100;
cfg_main.expected_fsample = cfg_main.resamplefs;

cfg_main.baseline_duration = 0.5; % seconds
cfg_main.min_baseline_samples = 10;

% 'rms'       = sqrt((g1.^2 + g2.^2)/2), calculated manually;
% 'magnitude' = FieldTrip ft_combineplanar(cfg.method='sum'), which for
%               raw/time-domain data returns sqrt(g1.^2 + g2.^2).
cfg_main.planar_combine_method = 'magnitude';

cfg_main.base_data_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton_T(-7_1)s';
cfg_main.base_behavior_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\behavior\analysis_data\behaviour_meg.mat';
cfg_main.behavior_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\pressPutton102(250hz)\data_out_T(-7_1)s.mat';

cfg_main.output_dir = ...
    'F:\xianliang\exp_data\OT_N2_MEG\sensor_level_ERF\planar_pressPutton_NoBaseline_100hz_T(-7_1)s';

cfg_main.condition_names = {'social_OT', 'social_PL', 'nonsocial_OT', 'nonsocial_PL'};
cfg_main.condition_dirs = {
    fullfile(cfg_main.base_data_path, 'social_OT_306(250hz)')
    fullfile(cfg_main.base_data_path, 'social_PL_306(250hz)')
    fullfile(cfg_main.base_data_path, 'nonsocial_OT_306(250hz)')
    fullfile(cfg_main.base_data_path, 'nonsocial_PL_306(250hz)')
    };

cfg_main.trialinfo_label = { ...
    'conflict','ischange','lo_s1','lo_o1','condition','run_id', ...
    'lo_chosen','lo_unchosen','lo_plus','lo_minus','lo_cplus','lo_cminus', ...
    'rtime2','osti2','oiti','lo_s2','ctime2','oconf2','answer1','answer2','osti'};

if ~exist(cfg_main.output_dir, 'dir')
    mkdir(cfg_main.output_dir);
end

% ============================== Load fixed data ===========================
data_tmp = load(cfg_main.base_behavior_file);
data_all = {data_tmp.ots_meg; data_tmp.pls_meg; data_tmp.otn_meg; data_tmp.pln_meg};

if ~exist('data_out', 'var')
    tmp = load(cfg_main.behavior_file, 'data_out');
    data_out = tmp.data_out;
end

all_subject_ids = cell(cfg_main.num_subjects, 1);
all_sensor_file = cell(cfg_main.num_subjects, 1);
all_quality = cell(cfg_main.num_subjects, 1);

% ============================== Main loop =================================
for i_sub = 1:cfg_main.num_subjects
    subject_id = normalize_subject_id(data_all{1}(i_sub).subsName);
    all_subject_ids{i_sub} = subject_id;

    fprintf('\n====================================================\n');
    fprintf('[%d/%d] Subject: %s\n', i_sub, cfg_main.num_subjects, subject_id);

    behavior_row = find(strcmp(data_out(:,1), subject_id));
    assert(isscalar(behavior_row), ...
        'data_out subject %s is not uniquely matched.', subject_id);

    all_run_data = load_subject_runs( ...
        i_sub, subject_id, behavior_row, data_out, data_all, cfg_main);
    all_run_data = all_run_data(~cellfun(@isempty, all_run_data));
    assert(~isempty(all_run_data), 'Subject %s has no usable runs.', subject_id);

    cfg_append = [];
    cfg_append.keepsampleinfo = 'no';
    data_allrun_250 = ft_appenddata(cfg_append, all_run_data{:});

    bhv = collect_run_bhv_table(all_run_data, cfg_main.trialinfo_label);
    assert(height(bhv) == numel(data_allrun_250.trial), ...
        'Behavior trial number does not match MEG trial number before resampling.');

    % Resample first, then combine planar and baseline-correct.
    cfg_rs = [];
    cfg_rs.resamplefs = cfg_main.resamplefs;
    cfg_rs.demean = 'no';
    cfg_rs.detrend = 'no';
    data_allrun = ft_resampledata(cfg_rs, data_allrun_250);

    assert(numel(data_allrun.trial) == height(bhv), ...
        'Behavior trial number does not match MEG trial number after resampling.');

    [sensor_data, sensor_time, sensor_label, sensor_baseline_mean, ...
        base_data, baseline_time, baseline_info, quality] = ...
        combine_planar_and_subtract_iti_baseline(data_allrun, bhv, cfg_main);

    valid_trial = baseline_info.valid_baseline;
    sensor_data = sensor_data(valid_trial,:,:);
    sensor_bhv = bhv(valid_trial,:);
    sensor_baseline_mean = sensor_baseline_mean(valid_trial,:);
    base_data = base_data(valid_trial,:,:);
    baseline_info = baseline_info(valid_trial,:);

    assert(size(sensor_data,1) == height(sensor_bhv) && ...
        size(base_data,1) == height(sensor_bhv), ...
        'Saved task, baseline, and behavior trial counts differ.');
    assert(size(base_data,2) == numel(sensor_label) && ...
        size(base_data,3) == numel(baseline_time), ...
        'Saved baseline data dimensions differ from labels/time.');

    subject_output = fullfile(cfg_main.output_dir, sprintf( ...
        '%s_sensor_planar_100Hz_noBaseline.mat', subject_id));

    save(subject_output, ...
        'sensor_data', 'sensor_time', 'sensor_label', 'sensor_bhv', ...
        'sensor_baseline_mean', 'base_data', 'baseline_time', ...
        'baseline_info', 'quality', 'cfg_main', 'subject_id', '-v7.3');

    all_sensor_file{i_sub} = subject_output;
    all_quality{i_sub} = quality;

    fprintf('Saved: %s\n', subject_output);
    fprintf('  trials kept after baseline check: %d/%d\n', ...
        size(sensor_data,1), height(bhv));
end


%% Sensor-level GLM
clear;
cfg_glm = struct();
cfg_glm.sensor_data_dir = 'F:\xianliang\exp_data\OT_N2_MEG\sensor_level_ERF\planar_pressPutton_NoBaseline_100hz_v2';
cfg_glm.output_file = 'F:\xianliang\exp_data\OT_N2_MEG\sensor_level_ERF\Sensor_vchan_GLM\pressPutton_task_and_ITI_glm_beta_OT_PL_noBaseline_labelAligned_v2.mat';
cfg_glm.file_pattern = '*_sensor_planar_100Hz_noBaseline.mat';

% Original condition codes:
%   1=social_OT, 2=social_PL, 3=nonsocial_OT, 4=nonsocial_PL.
% The following setting pools social and nonsocial within drug day.
cfg_glm.condition_names = {'OT','PL'};
cfg_glm.condition_codes = {[1, 3] ,[2 ,4]};
cfg_glm.min_trial = 10;
cfg_glm.iti_fsample = 100;

% 'complete' (recommended): use only trials with all ITI time samples, so
% every ITI time point is fitted with the same trials.
% 'available_at_each_time': retain partial ITIs; NaN trials are omitted at
% each time point, so sample composition can change over ITI time.
cfg_glm.baseline_trial_policy = 'complete';

files = dir(fullfile(cfg_glm.sensor_data_dir, cfg_glm.file_pattern));
assert(numel(files) == 19, 'Expected 19 subject files, found %d.', numel(files));

S0 = load(fullfile(files(1).folder, files(1).name), ...
    'sensor_label', 'sensor_time', 'baseline_time', 'base_data');
required0 = {'sensor_label','sensor_time','base_data'};
for k = 1:numel(required0)
    assert(isfield(S0,required0{k}), 'First file lacks variable: %s.', required0{k});
end

label_ref = S0.sensor_label(:);
task_time_ref = S0.sensor_time(:)';
n_baseline_time = size(S0.base_data,3);
baseline_time_ref = (0:(n_baseline_time-1)) / cfg_glm.iti_fsample;
if isfield(S0,'baseline_time')
    assert(numel(S0.baseline_time)==n_baseline_time && ...
        max(abs(S0.baseline_time(:)'-baseline_time_ref))<1e-10, ...
        'Saved ITI time axis differs from the configured 100-Hz 0-based axis.');
end

n_condition = numel(cfg_glm.condition_names);
n_subject = numel(files);
n_channel = numel(label_ref);
n_task_time = numel(task_time_ref);
n_beta = 3;
n_model = 5;

task_size = [n_condition, n_subject, n_channel, n_task_time, n_beta];
base_size = [n_condition, n_subject, n_channel, n_baseline_time, n_beta];

efr_glm = struct();
efr_glm.beta_stay = nan(task_size, 'single');
efr_glm.beta_change = nan(task_size, 'single');
efr_glm.beta_congruent = nan(task_size, 'single');
efr_glm.beta_chosen = nan(task_size, 'single');
efr_glm.beta_bias = nan(task_size, 'single');

efr_glm.beta_stay_base = nan(base_size, 'single');
efr_glm.beta_change_base = nan(base_size, 'single');
efr_glm.beta_congruent_base = nan(base_size, 'single');
efr_glm.beta_chosen_base = nan(base_size, 'single');
efr_glm.beta_bias_base = nan(base_size, 'single');

efr_glm.trial_count_task = zeros(n_condition,n_subject,n_model,'uint16');
efr_glm.trial_count_base = zeros(n_condition,n_subject,n_model,'uint16');
efr_glm.trial_count_label = {'conflict_stay','conflict_change', ...
    'congruent_stay','all_conflict_chosen','all_conflict_bias'};
efr_glm.subject_id = cell(n_subject,1);
efr_glm.subject_file = cell(n_subject,1);

for i_sub = 1:n_subject
    file = fullfile(files(i_sub).folder,files(i_sub).name);
    S = load(file, 'sensor_data','sensor_time','sensor_label', ...
        'base_data','baseline_time','baseline_info','sensor_bhv','subject_id');

    required = {'sensor_data','sensor_time','sensor_label','base_data', ...
        'sensor_bhv','subject_id'};
    for k = 1:numel(required)
        assert(isfield(S,required{k}), 'File %s lacks variable: %s.', ...
            files(i_sub).name,required{k});
    end

    n_trial = height(S.sensor_bhv);
    assert(size(S.sensor_data,1)==n_trial && size(S.base_data,1)==n_trial, ...
        'Task/ITI/behavior trial counts differ: %s.', files(i_sub).name);
    assert(size(S.sensor_data,2)==numel(S.sensor_label) && ...
        size(S.base_data,2)==numel(S.sensor_label), ...
        'Task/ITI channel dimensions differ from labels: %s.',files(i_sub).name);
    assert(size(S.sensor_data,3)==numel(S.sensor_time), ...
        'Task data/time dimensions differ: %s.',files(i_sub).name);
    assert(size(S.base_data,3)==n_baseline_time, ...
        'ITI length differs across subjects: %s.',files(i_sub).name);
    assert(numel(S.sensor_time)==n_task_time && ...
        max(abs(S.sensor_time(:)'-task_time_ref))<1e-10, ...
        'Task time axis differs across subjects: %s.',files(i_sub).name);
    if isfield(S,'baseline_time')
        assert(numel(S.baseline_time)==n_baseline_time && ...
            max(abs(S.baseline_time(:)'-baseline_time_ref))<1e-10, ...
            'ITI time axis differs from the configured 100-Hz axis: %s.', ...
            files(i_sub).name);
    end

    [present,loc_current] = ismember(label_ref,S.sensor_label(:));
    assert(all(present) && numel(S.sensor_label)==n_channel, ...
        'Channel set differs across subjects: %s.',files(i_sub).name);
    sensor_data = S.sensor_data(:,loc_current,:);
    base_data = S.base_data(:,loc_current,:);
    bhv = S.sensor_bhv;

    assert(all(ismember(unique(bhv.condition)',1:4)), ...
        'Invalid original condition code: %s.',files(i_sub).name);
    assert(all(ismember(unique(bhv.run_id)',1:3)), ...
        'Invalid run_id: %s.',files(i_sub).name);

    % For time-resolved ITI GLM, a constant trial set across time is the
    % cleanest choice. Partial ITI signals remain saved in the source files.
    switch lower(cfg_glm.baseline_trial_policy)
        case 'complete'
            complete_base = all(isfinite(reshape(base_data,n_trial,[])),2);
        case 'available_at_each_time'
            complete_base = true(n_trial,1);
        otherwise
            error('Unknown baseline_trial_policy: %s.',cfg_glm.baseline_trial_policy);
    end

    for i_condition = 1:n_condition
        in_condition = ismember(bhv.condition,cfg_glm.condition_codes{i_condition});
        masks_task = make_model_masks(bhv,in_condition);
        masks_base = masks_task;
        for i_model = 1:n_model
            masks_base{i_model} = masks_base{i_model} & complete_base;
        end

        efr_glm.trial_count_task(i_condition,i_sub,:) = ...
            reshape(uint16(cellfun(@nnz,masks_task)),1,1,[]);
        efr_glm.trial_count_base(i_condition,i_sub,:) = ...
            reshape(uint16(cellfun(@nnz,masks_base)),1,1,[]);

        task_beta = fit_five_glms(sensor_data,bhv,masks_task,cfg_glm.min_trial);
        base_beta = fit_five_glms(base_data,bhv,masks_base,cfg_glm.min_trial);

        efr_glm.beta_stay(i_condition,i_sub,:,:,:) = task_beta{1};
        efr_glm.beta_change(i_condition,i_sub,:,:,:) = task_beta{2};
        efr_glm.beta_congruent(i_condition,i_sub,:,:,:) = task_beta{3};
        efr_glm.beta_chosen(i_condition,i_sub,:,:,:) = task_beta{4};
        efr_glm.beta_bias(i_condition,i_sub,:,:,:) = task_beta{5};

        efr_glm.beta_stay_base(i_condition,i_sub,:,:,:) = base_beta{1};
        efr_glm.beta_change_base(i_condition,i_sub,:,:,:) = base_beta{2};
        efr_glm.beta_congruent_base(i_condition,i_sub,:,:,:) = base_beta{3};
        efr_glm.beta_chosen_base(i_condition,i_sub,:,:,:) = base_beta{4};
        efr_glm.beta_bias_base(i_condition,i_sub,:,:,:) = base_beta{5};
    end

    efr_glm.subject_id{i_sub} = char(string(S.subject_id));
    efr_glm.subject_file{i_sub} = file;
    fprintf('Subject %d/%d finished: %s\n', ...
        i_sub,n_subject,efr_glm.subject_id{i_sub});
end

efr_glm.sensor_time = task_time_ref;
efr_glm.baseline_time = baseline_time_ref;
efr_glm.iti_time = baseline_time_ref;
efr_glm.sensor_label = label_ref;
efr_glm.condition_names = cfg_glm.condition_names;
efr_glm.condition_codes = cfg_glm.condition_codes;
efr_glm.beta_label_self_other = {'intercept','self','other'};
efr_glm.beta_label_chosen = {'intercept','chosen','unchosen'};
efr_glm.beta_label_bias = {'intercept','cplus','cminus'};
efr_glm.cfg = cfg_glm;

output_dir = fileparts(cfg_glm.output_file);
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(cfg_glm.output_file,'efr_glm','-v7.3');










%% ============================== Functions =================================
function all_run_data = load_subject_runs( ...
    i_sub, subject_id, behavior_row, data_out, data_all, cfg_main)

all_run_data = cell(1, 12);
run_counter = 0;

for i_condition = 1:4
    bhv_bad = data_out{behavior_row, i_condition+1};
    bhv0 = data_all{i_condition}(i_sub);

    assert(isequal(normalize_subject_id(bhv0.subsName), subject_id), ...
        'Behavior subject order mismatch: %s, condition %s.', ...
        subject_id, cfg_main.condition_names{i_condition});

    meg_file = find_unique_subject_file(cfg_main.condition_dirs{i_condition}, subject_id);
    data_process = load_meg_structure(meg_file);

    bhv0.run_id = [ones(80,1); 2*ones(80,1); 3*ones(80,1)];
    bhv0.condition = i_condition * ones(numel(bhv0.answer1), 1);
    bhv0 = remove_bad_trials_from_behavior(bhv0, bhv_bad.bad_trials);

    cfg_crop = [];
    cfg_crop.toilim = cfg_main.analysis_window;
    data_long = ft_redefinetrial(cfg_crop, data_process);

    for i_run = 1:3
        run_counter = run_counter + 1;

        % Keep consistent with source/LCMV scripts: remove congruent-change trials.
        selected_trials = find( ...
            bhv0.run_id == i_run & ~(bhv0.conflict == 0 & bhv0.ischange == 1));

        if isempty(selected_trials)
            warning('Subject %s, %s run%d has no usable trials.', ...
                subject_id, cfg_main.condition_names{i_condition}, i_run);
            continue;
        end

        assert(max(selected_trials) <= numel(data_long.trial), ...
            'Behavior trial index exceeds MEG trial number: %s %s run%d.', ...
            subject_id, cfg_main.condition_names{i_condition}, i_run);

        cfg_select = [];
        cfg_select.trials = selected_trials;
        cfg_select.channel = 'MEGGRAD';
        run_data = ft_selectdata(cfg_select, data_long);

        % Manual planar pairing does not require grad. Removing these fields also
        % avoids appenddata warnings/errors from different head positions.
        run_data = rmfield_if_exists(run_data, {'grad','elec','hdr','sampleinfo'});
        run_data.bhv_info = select_behavior_fields(bhv0, selected_trials, cfg_main.trialinfo_label);

        all_run_data{run_counter} = run_data;
    end
end

assert(run_counter == 12, 'Internal run count error.');
end



function [sensor_data, sensor_time, pair_label, sensor_baseline_mean, ...
    base_data, baseline_time, baseline_info, quality] = ...
    combine_planar_and_subtract_iti_baseline(data_allrun, bhv, cfg_main)

time_full = data_allrun.time{2};

fsample = data_allrun.fsample;
if isempty(fsample) || ~isfinite(fsample)
    fsample = 1 / median(diff(time_full));
end

assert(abs(fsample - cfg_main.expected_fsample) < 1, ...
    'Data fsample %.3f Hz differs from expected %.3f Hz.', ...
    fsample, cfg_main.expected_fsample);

target_idx = time_full >= cfg_main.target_window(1) & ...
             time_full <= cfg_main.target_window(2);
sensor_time = time_full(target_idx);

[pair_idx, pair_label] = make_neuromag_planar_pairs(data_allrun.label);
n_trial = numel(data_allrun.trial);
n_pair = size(pair_idx, 1);
n_time = numel(sensor_time);

% ft_combineplanar is designed to operate on the complete FieldTrip raw
% structure. Run it once here (rather than once per trial) when the official
% FieldTrip magnitude is requested. The output is subsequently reordered to
% the same pair order used by pair_label.
ft_planar = [];
ft_pair_order = [];
if strcmpi(cfg_main.planar_combine_method, 'magnitude')
    cfg_cp = [];
    cfg_cp.method = 'sum';       % raw data: sqrt(g1.^2 + g2.^2)
    cfg_cp.demean = 'no';
    cfg_cp.updatesens = 'no';
    cfg_cp.feedback = 'none';
    ft_planar = ft_combineplanar(cfg_cp, data_allrun);

    % FieldTrip labels are e.g. 'MEG0112+0113', whereas the existing output
    % convention in this script is 'MEG0112+MEG0113'.
    expected_ft_label = regexprep(pair_label, '\+MEG', '+');
    [pair_found, ft_pair_order] = ismember(expected_ft_label, ft_planar.label);
    assert(all(pair_found), ...
        'ft_combineplanar output is missing %d expected planar pairs.', ...
        nnz(~pair_found));
    assert(numel(ft_planar.trial) == n_trial, ...
        'ft_combineplanar changed the number of trials.');
end

sensor_data = nan(n_trial, n_pair, n_time);

% baseline_start = bhv.oiti - bhv.osti2 - bhv.rtime2;
baseline_start = -(bhv.osti2 - bhv.osti + bhv.rtime2 + cfg_main.baseline_duration);
baseline_end = baseline_start + cfg_main.baseline_duration;
baseline_start_idx = nan(n_trial, 1);
baseline_end_idx = nan(n_trial, 1);
n_baseline_sample = nan(n_trial, 1);
valid_baseline = false(n_trial, 1);
trial_time_was_padded = false(n_trial, 1);
n_existing_time_sample = nan(n_trial, 1);

expected_base_samples = round(cfg_main.baseline_duration * fsample);
min_base_samples = cfg_main.min_baseline_samples;
baseline_time = (0:(expected_base_samples-1)) / fsample;
base_data = nan(n_trial, n_pair, expected_base_samples);
sensor_baseline_mean = nan(n_trial, n_pair);
baseline_expected_start_idx = nan(n_trial, 1);
baseline_expected_end_idx = nan(n_trial, 1);
baseline_first_saved_sample = nan(n_trial, 1);
baseline_last_saved_sample = nan(n_trial, 1);

for i_trial = 1:n_trial
    x = double(data_allrun.trial{i_trial});
    this_time = data_allrun.time{i_trial};
    n_existing_time_sample(i_trial) = numel(this_time);

    assert(size(x, 2) == numel(this_time), ...
        'trial %d data length does not match its time length.', i_trial);

    % Align each trial to the first trial's time axis. If this trial is
    % shorter, missing time points remain NaN.
    [is_match, loc] = ismembertol(this_time(:), time_full(:), 1e-8);
    assert(all(is_match), ...
        'trial %d contains time points not present in the first trial time axis.', i_trial);

    trial_time_was_padded(i_trial) = numel(this_time) < numel(time_full);

    switch lower(cfg_main.planar_combine_method)
        case 'rms'
            x_full = nan(size(x, 1), numel(time_full));
            x_full(:, loc) = x;
            x_pair = combine_planar_trial_rms(x_full, pair_idx);

        case 'magnitude'
            assert(numel(ft_planar.time{i_trial}) == numel(this_time) && ...
                max(abs(ft_planar.time{i_trial}(:) - this_time(:))) < 1e-8, ...
                'ft_combineplanar time axis differs for trial %d.', i_trial);
            x_ft = double(ft_planar.trial{i_trial}(ft_pair_order, :));
            assert(size(x_ft,2) == numel(this_time), ...
                'ft_combineplanar data/time length mismatch for trial %d.', i_trial);
            x_pair = nan(n_pair, numel(time_full));
            x_pair(:, loc) = x_ft;

        otherwise
            error('Unknown planar_combine_method: %s.', ...
                cfg_main.planar_combine_method);
    end

    b_start = baseline_start(i_trial);
    if ~isfinite(b_start)
        continue;
    end

    idx_start_raw = round((b_start - time_full(1)) * fsample) + 1;
    idx_end_raw = idx_start_raw + expected_base_samples - 1;
    baseline_expected_start_idx(i_trial) = idx_start_raw;
    baseline_expected_end_idx(i_trial) = idx_end_raw;

    % Reject only when the complete requested baseline lies outside the epoch.
    if idx_start_raw > numel(time_full) || idx_end_raw < 1
        continue;
    end

    % Construct a fixed-length baseline vector. Requested samples outside the
    % epoch, or absent because this trial was shorter, remain NaN at their
    % correct relative positions rather than being shifted to the beginning.
    requested_idx = idx_start_raw:idx_end_raw;
    within_epoch = requested_idx >= 1 & requested_idx <= numel(time_full);
    source_idx = requested_idx(within_epoch);
    destination_idx = find(within_epoch);

    if isempty(source_idx)
        continue;
    end

    source_has_data = any(isfinite(x_pair(:,source_idx)), 1);
    source_idx = source_idx(source_has_data);
    destination_idx = destination_idx(source_has_data);
    if isempty(source_idx)
        continue;
    end

    this_base = nan(n_pair, expected_base_samples);
    this_base(:,destination_idx) = x_pair(:,source_idx);
    n_sample_per_pair = sum(isfinite(this_base), 2);

    % A trial is valid only if every saved planar pair has enough baseline
    % samples. This avoids accepting a trial for which only some channels exist.
    n_baseline_sample(i_trial) = min(n_sample_per_pair);
    if any(n_sample_per_pair < min_base_samples)
        continue;
    end

    baseline_start_idx(i_trial) = source_idx(1);
    baseline_end_idx(i_trial) = source_idx(end);
    baseline_first_saved_sample(i_trial) = destination_idx(1);
    baseline_last_saved_sample(i_trial) = destination_idx(end);

    base_mean = mean(this_base, 2, 'omitnan');
    sensor_baseline_mean(i_trial,:) = base_mean(:)';
    base_data(i_trial,:,:) = reshape(this_base, ...
        [1, n_pair, expected_base_samples]);
    % sensor_data(i_trial,:,:) = x_pair(:, target_idx) - base_mean;
    sensor_data(i_trial,:,:) = x_pair(:, target_idx);
    valid_baseline(i_trial) = true;
end

baseline_info = table();
baseline_info.baseline_start_sec = baseline_start;
baseline_info.baseline_end_sec = baseline_end;
baseline_info.baseline_start_idx = baseline_start_idx;
baseline_info.baseline_end_idx = baseline_end_idx;
baseline_info.baseline_expected_start_idx = baseline_expected_start_idx;
baseline_info.baseline_expected_end_idx = baseline_expected_end_idx;
baseline_info.baseline_first_saved_sample = baseline_first_saved_sample;
baseline_info.baseline_last_saved_sample = baseline_last_saved_sample;
baseline_info.n_baseline_sample = n_baseline_sample;
baseline_info.valid_baseline = valid_baseline;
baseline_info.trial_time_was_padded = trial_time_was_padded;
baseline_info.n_existing_time_sample = n_existing_time_sample;

quality = struct();
quality.fsample = fsample;
quality.n_trial_total = n_trial;
quality.n_trial_valid_baseline = sum(valid_baseline);
quality.n_trial_invalid_baseline = sum(~valid_baseline);
quality.n_trial_time_padded = sum(trial_time_was_padded);
quality.n_pair = n_pair;
quality.n_time = n_time;
quality.time_full_start = time_full(1);
quality.time_full_end = time_full(end);
quality.target_window = cfg_main.target_window;
quality.analysis_window = cfg_main.analysis_window;
quality.baseline_duration = cfg_main.baseline_duration;
quality.expected_baseline_samples = expected_base_samples;
quality.min_baseline_samples = min_base_samples;
quality.planar_combine_method = cfg_main.planar_combine_method;

if quality.n_trial_invalid_baseline > 0
    warning('%d/%d trials have invalid baseline and will be excluded.', ...
        quality.n_trial_invalid_baseline, quality.n_trial_total);
end
end


function x_pair = combine_planar_trial_rms(x, pair_idx)
n_pair = size(pair_idx, 1);
n_time = size(x, 2);
x_pair = nan(n_pair, n_time);

for i_pair = 1:n_pair
    g1 = x(pair_idx(i_pair,1), :);
    g2 = x(pair_idx(i_pair,2), :);
    x_pair(i_pair,:) = sqrt((g1.^2 + g2.^2) / 2);
end
end


function [pair_idx, pair_label] = make_neuromag_planar_pairs(labels)
labels = labels(:);
n_label = numel(labels);
used = false(n_label, 1);
pair_idx = [];
pair_label = {};

for i = 1:n_label
    if used(i)
        continue;
    end

    lab = char(labels{i});
    tok = regexp(lab, '^(MEG\d{3})([23])$', 'tokens', 'once');
    if isempty(tok)
        continue;
    end

    stem = tok{1};
    lab2 = [stem '2'];
    lab3 = [stem '3'];
    idx2 = find(strcmp(labels, lab2), 1);
    idx3 = find(strcmp(labels, lab3), 1);

    if isempty(idx2) || isempty(idx3)
        continue;
    end

    pair_idx(end+1,:) = [idx2 idx3]; %#ok<AGROW>
    pair_label{end+1,1} = [lab2 '+' lab3]; %#ok<AGROW>
    used([idx2 idx3]) = true;
end

assert(~isempty(pair_idx), 'No Neuromag planar gradiometer pairs found.');
end


function bhv = collect_run_bhv_table(run_cells, labels)
bhv = table();
for i_run = 1:numel(run_cells)
    b = run_cells{i_run}.bhv_info;
    n_trial = numel(b.(labels{1}));
    this_tbl = table();
    for j = 1:numel(labels)
        field = labels{j};
        assert(isfield(b, field), 'bhv_info lacks field: %s.', field);
        value = b.(field);
        assert(numel(value) == n_trial, ...
            'bhv_info field %s has inconsistent trial number.', field);
        this_tbl.(field) = double(value(:));
    end
    bhv = [bhv; this_tbl]; %#ok<AGROW>
end
end


function bhv0 = remove_bad_trials_from_behavior(bhv0, bad_trials)
fields0 = fieldnames(bhv0);
for i = 1:numel(fields0)
    value = bhv0.(fields0{i});
    if isnumeric(value) || islogical(value) || iscell(value)
        if numel(value) > 100 && numel(value) < 250
            value(bad_trials) = [];
            bhv0.(fields0{i}) = value;
        end
    end
end
end


function bhv_info = select_behavior_fields(bhv0, selected_trials, labels)
bhv_info = struct();
for i = 1:numel(labels)
    field = labels{i};
    assert(isfield(bhv0, field), 'Behavior data lacks field: %s.', field);
    value = bhv0.(field);
    assert(numel(value) >= max(selected_trials), ...
        'Behavior field %s is too short.', field);
    bhv_info.(field) = value(selected_trials);
end
end


function subject_id = normalize_subject_id(subs_name)
if iscell(subs_name)
    subs_name = subs_name{1};
end
subject_id = char(string(subs_name));
end


function file_path = find_unique_subject_file(folder, subject_id)
files = dir(fullfile(folder, [subject_id, '*.mat']));
assert(isscalar(files), ...
    'Folder %s subject %s should match 1 MAT file, actual %d.', ...
    folder, subject_id, numel(files));
file_path = fullfile(files(1).folder, files(1).name);
end


function data = load_meg_structure(file_path)
loaded = load(file_path);
if isfield(loaded, 'data_combined_pre')
    data = loaded.data_combined_pre;
    return;
end

names = fieldnames(loaded);
is_raw = false(size(names));
for i = 1:numel(names)
    value = loaded.(names{i});
    is_raw(i) = isstruct(value) && ...
        isfield(value, 'trial') && isfield(value, 'time') && isfield(value, 'label');
end
assert(sum(is_raw) == 1, ...
    'File %s cannot uniquely identify a FieldTrip raw structure.', file_path);
data = loaded.(names{find(is_raw,1)});
end


function s = rmfield_if_exists(s, fields)
for i = 1:numel(fields)
    if isfield(s, fields{i})
        s = rmfield(s, fields{i});
    end
end
end


%% glm functions
function masks = make_model_masks(bhv,in_condition)
masks = cell(5,1);
finite_self_other = isfinite(bhv.lo_s1) & isfinite(bhv.lo_o1);
finite_chosen = isfinite(bhv.lo_chosen) & isfinite(bhv.lo_unchosen);
finite_bias = isfinite(bhv.lo_cplus) & isfinite(bhv.lo_cminus);
masks{1} = in_condition & bhv.conflict==1 & bhv.ischange==0 & finite_self_other;
masks{2} = in_condition & bhv.conflict==1 & bhv.ischange==1 & finite_self_other;
masks{3} = in_condition & bhv.conflict==0 & bhv.ischange==0 & finite_self_other;
masks{4} = in_condition & bhv.conflict==1 & finite_chosen;
masks{5} = in_condition & bhv.conflict==1 & finite_bias;
end


function beta_all = fit_five_glms(data,bhv,masks,min_trial)
% data: trial x channel x time
n_channel = size(data,2);
n_time = size(data,3);
n_beta = 3;
beta_all = cell(5,1);
for i_model = 1:5
    beta_all{i_model} = nan(n_channel,n_time,n_beta,'single');
end

X = cell(5,1);
X{1} = [bhv.lo_s1,bhv.lo_o1];
X{2} = X{1};
X{3} = X{1};
X{4} = [bhv.lo_chosen,bhv.lo_unchosen];
X{5} = [bhv.lo_cplus,bhv.lo_cminus];

for i_channel = 1:n_channel
    y_channel = double(squeeze(data(:,i_channel,:)));
    beta1 = nan(n_time,n_beta);
    beta2 = nan(n_time,n_beta);
    beta3 = nan(n_time,n_beta);
    beta4 = nan(n_time,n_beta);
    beta5 = nan(n_time,n_beta);
    X1 = X{1}; X2 = X{2}; X3 = X{3}; X4 = X{4}; X5 = X{5};
    mask1 = masks{1}; mask2 = masks{2}; mask3 = masks{3};
    mask4 = masks{4}; mask5 = masks{5};

    parfor i_time = 1:n_time
        y = y_channel(:,i_time);
        beta1(i_time,:) = run_glm_safe_checked(y,X1,mask1,min_trial);
        beta2(i_time,:) = run_glm_safe_checked(y,X2,mask2,min_trial);
        beta3(i_time,:) = run_glm_safe_checked(y,X3,mask3,min_trial);
        beta4(i_time,:) = run_glm_safe_checked(y,X4,mask4,min_trial);
        beta5(i_time,:) = run_glm_safe_checked(y,X5,mask5,min_trial);
    end

    beta_all{1}(i_channel,:,:) = beta1;
    beta_all{2}(i_channel,:,:) = beta2;
    beta_all{3}(i_channel,:,:) = beta3;
    beta_all{4}(i_channel,:,:) = beta4;
    beta_all{5}(i_channel,:,:) = beta5;
end
end


function beta = run_glm_safe_checked(y,X,valid_base,min_trial)
beta = nan(1,size(X,2)+1);
valid = valid_base(:) & isfinite(y) & all(isfinite(X),2);
if nnz(valid)<min_trial, return; end

Xv = double(X(valid,:));
yv = double(y(valid));
if std(yv)<=eps(max(abs(yv))) || ...
        any(std(Xv,0,1)<=eps(max(abs(Xv),[],1)))
    return;
end

Xz = zscore(Xv,0,1);
% yz = zscore(yv);
yz = yv;
if rank([ones(size(Xz,1),1),Xz])<size(Xz,2)+1
    return;
end

try
    b = glmfit(Xz,yz,'normal','constant','on');
    beta = b(:)';
catch ME
    warning('run_glm_safe_checked:Failure','%s',ME.message);
end
end


