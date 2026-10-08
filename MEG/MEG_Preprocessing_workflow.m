% The Complete Preprocessing Workflow for Task-Based MEG, 
% Reconstructed Based on the Core Logic of the FLUX Pipeline

%% set path
base_input_path_ot = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\OT_data';  %Participant folder path
base_output_path_ot = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preprocess\OT_data';  %Participant folder path
[ots_names,otn_names]=getPath(base_input_path_ot);

base_input_path_pl = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\PL_data';  %Participant folder path
base_output_path_pl = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preprocess\PL_data';  %Participant folder path
[pls_names,pln_names]=getPath(base_input_path_pl);

load('E:\xianliang\matlab_m\social_decision_m\data\OT_all\BHV_MEG_MRS_data_all\bhv_raw_data\bhv_runData_MEG_19.mat','sub_ots','sub_pls','sub_otn','sub_pln');
%% lock to choice made
trigger_infor=[];
trigger_infor.number = 6;
trigger_infor.prestim = 3.5; 
trigger_infor.poststim = 3;
trigger_infor.resamplefs = 250;
trigger_infor.do_demean  = true;
trigger_infor.phaseName='pressPutton_noMean';
trigger_infor.shift_time='bh_tmp.oconf2-bh_tmp.osti2-bh_tmp.rtime2';
base_path='F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2';
MEGData_Optimal_Pipeline(sub_ots, ots_names, trigger_infor, base_path, 'social', 'OT');
MEGData_Optimal_Pipeline(sub_pls, pls_names, trigger_infor, base_path, 'social', 'PL');
MEGData_Optimal_Pipeline(sub_otn, otn_names, trigger_infor, base_path, 'nonsocial', 'OT');
MEGData_Optimal_Pipeline(sub_pln, pln_names, trigger_infor, base_path, 'nonsocial', 'PL');


%% get valid trial index
base_data_path = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
behavior_file = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\behavior\analysis_data\behaviour_meg.mat';
output_file = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\pressPutton102(250hz)\data_vdx.mat';
condition_names = {
    'social_OT_306(250hz)'
    'social_PL_306(250hz)'
    'nonsocial_OT_306(250hz)'
    'nonsocial_PL_306(250hz)'};
n_run = 3;
trials_per_run = 80;
% Load behavioral data
B = load(behavior_file, 'ots_meg','pls_meg','otn_meg','pln_meg');

data_bhv = {B.ots_meg,B.pls_meg,B.otn_meg,B.pln_meg};
n_sub = numel(data_bhv{1});

assert(all(cellfun(@numel,data_bhv)==n_sub), ...
    '四个 condition 的被试数量不一致。');

data_vdx = cell(n_sub,5);

% Align behavioral data with preprocessed MEG trial
for i_sub = 1:n_sub

    sub_name = data_bhv{1}(i_sub).subsName;
    if iscell(sub_name), sub_name = sub_name{1}; end
    sub_name = char(sub_name);

    data_vdx{i_sub,1} = sub_name;

    for i_cond = 1:4

        bhv = data_bhv{i_cond}(i_sub);

        this_name = bhv.subsName;
        if iscell(this_name), this_name = this_name{1}; end

        assert(strcmp(sub_name,char(this_name)), ...
            '被试顺序不一致：subject %d，condition %d。', ...
            i_sub,i_cond);

        % Find the current participant's MEG file
        cond_dir = fullfile(base_data_path,condition_names{i_cond});
        files = dir(fullfile(cond_dir,[sub_name '*.mat']));

        assert(isscalar(files), ...
            '%s，condition %d：找到 %d 个 MEG 文件，要求恰好一个。', ...
            sub_name,i_cond,numel(files));

        S = load(fullfile(files.folder,files.name));
        fields_main = fieldnames(S);

        % Locate the structure containing trial and rejection_info data
        is_main = cellfun(@(f) ...
            isstruct(S.(f)) && isscalar(S.(f)) && ...
            isfield(S.(f),'trial') && ...
            isfield(S.(f),'rejection_info'),fields_main);

        assert(sum(is_main)==1, '%s，condition %d：无法唯一定位 MEG 数据结构。', sub_name,i_cond);

        meg = S.(fields_main{find(is_main,1)});

        assert(numel(meg.rejection_info)==n_run, ...
            '%s，condition %d：rejection_info 的 run 数不正确。', ...
            sub_name,i_cond);

        % Collect original behavioral trial indices to delete
        bad_trials = [];

        for i_run = 1:n_run

            bad = meg.rejection_info{i_run}.bad_trials(:);

            assert(all(isfinite(bad) & bad==fix(bad) & ...
                       bad>=1 & bad<=trials_per_run), ...
                '%s，condition %d，run %d：bad_trials 索引无效。', ...
                sub_name,i_cond,i_run);

            bad_trials = [bad_trials; ...
                bad + trials_per_run*(i_run-1)];
        end
        bad_trials = unique(bad_trials);
        n_trial = n_run*trials_per_run;
        assert(all([numel(bhv.answer2),numel(bhv.conflict), ...
                    numel(bhv.ischange)]==n_trial), ...
            '%s，condition %d：原始行为 trial 数不是 %d。', ...
            sub_name,i_cond,n_trial);
        keep = true(n_trial,1);
        keep(bad_trials) = false;
        assert(sum(keep)==numel(meg.trial), ...
            '%s，condition %d：删除后行为 trial 数与 MEG 不一致。', ...
            sub_name,i_cond);
        % Apply the same retained indices to all fields
        is_right = bhv.answer2(:)==2;
        is_conflict = bhv.conflict(:);
        ischange = bhv.ischange(:);
        run_idx = repelem((1:n_run)',trials_per_run);
        out_bhv1=true(240,1);out_bhv1(is_conflict==0 & ischange==1)=false;
        data_vdx{i_sub,i_cond+1} = struct( ...
            'is_right',is_right(keep), ...
            'run_idx',run_idx(keep), ...
            'is_conflict',is_conflict(keep), ...
            'ischange',ischange(keep), ...
            'is_valid',out_bhv1);
    end
end
% Save only data_vdx
output_dir = fileparts(output_file);
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(output_file,'data_vdx');













%% data preprocess
function MEGData_Optimal_Pipeline(sub_ots, ots_names, trigger_infor, base_path, context_name, drug_name)
    % Based on FLUX Pipeline Reconstructed task-based MEG complete preprocessing workflow
    % Workflow: Continuous filtering -> Extract valid epochs -> Visual rejection -> Downsampling ICA -> Back-project to the original high-sampling-rate data -> Baseline correction/Saving

    % --- 1. Parameter and path settings ---
    cfg_preproc = [];
    cfg_preproc.hpfreq                = 0.5;
    cfg_preproc.eog_channels          = {'EOG001', 'EOG002'};
    cfg_preproc.ecg_channels          = {'ECG003'};
    cfg_preproc.ica_corr_threshold    = 0.4;
    cfg_preproc.ica_numcomponent      = 60;
    cfg_preproc.resamplefs_ica        = 200; % For accelerating ICA downsampling rate

    % Output path settings
    output_dir = fullfile(base_path, trigger_infor.phaseName,[context_name,'_',drug_name,'_', '306(',mat2str(trigger_infor.resamplefs),'hz)']);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % Get participant names(For matching behavioral data)
    sub_names = cell(length(sub_ots),1);
    for i = 1:length(sub_ots)
        if ~ischar(sub_ots{i}(1).expdata.subsName)
            sub_names{i,1} = regexprep(sub_ots{i}(1).expdata.subsName{1},'\d','');
        else
            sub_names{i,1} = regexprep(sub_ots{i}(1).expdata.subsName,'\d','');
        end
    end

    %% --- 2. Iterate over participants and Run ---
    for i_sub = 1:size(ots_names, 1)
        % From ots_names Get subject_id (Assume extraction by a file-path regular expression, Or use the naming convention)
        % Use simple splitting here, Adjust to the actual ots_names format
        input_filepath = ots_names{i_sub, 1};
        path_parts = strsplit(input_filepath, {'\', '/'});
        subject_id = path_parts{end-2};
        if isempty(subject_id), subject_id = sprintf('sub%02d', i_sub); end
        
        fprintf('\n=================================================\n');
        fprintf('===  开始处理被试: %s (Context: %s_%s)  ===\n', subject_id, context_name,drug_name);
        fprintf('=================================================\n');

        all_runs_data_for_subject = cell(size(ots_names, 2), 1);
        all_runs_rej_info = cell(size(ots_names, 2), 1);

        for i_run = 1:size(ots_names, 2)
            input_filepath = ots_names{i_sub, i_run};
            if isempty(input_filepath) || ~exist(input_filepath, 'file')
                continue;
            end
            
            fprintf('\n--- 处理 Run %d ---\n', i_run);

            %% Step 2.1: Read continuous data and apply basic filtering (Prevent edge effects)
            fprintf('读取数据与连续波滤波 (0.5Hz HPF & 50Hz Notch)...\n');
            cfg = [];
            cfg.dataset = input_filepath;
            cfg.channel = [{'MEG'}, cfg_preproc.eog_channels(:)', cfg_preproc.ecg_channels(:)'];
            data_continuous = ft_preprocessing(cfg);
            
            % Adaptive band-stop filtering (Use the provided function)
            cfg_bsf = adaptive_notch_filter(data_continuous, 50);
            cfg_filt = [];
            cfg_filt.hpfilter   = 'yes';
            cfg_filt.hpfreq     = cfg_preproc.hpfreq;
            cfg_filt.hpfilttype = 'fir';
            cfg_filt.hpfiltdir  = 'twopass';
            cfg_filt.hpfiltord  = 4;
            cfg_filt.bsfilter   = cfg_bsf.bsfilter;
            if isequal(cfg_bsf.bsfilter,'yes')        
                cfg_filt.bsfilttype = cfg_bsf.bsfilttype;
                cfg_filt.bsfiltord  = cfg_bsf.bsfiltord;
                cfg_filt.bsfiltdir  = cfg_bsf.bsfiltdir;
                cfg_filt.bsfreq     = cfg_bsf.bsfreq;
            end
            data_filtered = ft_preprocessing(cfg_filt, data_continuous);

            %% Step 2.2: Based on behavioral data and Trigger epoch the data (Epoching)
            fprintf('依据 Trigger 进行分段提取 (Epoching)...\n');
            cfg_trial = [];
            cfg_trial.dataset             = input_filepath;
            cfg_trial.trialfun            = 'ft_trialfun_general';
            cfg_trial.trialdef.eventtype  = 'STI101';
            cfg_trial.trialdef.eventvalue = trigger_infor.number;
            cfg_trial.trialdef.prestim    = trigger_infor.prestim;
            cfg_trial.trialdef.poststim   = trigger_infor.poststim;
            cfg_trials_tmp = ft_definetrial(cfg_trial);
            
            % Extract all events to match the getTriggerIndex Logic
            event_values = [cfg_trials_tmp.event(strcmp({cfg_trials_tmp.event.type}, 'STI101')).value]';
            
            % Get behavioral-data indices
            bh_idx = find(strcmp(sub_names, regexprep(subject_id,'-\d','')));
            if sum(bh_idx)==0, error('未找到对应被试行为数据，请检查索引'); end
            bh_tmp = sub_ots{bh_idx}(i_run).expdata;
            
            % Get valid trial indices
            [~, is_v4] = getTriggerIndex(subject_id, event_values, bh_tmp, i_run, trigger_infor.number);
            trigger4_samples = [cfg_trials_tmp.event([cfg_trials_tmp.event.value] == trigger_infor.number & strcmp({cfg_trials_tmp.event.type}, 'STI101')).sample]';
            
            % Construct custom trl matrix
            Fs = data_filtered.fsample;
            cfg_trl = [];
            cfg_trl.trl = [];
            
            % Process shift_time (Compatibility handling)
            if isfield(trigger_infor, 'shift_time')
                shift_time = eval(trigger_infor.shift_time);
            else
                shift_time = zeros(length(trigger4_samples),1);
            end
            if isscalar(shift_time), shift_time = repmat(shift_time, length(trigger4_samples), 1); end
            shift_sample = round(shift_time * Fs);
            
            prestim_samples = round(trigger_infor.prestim * Fs);
            poststim_samples = round(trigger_infor.poststim * Fs);
            
            n = 1;
            for i = 1:length(trigger4_samples)
                if is_v4(i) == 1
                    start_sample = max(1, trigger4_samples(i) - (prestim_samples + shift_sample(n)));
                    end_sample   = min(data_filtered.hdr.nSamples, trigger4_samples(i) + (poststim_samples - shift_sample(n)));
                    offset       = -prestim_samples;
                    
                    % If required downstream, store additional information in trialinfo
                    % cfg_trl.trl(end+1, :) = [start_sample, end_sample, offset, trigger_infor.number];
                    % n is the current run valid within trial behavioral trial_id
                    cfg_trl.trl(end+1, :) = [ ...
                        start_sample, ...
                        end_sample, ...
                        offset, ...
                        trigger_infor.number, ...
                        n, ...        % trial_id within run
                        i_run];       % run_id
                    n = n + 1;
                end
            end
            
            % Based on the generated trl epoch the filtered continuous data
            data_epoched = ft_redefinetrial(cfg_trl, data_filtered);

            % % ==========================================
            % % Step 2.3: Manually reject severe artifacts by visual inspection (Protected data-isolation version)
            % % ==========================================
            % fprintf('\n Isolating data, Protect electrooculography (EOG)/electrocardiography (ECG) from accidental deletion...\n');
            % 
            % % --- 1. Split data: Extract only MEG and Non- MEG ---
            % cfg_sel = [];
            % cfg_sel.channel = 'MEG';
            % data_meg_only = ft_selectdata(cfg_sel, data_epoched);
            % 
            % cfg_sel.channel = setdiff(data_epoched.label, data_meg_only.label);
            % data_nonmeg_only = ft_selectdata(cfg_sel, data_epoched); % retain separately
            % 
            % % --- 2. Open the visual rejection interface only for MEG data ---
            % fprintf('Open visual rejection interface: Reject bad channels or bad epochs safely Trial (Non- MEG channels are protected in the background)...\n');
            % cfg_rej = [];
            % cfg_rej.method = 'summary';
            % cfg_rej.keepchannel = 'no'; % Reject bad channels
            % % The interface will not contain EOG and ECG!
            % data_meg_clean = ft_rejectvisual(cfg_rej, data_meg_only);
            % 
            % % --- 3. Extract rejection statistics, and update non- MEG data ---
            % % Find retained trial identifier
            % original_trials = 1:length(data_meg_only.trial);
            % kept_trials = find(ismember(data_meg_only.sampleinfo(:,1), data_meg_clean.sampleinfo(:,1)));
            % bad_trials = setdiff(original_trials, kept_trials);
            % bad_channels = setdiff(data_meg_only.label, data_meg_clean.label);
            % 
            % fprintf('\n--- Rejection statistics ---\n');
            % fprintf('Rejected bad Trials: %s\n', mat2str(bad_trials));
            % fprintf('Rejected bad Channels: %s\n', strjoin(bad_channels, ', '));
            % 
            % % ==========================================
            % % [Critical fix]: Store rejection information in an external variable
            % % ==========================================
            % rej_info = struct();
            % rej_info.bad_trials = bad_trials;
            % rej_info.bad_channels = bad_channels;
            % rej_info.num_original_trials = length(original_trials);
            % rej_info.num_kept_trials = length(kept_trials);
            % 
            % % Store in the externally defined cell array
            % all_runs_rej_info{i_run} = rej_info;
            % 
            % % Critical: Make the separately retained non- MEG data drop the same bad trial, Keep matrix dimensions exactly aligned!
            % cfg_sync = [];
            % cfg_sync.trials = kept_trials;
            % data_nonmeg_clean = ft_selectdata(cfg_sync, data_nonmeg_only);
            % % --- 4. Repair rejected MEG bad channels ---
            % if ~isempty(bad_channels)
            %     fprintf('Detected rejected MEG channels, Automatically interpolate to repair...\n');
            %     cfg_neighb = [];
            %     cfg_neighb.method    = 'triangulation';
            %     cfg_neighb.layout    = 'neuromag306all.lay';
            %     neighbours = ft_prepare_neighbours(cfg_neighb, data_meg_only);
            % 
            %     cfg_repair = [];
            %     cfg_repair.badchannel = bad_channels;
            %     cfg_repair.neighbours = neighbours;
            %     cfg_repair.method     = 'spline';
            %     data_meg_clean = ft_channelrepair(cfg_repair, data_meg_clean);
            % end
            % 
            % % --- 5. Lossless concatenation ---
            % cfg_append = [];
            % data_epoched_clean = ft_appenddata(cfg_append, data_meg_clean, data_nonmeg_clean);
            % 
            % % Restore original channel order, Prevent downstream errors caused by channel reordering
            % cfg_sort = [];
            % cfg_sort.channel = data_epoched.label;
            % data_epoched_clean = ft_selectdata(cfg_sort, data_epoched_clean);
            % 
            % fprintf('Data cleaning and merging completed, All channels proceed to the next step!\n');

            % ==========================================
            % Step 2.3 (Automated version): Variance-based Z-score Automatically reject severe artifacts (Protected data-isolation version)
            % ==========================================
            fprintf('\n正在执行自动化坏通道与坏 Trial 检测 (已在后台隔离保护非MEG通道)...\n');

            % Set the anomaly-detection Z-score threshold (Suggested value: 3.5 Up to 4)
            % If too much is rejected, increase to 4.5; If artifacts remain, decrease to 3
            z_threshold_chan  = 3.5;
            z_threshold_trial = 3.5;

            % --- 0. Isolate data first: Extract only MEG and Non- MEG ---
            cfg_sel = [];
            cfg_sel.channel = 'MEG';
            data_meg_only = ft_selectdata(cfg_sel, data_epoched);

            cfg_sel.channel = setdiff(data_epoched.label, data_meg_only.label);
            data_nonmeg_only = ft_selectdata(cfg_sel, data_epoched); % retain separately

            % --- 1. Get GRAD and MAG indices (Must compute separately in pure MEG data) ---
            grad_idx = match_str(data_meg_only.label, ft_channelselection('MEGGRAD', data_meg_only.label));
            mag_idx  = match_str(data_meg_only.label, ft_channelselection('MEGMAG', data_meg_only.label));

            num_trials = length(data_meg_only.trial);

            % Compute each Trial in GRAD and MAG channel's variance matrix
            var_grad = zeros(length(grad_idx), num_trials);
            var_mag  = zeros(length(mag_idx), num_trials);

            for trl = 1:num_trials
                trial_data = data_meg_only.trial{trl};
                % Across the time dimension (Dimension 2 index) compute variance
                var_grad(:, trl) = var(trial_data(grad_idx, :), 0, 2);
                var_mag(:, trl)  = var(trial_data(mag_idx, :), 0, 2);
            end

            % --- A. Automatically detect bad channels (Across Trial excessive mean variance) ---
            chan_mean_var_grad = mean(var_grad, 2);
            chan_mean_var_mag  = mean(var_mag, 2);

            bad_channels = {};
            if ~isempty(grad_idx)
                z_grad_chan = zscore(chan_mean_var_grad);
                bad_channels = [bad_channels; data_meg_only.label(grad_idx(z_grad_chan > z_threshold_chan))];
            end
            if ~isempty(mag_idx)
                z_mag_chan = zscore(chan_mean_var_mag);
                bad_channels = [bad_channels; data_meg_only.label(mag_idx(z_mag_chan > z_threshold_chan))];
            end

            % --- B. Automatically detect bad Trial (transient extremes) ---
            % Exclude effects of already identified bad channels (When calculating Trial exclude them from the mean)
            valid_grad_local_idx = find(~ismember(data_meg_only.label(grad_idx), bad_channels));
            valid_mag_local_idx  = find(~ismember(data_meg_only.label(mag_idx), bad_channels));

            bad_trials = [];
            if ~isempty(valid_grad_local_idx)
                trial_mean_var_grad = mean(var_grad(valid_grad_local_idx, :), 1);
                z_trial_grad = zscore(trial_mean_var_grad);
                bad_trials = [bad_trials, find(z_trial_grad > z_threshold_trial)];
            end
            if ~isempty(valid_mag_local_idx)
                trial_mean_var_mag = mean(var_mag(valid_mag_local_idx, :), 1);
                z_trial_mag = zscore(trial_mean_var_mag);
                bad_trials = [bad_trials, find(z_trial_mag > z_threshold_trial)];
            end
            bad_trials = unique(bad_trials); % Remove duplicates

            % --- C. Apply rejection ---
            kept_trials   = setdiff(1:num_trials, bad_trials);
            kept_channels = setdiff(data_meg_only.label, bad_channels);

            cfg_rej_auto = [];
            cfg_rej_auto.channel = kept_channels;
            cfg_rej_auto.trials  = kept_trials;
            data_meg_clean = ft_selectdata(cfg_rej_auto, data_meg_only);

            % [Critical synchronization]: Make the separately retained non- MEG data also reject the same bad trial
            cfg_sync = [];
            cfg_sync.trials = kept_trials;
            data_nonmeg_clean = ft_selectdata(cfg_sync, data_nonmeg_only);

            % --- D. Record rejection log ---
            rej_info = struct();
            rej_info.bad_trials = bad_trials;
            rej_info.bad_channels = bad_channels;
            rej_info.num_original_trials = num_trials;
            rej_info.num_kept_trials = length(kept_trials);
            all_runs_rej_info{i_run} = rej_info;

            fprintf('\n--- 自动化剔除统计 (Run %d) ---\n', i_run);
            fprintf('检测到坏 Trials: %d 个 (索引: %s)\n', length(bad_trials), mat2str(bad_trials));
            fprintf('检测到坏 Channels: %d 个 (%s)\n', length(bad_channels), strjoin(bad_channels, ', '));

            % --- E. Interpolate to repair rejected MEG bad channels ---
            if ~isempty(bad_channels)
                fprintf('\n正在使用周围通道对剔除的 MEG 坏通道进行自动插值修复...\n');
                cfg_neighb = [];
                cfg_neighb.method    = 'triangulation';
                cfg_neighb.layout    = 'neuromag306all.lay';
                neighbours = ft_prepare_neighbours(cfg_neighb, data_meg_only); % Based on the original complete MEG compute neighbors

                cfg_repair = [];
                cfg_repair.badchannel = bad_channels;
                cfg_repair.neighbours = neighbours;
                cfg_repair.method     = 'spline';
                data_meg_clean = ft_channelrepair(cfg_repair, data_meg_clean);
            end

            % --- F. Lossless merging ---
            cfg_append = [];
            data_epoched_clean = ft_appenddata(cfg_append, data_meg_clean, data_nonmeg_clean);

            % Restore original channel order, Prevent downstream ICA operation errors
            cfg_sort = [];
            cfg_sort.channel = data_epoched.label;
            data_epoched_clean = ft_selectdata(cfg_sort, data_epoched_clean);

            fprintf('自动化清洗与合并完成，所有通道已安全进入下一步！\n');

            %% Step 2.4: Downsample and run ICA
            fprintf('降采样到 %d Hz 以加速 ICA...\n', cfg_preproc.resamplefs_ica);
            cfg_ds = [];
            cfg_ds.resamplefs = cfg_preproc.resamplefs_ica;
            data_epoched_ds   = ft_resampledata(cfg_ds, data_epoched_clean);
            
            fprintf('正在运行 ICA 分析...\n');
            cfg_ica = [];
            cfg_ica.method       = 'runica';
            cfg_ica.channel      = 'MEG';
            cfg_ica.numcomponent = cfg_preproc.ica_numcomponent;
            comp = ft_componentanalysis(cfg_ica, data_epoched_ds);

            % Automatically identify correlations for EOG/ECG bad components
            bad_components = [];
            artifact_channels = [cfg_preproc.eog_channels, cfg_preproc.ecg_channels];

            % Extract and concatenate all ICA components' trial data (Place outside the loop, Execute once only)
            % comp.trial is a cell array, cell2mat After concatenation becomes: [component count x total samples]
            comp_data_all = cell2mat(comp.trial);

            for i_art = 1:length(artifact_channels)
                % Extract the current artifact channel's data
                art_ts = ft_selectdata(struct('channel', artifact_channels{i_art}), data_epoched_ds);

                % Extract and concatenate the current artifact channel's trial data
                % art_ts.trial After concatenation becomes: [1 x total samples]
                art_data_all = cell2mat(art_ts.trial);
                art_cat = art_data_all(1, :);

                for i_comp = 1:numel(comp.label)
                    % Index continuous data for the current component using an intermediate variable
                    comp_cat = comp_data_all(i_comp, :);

                    % Ensure equal lengths(Normally lengths match after concatenation)
                    len = min(length(art_cat), length(comp_cat));

                    % Compute Pearson correlation coefficients
                    R = corrcoef(art_cat(1:len), comp_cat(1:len));

                    % Extract off-diagonal correlations and compare with the threshold
                    if abs(R(1,2)) > cfg_preproc.ica_corr_threshold
                        bad_components = [bad_components, i_comp];
                    end
                end
            end

            % Remove duplicates, A component may correlate strongly with multiple artifact channels
            bad_components = unique(bad_components);
            fprintf('\n==================================================\n');
            fprintf('算法自动检测到的可疑坏成分 (EOG/ECG) 为: %s\n', mat2str(bad_components));
            fprintf('==================================================\n');
            
            % Retain the existing plotting logic here (ft_topoplotIC etc.) For visual confirmation...
            fprintf('正在打开 databrowser 进行视觉确认...\n');
            fprintf('请重点检查上述算法推荐的成分。\n');
            fprintf('【操作提示】：检查完毕后请关闭 databrowser 窗口，代码将继续运行。\n');

            cfg_browse = [];
            cfg_browse.layout   = 'neuromag306all.lay'; % Or use neuromag306mag.lay / neuromag306planar.lay as needed
            cfg_browse.viewmode = 'component';
            ft_databrowser(cfg_browse, comp);

            % % 4. Allow command-line override of automatic detection(Strongly recommended!)
            % user_input = input('Confirm indices of components to reject (Press Enter to accept the automatic result, Or enter [1 2 5] to override): ', 's');
            % if ~isempty(user_input)
            %     % If a new array is entered, override automatic detection
            %     bad_components = str2num(user_input);
            % end
            fprintf('最终决定剔除的成分为: %s\n', mat2str(bad_components));
            %% Step 2.5: Restrict ICA Back-project results to high-sampling-rate data (Core optimization)
            fprintf('将 ICA 坏成分从高采样率分段数据中剔除...\n');
            cfg_rej_comp = [];
            cfg_rej_comp.component = bad_components;
            % Key point: Input is data_epoched_clean (not downsampled), Output is artifact-cleaned high-quality data
            data_clean = ft_rejectcomponent(cfg_rej_comp, comp, data_epoched_clean);

            %% Step 2.6: Baseline correction and final downsampling
            % Apply baseline correction
            if isfield(trigger_infor, 'do_demean') && trigger_infor.do_demean == false
                fprintf('Skipping epoch-wise demean because trigger_infor.do_demean = false.\n');
            else
                if isfield(trigger_infor, 'baselinewindow') && ~isempty(trigger_infor.baselinewindow)
                    cfg_base = [];
                    cfg_base.demean = 'yes';
                    cfg_base.baselinewindow = trigger_infor.baselinewindow;
                    data_clean = ft_preprocessing(cfg_base, data_clean);
                else
                    cfg_base = [];
                    cfg_base.demean = 'yes';
                    data_clean = ft_preprocessing(cfg_base, data_clean);
                end
            end
            
            % Final downsampling to the target analysis rate (For example 250Hz)
            cfg_final_ds = [];
            cfg_final_ds.resamplefs = trigger_infor.resamplefs;
            data_final = ft_resampledata(cfg_final_ds, data_clean);
            
            all_runs_data_for_subject{i_run, 1} = data_final;
            close all;
        end
        
        %% --- 3. Merge Runs and save ---
        % Remove empty run slots (Contains data and rejection information)
        valid_runs = ~cellfun('isempty', all_runs_data_for_subject);
        all_runs_data_for_subject = all_runs_data_for_subject(valid_runs);
        all_runs_rej_info = all_runs_rej_info(valid_runs); % [Additional]Remove empty rejection information concurrently

        if ~isempty(all_runs_data_for_subject)
            fprintf('\n--- 正在合并被试 %s 的所有有效 runs ---\n', subject_id);
            cfg_app = [];
            cfg_app.keepsampleinfo = 'no';
            data_combined_pre = ft_appenddata(cfg_app, all_runs_data_for_subject{:});
            data_combined_pre.grad = all_runs_data_for_subject{1}.grad;

            % [Core modification]: Attach all run rejection information as a new field in the aggregate data
            data_combined_pre.rejection_info = all_runs_rej_info;
 
            output_filename = sprintf('%s_%s_%s_%s_306.mat', ...
                subject_id, trigger_infor.phaseName, context_name, drug_name);
            output_filepath = fullfile(output_dir, output_filename);

            fprintf('正在保存包含完整剔除履历的合并文件到: %s\n', output_filepath);
            save(output_filepath, 'data_combined_pre', '-v7.3');
        end
    end
end
%%
function cfg = adaptive_notch_filter(data, base_frequency)
    % Adapt band-stop filtering to the observed line noise
    
    % 1. Analyze spectral features
    cfg_spectrum = [];
    cfg_spectrum.method = 'mtmfft';
    cfg_spectrum.taper = 'hanning';
    cfg_spectrum.foi = 1:300;
    cfg_spectrum.channel = 'MEG';
    spectrum = ft_freqanalysis(cfg_spectrum, data);
    
    % 2. Detect line-noise peaks
    noise_peaks = detect_powerline_peaks(spectrum, base_frequency);
    
    % 3. Choose filtering parameters by noise strength
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfilttype = 'fir';
    cfg.bsfiltdir = 'twopass';
    
    if isempty(noise_peaks)
        % No obvious line noise, Do not use band-stop filtering
        cfg.bsfilter = 'no';
        fprintf('未检测到明显线噪声，跳过带阻滤波\n');
    else
        % Set filtering parameters by noise strength
        cfg.bsfreq = [];
        for i = 1:length(noise_peaks)
            freq = noise_peaks(i).frequency;
            strength = noise_peaks(i).strength;
            
            % Choose stopband width by noise strength
            if strength > 5  % Strong noise
                bandwidth = 2;  % Wide stopband
                filt_order = 6;
            elseif strength > 2  % Moderate noise
                bandwidth = 1.5;
                filt_order = 4;
            else  % Weak noise
                bandwidth = 1;
                filt_order = 2;
            end
            
            cfg.bsfreq = [cfg.bsfreq; freq-bandwidth, freq+bandwidth];
            
            if i == 1
                cfg.bsfiltord = filt_order;
            end
        end
        
        fprintf('检测到线噪声峰值: ');
        for i = 1:length(noise_peaks)
            fprintf('%.1fHz(强度=%.2f) ', noise_peaks(i).frequency, noise_peaks(i).strength);
        end
        fprintf('\n');
    end
end

function peaks = detect_powerline_peaks(spectrum, base_freq)
    % Detect peaks at line frequency and its harmonics
    peaks = [];
    
    harmonic_freqs = base_freq:base_freq:base_freq*5;  % Fundamental and 2 order, 3 harmonic
    
    for i = 1:length(harmonic_freqs)
        freq = harmonic_freqs(i);
        % Find peaks near the target frequency
        freq_range = [freq-2, freq+2];
        idx = (spectrum.freq >= freq_range(1)) & (spectrum.freq <= freq_range(2));
        
        if any(idx)
            power_in_band = mean(mean(abs(spectrum.powspctrm(:, idx)), 1));
            power_around = mean(mean(abs(spectrum.powspctrm(:, ...
                (spectrum.freq >= freq-10) & (spectrum.freq <= freq-5) | ...
                (spectrum.freq >= freq+5) & (spectrum.freq <= freq+10))), 1));
            
            % Compute signal-to-noise ratio
            snr_ratio = power_in_band / power_around;
            
            if snr_ratio > 1.2  % If the peak is prominent
                peaks(end+1).frequency = freq;
                peaks(end).strength = snr_ratio;
                peaks(end).bandwidth = 2;  % Initial bandwidth estimate
            end
        end
    end
end

function [ots_names,otn_names]=getPath(base_input_path)
% --- path name ---
sub_name = dir(base_input_path);
% Recommended to use tsss file
ots_names = cell(length(sub_name)-2,3);
otn_names = cell(length(sub_name)-2,3);
for i=1:length(sub_name)-2
    sub_path = [base_input_path,'\',sub_name(i+2).name];
    sub_name1 = dir(sub_path);
    sub_path0 = [sub_path,'\',sub_name1(3).name,'\'];    
    for j=1:3
        % ots_names{i,j}= [sub_path0,'s',mat2str(j),'_quat_trans_tsss.fif'];
        % otn_names{i,j}= [sub_path0,'n',mat2str(j),'_quat_trans_tsss.fif'];
        ots_names{i,j}= [sub_path0,'s',mat2str(j),'_tsss_mc.fif'];
        otn_names{i,j}= [sub_path0,'n',mat2str(j),'_tsss_mc.fif'];        
    end   
end
end


function [is_v8,is_v4]=getTriggerIndex(subject_id,event_values,bh_tmp,i_run,t_number)
        % get behaviour data
        % bh_idx=find(strcmp(sub_names,regexprep(subject_id,'-\d','')));
        % if sum(bh_idx)==0
        % error('Matching participant not found, Check indices');
        % end
        % bh_tmp=sub_ots{bh_idx}(i_run).expdata;
        if length(bh_tmp.retry)~=length(bh_tmp.oretry)
            bh_tmp.retry(length(bh_tmp.retry)+1:length(bh_tmp.oretry))=nan;
            bh_tmp.retry(bh_tmp.retry==0)=nan;
        end
        is_retry=bh_tmp.retry;is_retry(isnan(bh_tmp.oretry))=[];

        if ~isequal(unique(event_values),[1 2 3 4 6 8]')
            fprintf('-----------------------------------------\n');
            tabulate(event_values);
            ev_unique=unique(event_values);idx=ismember(unique(event_values),[1 2 3 4 6 8]');
            % input('>>> Check the report above.By[Enter]key to continue, Or press[Ctrl+C]to terminate the script.');
            if sum(event_values==8)~=length(is_retry)
                event_values(event_values==ev_unique(idx==0))=8;
            end
        end
        % combine behaviore data to get trial isvalid idx
        t8_pos=find(event_values==8);
        t8_idx=ones(sum(~isnan(bh_tmp.oretry)),1);
        t4_idx=ones(sum(~isnan(bh_tmp.oretry)),1);
        if t8_pos(1)==1
        t8_pos0=[t8_pos;length(event_values)];
        else
        t8_pos0=[1;t8_pos;length(event_values)];
        end
        t6_pos=find(event_values==6);
        mn=0;
        mn6=1;
        for n=1:sum(~isnan(bh_tmp.oretry))
            if isnan(is_retry(n))
                t8_idx(n)=sum(event_values(mn+1:t6_pos(mn6))==8);
                t4_idx(n)=sum(event_values(mn+1:t6_pos(mn6))==t_number);
                mn=t6_pos(mn6);
                mn6=mn6+1;
            elseif is_retry(n)==1
                if event_values(mn+1)==8
                mn_x=1;mn=mn+1;
                while event_values(mn+mn_x)~=8
                    mn_x=mn_x+1;
                end
                t8_idx(n)=sum(event_values(mn:mn+mn_x-1)==8);
                t4_idx(n)=sum(event_values(mn:mn+mn_x-1)==t_number);                
                if sum(event_values(mn:mn+mn_x-1)==6)==1
                    mn6=mn6+1;
                end
                else
                    error(sprintf('something in sub%s run%d',subject_id,i_run));
                end
                mn=mn+mn_x-1;
            end
        end
        if sum(t8_idx)~=sum(event_values==8) || sum(t4_idx)~=sum(event_values==t_number)
            error('trigger number is wrong, 请检查trial valid idx');
        end
        if length(t8_idx)~=sum(~isnan(bh_tmp.oretry)) || length(t4_idx)~=sum(~isnan(bh_tmp.oretry))
            error('trigger trial index is wrong, 请检查trial valid idx');
        end
        is_v4=getInvalidIndix(bh_tmp,t4_idx);
        is_v8=getInvalidIndix(bh_tmp,t8_idx);

end

%%
function is_v4=getInvalidIndix(datas,t4_idx)
is_retry=datas.retry(~isnan(datas.oretry));
is_valid=isnan(is_retry);
is_retry4=is_valid(t4_idx~=0);
is_t4=t4_idx(t4_idx~=0);
n=1;is_v4=[];
for j=1:length(is_t4)
    if is_t4(j)==1
        is_v4(n,1)=is_retry4(j);
        n=n+1;
    elseif is_t4(j)>1
        for m=1:is_t4(j)
            if m==1
                is_v4(n,1)=is_retry4(j);
                n=n+1;
            else
                is_v4(n,1)=0;
                n=n+1;
            end
        end
    end
end
end


