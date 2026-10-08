% KH06 fMRS batch processing. Set input paths in kh06_config.m.

rehash;
clc; clear; close all;

cfg = kh06_config();
base_dir = cfg.data_dir;
beh_path = cfg.behavior_file;

if exist(beh_path, 'file')
    fprintf('>>> [Setup] Loading behavioral data: %s ...\n', beh_path);
    load(beh_path, 'sub_data');
else
    error('Behavioral data file not found: %s', beh_path);
end

[sub_folders, valid_sub_indices, sub_ids] = find_subject_folders(base_dir);
if isempty(valid_sub_indices)
    error('No subject folder containing S followed by digits was found under %s', base_dir);
end

fprintf('>>> Found %d subjects. Starting batch processing.\n', length(valid_sub_indices));
pause(1);

for s_idx = 1:length(valid_sub_indices)
    this_folder_struct = sub_folders(valid_sub_indices(s_idx));
    subj_folder_name = this_folder_struct.name;
    subj_full_path = fullfile(base_dir, subj_folder_name);
    target_Subject = sub_ids(s_idx);

    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('>>> Processing subject: %s (ID: %d) [%d/%d]\n', ...
        subj_folder_name, target_Subject, s_idx, length(valid_sub_indices));
    fprintf('%s\n', repmat('=', 1, 60));

    for target_Run = 1:4
        fprintf('\n    --------------------------------------------------\n');
        fprintf('    >>> Checking Run %d ...\n', target_Run);
        try
            process_single_run(target_Subject, target_Run, subj_full_path, ...
                subj_folder_name, base_dir, sub_data);
        catch ME
            fprintf('\n    !!! [Error] Run %d failed !!!\n', target_Run);
            fprintf('    Error message: %s\n', ME.message);
            fprintf('    >>> Skipping this run and continuing...\n');
        end
        close all;
    end
end

fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('>>> Batch processing completed.\n');
fprintf('%s\n', repmat('=', 1, 60));

function [sub_folders, valid_sub_indices, sub_ids] = find_subject_folders(base_dir)
    % Keep the final S-number match in each eligible folder name.
    all_contents = dir(base_dir);
    dir_flags = [all_contents.isdir];
    sub_folders = all_contents(dir_flags);
    valid_sub_indices = [];
    sub_ids = [];

    fprintf('\n>>> [Setup] Scanning subject folders...\n');
    for k = 1:length(sub_folders)
        fname = sub_folders(k).name;
        tokens = regexp(fname, 'S(\d+)', 'tokens');
        if ~isempty(tokens)
            sid_str = tokens{end}{1};
            sid = str2double(sid_str);
            if ~contains(fname, 'Results', 'IgnoreCase', true)
                valid_sub_indices = [valid_sub_indices, k];
                sub_ids = [sub_ids, sid];
                fprintf('    Found subject folder: %-25s (ID: %d)\n', fname, sid);
            end
        end
    end
end

function process_single_run(target_Subject, target_Run, subj_full_path, subj_folder_name, base_dir, sub_data)
    subj_str_filename = sprintf('S%02d', target_Subject);
    [mrs_filename, mrs_file] = find_metabolite_file(subj_full_path, target_Run);

    if contains(mrs_filename, 'nonsocial', 'IgnoreCase', true)
        cond_name = 'nonsocial';
    elseif contains(mrs_filename, 'social', 'IgnoreCase', true)
        cond_name = 'social';
    else
        cond_name = 'unknown';
    end

    target_subj_str = sprintf('S0%d', target_Subject);
    row_idx = find(strcmp(sub_data(:, 2), target_subj_str));
    if isempty(row_idx)
        fallback_str = sprintf('S%02d', target_Subject);
        row_idx = find(strcmp(sub_data(:, 2), fallback_str));
        if ~isempty(row_idx)
            fprintf('    [Note] "%s" was not found; using "%s" instead.\n', ...
                target_subj_str, fallback_str);
            target_subj_str = fallback_str;
        end
    end
    if isempty(row_idx)
        error('Subject "%s" not found in sub_data (folder ID: %d)', target_subj_str, target_Subject);
    end

    fprintf('    [Match] Behavioral data ID: %s (row %d)\n', target_subj_str, row_idx);
    run_data = sub_data{row_idx, 1}(target_Run);
    if isfield(run_data, 'context')
        beh_context = run_data.context;
        if ~strcmpi(cond_name, beh_context)
            cond_name = beh_context;
        end
    end

    current_exp_order = 1;
    try
        if isfield(sub_data{row_idx, 1}, 'exp_order')
            current_exp_order = sub_data{row_idx, 1}(1).exp_order;
        end
    catch
    end

    split_name = strsplit(subj_folder_name, '_');
    if length(split_name) >= 3
        res_folder_name = sprintf('%s_%s_Results', split_name{1}, split_name{2});
    else
        res_folder_name = [subj_folder_name '_Results'];
    end
    output_dir = fullfile(base_dir, res_folder_name);
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    output_mat_name = fullfile(output_dir, sprintf('%s_Run%d_%s_exp_order%d_Results.mat', ...
        subj_str_filename, target_Run, cond_name, current_exp_order));
    if exist(output_mat_name, 'file')
        [~, fname_short, ~] = fileparts(output_mat_name);
        fprintf('    >>> [Skip] Result file already exists: %s\n', fname_short);
        return;
    end

    fprintf('    [Process] Starting analysis: %s\n', mrs_filename);
    water_file = find_water_file(subj_full_path);
    nii_file = find_t1_file(subj_full_path);
    [ranges_conf_p1, ranges_conf_p2, ranges_cons_p1, ranges_cons_p2] = ...
        build_trial_ranges(run_data);

    MRS_conf_p1 = []; MRS_conf_p2 = []; MRS_cons_p1 = []; MRS_cons_p2 = [];
    if ~isempty(ranges_conf_p1)
        MRS_conf_p1 = run_gannet_internal(mrs_file, water_file, nii_file, ranges_conf_p1, 'Conflict_Phase1');
    end
    if ~isempty(ranges_conf_p2)
        MRS_conf_p2 = run_gannet_internal(mrs_file, water_file, nii_file, ranges_conf_p2, 'Conflict_Phase2');
    end
    if ~isempty(ranges_cons_p1)
        MRS_cons_p1 = run_gannet_internal(mrs_file, water_file, nii_file, ranges_cons_p1, 'Consistent_Phase1');
    end
    if ~isempty(ranges_cons_p2)
        MRS_cons_p2 = run_gannet_internal(mrs_file, water_file, nii_file, ranges_cons_p2, 'Consistent_Phase2');
    end

    vars_to_save = {'ranges_conf_p1', 'ranges_conf_p2', 'ranges_cons_p1', 'ranges_cons_p2', ...
        'cond_name', 'subj_folder_name', 'mrs_filename', 'nii_file'};
    if ~isempty(MRS_conf_p1) && isfield(MRS_conf_p1, 'out')
        MRS_conf_p1 = MRS_conf_p1.out;
        vars_to_save{end+1} = 'MRS_conf_p1';
    end
    if ~isempty(MRS_conf_p2) && isfield(MRS_conf_p2, 'out')
        MRS_conf_p2 = MRS_conf_p2.out;
        vars_to_save{end+1} = 'MRS_conf_p2';
    end
    if ~isempty(MRS_cons_p1) && isfield(MRS_cons_p1, 'out')
        MRS_cons_p1 = MRS_cons_p1.out;
        vars_to_save{end+1} = 'MRS_cons_p1';
    end
    if ~isempty(MRS_cons_p2) && isfield(MRS_cons_p2, 'out')
        MRS_cons_p2 = MRS_cons_p2.out;
        vars_to_save{end+1} = 'MRS_cons_p2';
    end

    save(output_mat_name, vars_to_save{:});
    fprintf('    >>> Saved: %s\n', output_mat_name);
end

function [mrs_filename, mrs_file] = find_metabolite_file(subj_full_path, target_Run)
    % Preserve the original first-match selection after reference filtering.
    mrs_pattern = sprintf('*run%d*.dat', target_Run);
    mrs_files = dir(fullfile(subj_full_path, mrs_pattern));
    valid_mrs_idx = [];
    for k = 1:length(mrs_files)
        fname = mrs_files(k).name;
        is_water_ref = contains(fname, 'Water', 'IgnoreCase', true) && ...
            ~contains(fname, 'NOWater', 'IgnoreCase', true);
        is_other_ref = contains(fname, 'ref', 'IgnoreCase', true);
        if ~is_water_ref && ~is_other_ref
            valid_mrs_idx = [valid_mrs_idx, k];
        end
    end
    if isempty(valid_mrs_idx)
        error('Metabolite .dat file not found for run%d', target_Run);
    end
    mrs_filename = mrs_files(valid_mrs_idx(1)).name;
    mrs_file = fullfile(subj_full_path, mrs_filename);
end

function water_file = find_water_file(subj_full_path)
    water_files = dir(fullfile(subj_full_path, '*Water*.dat'));
    real_water_idx = [];
    for k = 1:length(water_files)
        if ~contains(water_files(k).name, 'NOWater', 'IgnoreCase', true)
            real_water_idx = [real_water_idx, k];
        end
    end
    if isempty(real_water_idx), error('Water reference file not found'); end
    water_file = fullfile(subj_full_path, water_files(real_water_idx(1)).name);
end

function nii_file = find_t1_file(subj_full_path)
    t1_pattern = '*t1_mprage_sag_p2*.nii';
    all_t1_candidates = dir(fullfile(subj_full_path, '**', t1_pattern));
    raw_t1_list = [];
    for k = 1:length(all_t1_candidates)
        fname = all_t1_candidates(k).name;
        is_seg_file = (length(fname) > 1 && fname(1) == 'c' && isstrprop(fname(2), 'digit'));
        is_norm_file = strncmpi(fname, 'y_', 2);
        if ~is_seg_file && ~is_norm_file
            raw_t1_list = [raw_t1_list; all_t1_candidates(k)];
        end
    end
    if isempty(raw_t1_list), error('Original T1 file not found'); end

    series_nums = zeros(length(raw_t1_list), 1);
    for k = 1:length(raw_t1_list)
        [~, fn, ~] = fileparts(raw_t1_list(k).name);
        parts = strsplit(fn, '_');
        val = str2double(parts{end});
        if isnan(val), val = 9999; end
        series_nums(k) = val;
    end
    [~, sort_idx] = sort(series_nums, 'ascend');
    t1_target = raw_t1_list(sort_idx(1));
    nii_file = fullfile(t1_target.folder, t1_target.name);
end

function [ranges_conf_p1, ranges_conf_p2, ranges_cons_p1, ranges_cons_p2] = build_trial_ranges(run_data)
    % Preserve the original time-to-frame conversion and phase windows.
    starts = run_data.trial_start;
    conflict = run_data.isconflict;
    valid = run_data.isvalid;
    TR = 1.6;
    TRS_PER_TRIAL_DEFAULT = 4;
    TOTAL_POINTS = 400;
    ranges_conf_p1 = []; ranges_conf_p2 = []; ranges_cons_p1 = []; ranges_cons_p2 = [];

    for i = 1:length(starts)
        if valid(i) == 1
            start_idx = round(starts(i) / TR);
            if i < length(starts)
                duration_sec = starts(i+1) - starts(i);
                current_trs_count = round(duration_sec / TR);
            else
                current_trs_count = TRS_PER_TRIAL_DEFAULT;
            end
            if current_trs_count > 4 || current_trs_count < 2, continue; end
            range_p1 = [start_idx, start_idx + 1];
            if range_p1(2) > TOTAL_POINTS, continue; end
            if conflict(i) == 1
                ranges_conf_p1 = [ranges_conf_p1; range_p1];
            else
                ranges_cons_p1 = [ranges_cons_p1; range_p1];
            end
            if current_trs_count >= 4
                range_p2 = [start_idx + 2, start_idx + 3];
                if range_p2(2) <= TOTAL_POINTS
                    if conflict(i) == 1
                        ranges_conf_p2 = [ranges_conf_p2; range_p2];
                    else
                        ranges_cons_p2 = [ranges_cons_p2; range_p2];
                    end
                end
            end
        end
    end
end

function MRS = run_gannet_internal(mrs_f, wat_f, nii_f, ranges, label_name)
    fprintf('      Processing: %s (%d trials)...\n', label_name, size(ranges,1));
    MRS = GannetLoad_MultiSite({mrs_f}, {wat_f}, ranges);
    MRS = GannetFit(MRS);
    MRS = GannetCoRegister(MRS, {nii_f});
    try MRS = GannetSegment(MRS); catch, end
    try MRS = GannetQuantify(MRS); catch, end
    h = findobj('Type','figure', 'Name', 'GannetFit Output');
    if ~isempty(h), set(h(1), 'Name', ['Result: ' label_name]); end
end
