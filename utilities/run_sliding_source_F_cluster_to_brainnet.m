function summary_table = run_sliding_source_F_cluster_to_brainnet(user_cfg)
%RUN_SLIDING_SOURCE_F_CLUSTER_TO_BRAINNET
% For sliding time windows' DICS source power perform paired tests window by window F cluster-permutation, 
% and output BrainNet Viewer readable NIfTI map files.
%
% Default call: 
%   summary_table = run_sliding_source_F_cluster_to_brainnet();
%
% Primary input file(Within each sliding-window folder): 
%   all_source.mat
%     all_source_conf : 19 x 1 cell
%     all_source_cong : 19 x 1 cell
%
% Each time window outputs four types of NIfTI: 
%   1) F_all                       Unthresholded F value
%   2) F_clusterSig                cluster Significant after correction F value
%   3) SignedLogRatio_clusterSig   Group mean at significant locations log(conflict)-log(congruent)
%   4) ClusterMask                 Binary significance mask
%
% Notes: 
%   - Each window receives a separate spatial cluster correction, No multiple-comparison correction across time windows.
%   - F values have no direction; Direction is determined by SignedLogRatio_clusterSig to represent.
%   - Continuous statistics use linear Interpolation, Binary mask Use nearest Interpolation.

if nargin < 1 || isempty(user_cfg)
    user_cfg = struct();
end

%% ----------------------------- parameters -----------------------------
cfg = struct();
cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg.root_dir = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MEG_source_408IndependDay\decision_30-100hz'];
cfg.output_dir = fullfile(cfg.root_dir, ...
    'BrainNet_pairedF_cluster_AALcortex_5000perm');

cfg.expected_subjects = 19;
cfg.expected_starts = -1.00:0.05:-0.10;
cfg.window_length = 0.20;
cfg.window_tolerance = 1e-8;

cfg.numrandomization = 5000;
cfg.clusteralpha = 0.01;
cfg.alpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0;
cfg.randomseed = 20260814;

cfg.resume = true;             % Resume completed statistical windows from checkpoint continue
cfg.overwrite_nifti = true;    % checkpoint Regenerate when present NIfTI
cfg.dry_run = false;           % true: Check files and data only, Do not run statistics

% Construct"Cortical only"AAL mask; Lingual_L/R Not in the exclusion list, Therefore retained.
cfg.exclude_aal_keywords = {
    'Cerebellum', 'Vermis', 'Thalamus', 'Caudate', 'Putamen', ...
    'Pallidum', 'Amygdala', 'Hippocampus', 'Olfactory'};

cfg = overwrite_struct(cfg, user_cfg);

assert(isfolder(cfg.ft_path), 'FieldTrip 路径不存在：%s', cfg.ft_path);
assert(isfolder(cfg.root_dir), '源结果路径不存在：%s', cfg.root_dir);

addpath(cfg.ft_path);
ft_defaults;

%% ------------------------- Find and sort windows -------------------------
window_table = discover_source_windows(cfg.root_dir);
validate_window_sequence(window_table, cfg);

fprintf('\n============================================================\n');
fprintf('Sliding-window source paired-F cluster analysis\n');
fprintf('Root             : %s\n', cfg.root_dir);
fprintf('Number of windows: %d\n', height(window_table));
fprintf('Subjects/window  : %d\n', cfg.expected_subjects);
fprintf('Randomizations   : %d\n', cfg.numrandomization);
fprintf('============================================================\n\n');

%% -------------------- Standard grid, AAL and cortical mask --------------------
grid_file = fullfile(cfg.ft_path, 'template', 'sourcemodel', ...
    'standard_sourcemodel3d5mm.mat');
grid_loaded = load(grid_file);
assert(isfield(grid_loaded, 'sourcemodel'), ...
    '模板文件中没有变量 sourcemodel：%s', grid_file);
template_grid = grid_loaded.sourcemodel;
template_grid = ft_convert_units(template_grid, 'mm');

n_pos = size(template_grid.pos, 1);
template_inside = inside_to_mask(template_grid.inside, n_pos);

atlas_file = fullfile(cfg.ft_path, 'template', 'atlas', 'aal', ...
    'ROI_MNI_V4.nii');
atlas = ft_read_atlas(atlas_file);

cfg_interp = [];
cfg_interp.interpmethod = 'nearest';
cfg_interp.parameter = 'tissue';
atlas_grid = ft_sourceinterpolate(cfg_interp, atlas, template_grid);
assert(numel(atlas_grid.tissue) == n_pos, ...
    'AAL 插值后体素数与模板网格不一致。');

cortical_mask = false(n_pos, 1);
for i_roi = 1:numel(atlas.tissuelabel)
    roi_name = atlas.tissuelabel{i_roi};
    is_excluded = any(contains(roi_name, cfg.exclude_aal_keywords, ...
        'IgnoreCase', true));
    if ~is_excluded
        cortical_mask(atlas_grid.tissue(:) == i_roi) = true;
    end
end
final_inside = template_inside & cortical_mask;

fprintf('Template inside points : %d\n', nnz(template_inside));
fprintf('AAL cortical points    : %d\n', nnz(final_inside));
fprintf('Lingual_L/R retained   : yes\n\n');
assert(nnz(final_inside) > 0, 'AAL 皮层 mask 为空。');

% BrainNet NIfTI target space.
template_mri_file = fullfile(cfg.ft_path, 'template', 'anatomy', ...
    'single_subj_T1_1mm.nii');
template_mri = ft_read_mri(template_mri_file);
template_mri = ft_convert_units(template_mri, 'mm');

%% --------------------------- Output directory ---------------------------
output_dirs = struct();
output_dirs.mat = fullfile(cfg.output_dir, 'statistics_mat');
output_dirs.F_all = fullfile(cfg.output_dir, 'NIfTI_F_all');
output_dirs.F_sig = fullfile(cfg.output_dir, 'NIfTI_F_clusterSig');
output_dirs.effect_sig = fullfile(cfg.output_dir, ...
    'NIfTI_SignedLogRatio_clusterSig');
output_dirs.mask = fullfile(cfg.output_dir, 'NIfTI_ClusterMask');

dir_names = fieldnames(output_dirs);
for i_dir = 1:numel(dir_names)
    if ~isfolder(output_dirs.(dir_names{i_dir}))
        mkdir(output_dirs.(dir_names{i_dir}));
    end
end

analysis_signature = sprintf([...
    'pairedF|cluster|maxsum|nperm=%d|clusteralpha=%.8g|' ...
    'alpha=%.8g|minnbchan=%d|AALcortex_keepLingual|v1'], ...
    cfg.numrandomization, cfg.clusteralpha, cfg.alpha, cfg.minnbchan);

%% ------------------------- Analyze each time window -------------------------
summary_table = table();

for i_win = 1:height(window_table)
    t_start = window_table.start_time(i_win);
    t_end = window_table.end_time(i_win);
    t_center = mean([t_start, t_end]);
    input_file = window_table.mat_file{i_win};

    frame_tag = sprintf('frame_%02d_t_%+.2f_%+.2f', ...
        i_win, t_start, t_end);
    frame_tag = strrep(frame_tag, '+', 'p');
    frame_tag = strrep(frame_tag, '-', 'm');
    % The filename stem must not contain decimal points, Otherwise ft_volumewrite/fileparts interprets the final
    % decimal segment as an extension, Causing NIfTI format detection to fail.
    frame_tag = strrep(frame_tag, '.', 'p');

    stat_file = fullfile(output_dirs.mat, [frame_tag '_stat.mat']);

    fprintf('\n[%02d/%02d] window [%+.2f, %+.2f] s\n', ...
        i_win, height(window_table), t_start, t_end);
    fprintf('Input: %s\n', input_file);

    if cfg.resume && isfile(stat_file)
        checkpoint = load(stat_file, 'stat_res', 'mean_logratio', ...
            'binary_mask_source', ...
            'analysis_signature_saved');
        assert(strcmp(checkpoint.analysis_signature_saved, analysis_signature), ...
            ['Checkpoint 参数与当前参数不一致：%s\n' ...
             '请更换输出目录或关闭 resume。'], stat_file);

        stat_res = checkpoint.stat_res;
        mean_logratio = checkpoint.mean_logratio;
        binary_mask_source = checkpoint.binary_mask_source;
        fprintf('Loaded checkpoint: %s\n', stat_file);
    else
        data_loaded = load(input_file, ...
            'all_source_conf', 'all_source_cong');
        assert(isfield(data_loaded, 'all_source_conf') && ...
               isfield(data_loaded, 'all_source_cong'), ...
            '文件缺少 all_source_conf 或 all_source_cong：%s', input_file);

        all_source_conf = data_loaded.all_source_conf;
        all_source_cong = data_loaded.all_source_cong;
        validate_source_cells(all_source_conf, all_source_cong, ...
            template_grid, cfg.expected_subjects, input_file);

        log_source_conf = cell(cfg.expected_subjects, 1);
        log_source_cong = cell(cfg.expected_subjects, 1);
        logratio_subject = nan(cfg.expected_subjects, n_pos);

        for i_sub = 1:cfg.expected_subjects
            conf_pow = double(all_source_conf{i_sub}.avg.pow(:));
            cong_pow = double(all_source_cong{i_sub}.avg.pow(:));

            invalid_inside = final_inside & (...
                ~isfinite(conf_pow) | ~isfinite(cong_pow) | ...
                conf_pow <= 0 | cong_pow <= 0);
            assert(~any(invalid_inside), ...
                ['window [%+.2f,%+.2f], subject %d：' ...
                 '皮层 mask 内存在非正值或非有限功率。'], ...
                t_start, t_end, i_sub);

            conf_log = nan(n_pos, 1);
            cong_log = nan(n_pos, 1);
            conf_log(final_inside) = log(conf_pow(final_inside))-mean(log(conf_pow(final_inside)),'omitnan');
            cong_log(final_inside) = log(cong_pow(final_inside))-mean(log(cong_pow(final_inside)),'omitnan');

            logratio_subject(i_sub, :) = (conf_log - cong_log).';

            log_source_conf{i_sub} = make_stat_source(...
                all_source_conf{i_sub}, conf_log, template_grid, final_inside);
            log_source_cong{i_sub} = make_stat_source(...
                all_source_cong{i_sub}, cong_log, template_grid, final_inside);
        end

        mean_logratio = mean(logratio_subject, 1, 'omitnan').';

        % Paired design: 19 participants x Two conditions.
        design = zeros(2, 2 * cfg.expected_subjects);
        design(1, :) = [1:cfg.expected_subjects, 1:cfg.expected_subjects];
        design(2, :) = [ones(1, cfg.expected_subjects), ...
                        2 * ones(1, cfg.expected_subjects)];

        cfg_stat = [];
        cfg_stat.dim = template_grid.dim;
        cfg_stat.method = 'montecarlo';
        cfg_stat.statistic = 'ft_statfun_depsamplesFunivariate';
        cfg_stat.parameter = 'avg.pow';
        cfg_stat.correctm = 'cluster';
        cfg_stat.clusteralpha = cfg.clusteralpha;
        cfg_stat.clusterstatistic = cfg.clusterstatistic;
        cfg_stat.minnbchan = cfg.minnbchan;
        cfg_stat.tail = 1;
        cfg_stat.clustertail = 1;
        cfg_stat.alpha = cfg.alpha;
        cfg_stat.numrandomization = cfg.numrandomization;
        cfg_stat.randomseed = cfg.randomseed;
        cfg_stat.design = design;
        cfg_stat.uvar = 1;
        cfg_stat.ivar = 2;

        if cfg.dry_run
            fprintf('Dry run passed; statistics were not executed.\n');
            continue;
        end

        stat_res = ft_sourcestatistics(cfg_stat, ...
            log_source_conf{:}, log_source_cong{:});

        assert(numel(stat_res.stat) == n_pos, ...
            '统计输出点数与模板网格不一致。');

        if isfield(stat_res, 'mask') && ~isempty(stat_res.mask)
            raw_mask = double(stat_res.mask(:));
            % FieldTrip Typically outside source grid returns NaN; Only finite nonzero points
            % can be considered significant, Cannot directly logical(NaN).
            binary_mask_source = isfinite(raw_mask) & raw_mask > 0.5;
        else
            binary_mask_source = false(n_pos, 1);
        end
        binary_mask_source = binary_mask_source & final_inside;

        window_meta = struct();
        window_meta.index = i_win;
        window_meta.start_time = t_start;
        window_meta.end_time = t_end;
        window_meta.center_time = t_center;
        window_meta.input_file = input_file;

        analysis_signature_saved = analysis_signature; %#ok<NASGU>
        save(stat_file, 'stat_res', 'mean_logratio', ...
            'binary_mask_source', 'window_meta', ...
            'analysis_signature_saved', 'cfg', '-v7.3');
        fprintf('Saved checkpoint: %s\n', stat_file);
    end

    if cfg.dry_run
        continue;
    end

    %% Prepare source-space parameters
    F_all_source = double(stat_res.stat(:));
    F_all_source(~isfinite(F_all_source) | ~final_inside) = 0;

    mean_logratio = double(mean_logratio(:));
    mean_logratio(~isfinite(mean_logratio) | ~final_inside) = 0;

    binary_mask_source = isfinite(double(binary_mask_source(:))) & ...
        double(binary_mask_source(:)) > 0.5 & final_inside;

    %% Continuous parameters linear Interpolation; Binary mask nearest Interpolation
    source_cont = make_interpolation_source(template_grid, template_inside);
    source_cont.F_all = F_all_source;
    source_cont.mean_logratio = mean_logratio;

    cfg_cont = [];
    cfg_cont.parameter = {'F_all', 'mean_logratio'};
    cfg_cont.interpmethod = 'linear';
    interp_cont = ft_sourceinterpolate(cfg_cont, source_cont, template_mri);

    source_mask = make_interpolation_source(template_grid, template_inside);
    source_mask.binary_mask = double(binary_mask_source);

    cfg_mask = [];
    cfg_mask.parameter = 'binary_mask';
    cfg_mask.interpmethod = 'nearest';
    interp_mask = ft_sourceinterpolate(cfg_mask, source_mask, template_mri);

    F_all_vol = double(interp_cont.F_all);
    mean_logratio_vol = double(interp_cont.mean_logratio);
    mask_vol = logical(interp_mask.binary_mask >= 0.5);

    F_all_vol(~isfinite(F_all_vol)) = 0;
    mean_logratio_vol(~isfinite(mean_logratio_vol)) = 0;

    F_sig_vol = F_all_vol;
    F_sig_vol(~mask_vol) = 0;

    signed_logratio_sig_vol = mean_logratio_vol;
    signed_logratio_sig_vol(~mask_vol) = 0;

    %% Write BrainNet NIfTI
    file_F_all = fullfile(output_dirs.F_all, [frame_tag '_F_all']);
    file_F_sig = fullfile(output_dirs.F_sig, [frame_tag '_F_clusterSig']);
    file_effect = fullfile(output_dirs.effect_sig, ...
        [frame_tag '_SignedLogRatio_clusterSig']);
    file_mask = fullfile(output_dirs.mask, [frame_tag '_ClusterMask']);

    write_brainnet_nifti(template_mri, F_all_vol, 'F_all', ...
        file_F_all, 'single', cfg.overwrite_nifti);
    write_brainnet_nifti(template_mri, F_sig_vol, 'F_clusterSig', ...
        file_F_sig, 'single', cfg.overwrite_nifti);
    write_brainnet_nifti(template_mri, signed_logratio_sig_vol, ...
        'SignedLogRatio_clusterSig', file_effect, 'single', ...
        cfg.overwrite_nifti);
    write_brainnet_nifti(template_mri, uint8(mask_vol), 'ClusterMask', ...
        file_mask, 'uint8', cfg.overwrite_nifti);

    %% Aggregate statistics
    row = summarize_window(stat_res, binary_mask_source, mean_logratio, ...
        final_inside, template_grid, atlas_grid, atlas, ...
        i_win, t_start, t_end, input_file, cfg.alpha);
    summary_table = [summary_table; row]; %#ok<AGROW>

    fprintf('Significant source points: %d; significant clusters: %d\n', ...
        row.n_sig_source_points, row.n_sig_clusters);
end

if cfg.dry_run
    fprintf('\nDry run completed. No statistics or NIfTI files were generated.\n');
    return;
end

summary_file = fullfile(cfg.output_dir, ...
    'sliding_window_cluster_summary.csv');
writetable(summary_table, summary_file);
save(fullfile(cfg.output_dir, 'analysis_settings.mat'), ...
    'cfg', 'window_table', 'final_inside', 'analysis_signature', '-v7.3');

fprintf('\n============================================================\n');
fprintf('All windows finished.\n');
fprintf('Summary: %s\n', summary_file);
fprintf('NIfTI root: %s\n', cfg.output_dir);
fprintf('============================================================\n');
end


%% ========================================================================
function out = overwrite_struct(defaults, user_cfg)
out = defaults;
names = fieldnames(user_cfg);
for i = 1:numel(names)
    out.(names{i}) = user_cfg.(names{i});
end
end


function window_table = discover_source_windows(root_dir)
D = dir(fullfile(root_dir, 'Source_*'));
D = D([D.isdir]);

start_time = [];
end_time = [];
folder = {};
mat_file = {};

% pattern = '65hz\(([-+]?\d*\.?\d+)_([-+]?\d*\.?\d+)\)';
pattern = '\(([-+]?\d*\.?\d+)_([-+]?\d*\.?\d+)\)$';
for i = 1:numel(D)
    token = regexp(D(i).name, pattern, 'tokens', 'once');
    if isempty(token)
        continue;
    end
    this_file = fullfile(D(i).folder, D(i).name, 'all_source.mat');
    assert(isfile(this_file), ...
        '时间窗文件夹中缺少 all_source.mat：%s', ...
        fullfile(D(i).folder, D(i).name));

    start_time(end+1, 1) = str2double(token{1}); %#ok<AGROW>
    end_time(end+1, 1) = str2double(token{2}); %#ok<AGROW>
    folder{end+1, 1} = fullfile(D(i).folder, D(i).name); %#ok<AGROW>
    mat_file{end+1, 1} = this_file; %#ok<AGROW>
end

assert(~isempty(start_time), '未识别到任何滑动时间窗文件夹。');
window_table = table(start_time, end_time, folder, mat_file);
window_table = sortrows(window_table, {'start_time', 'end_time'});
end


function validate_window_sequence(window_table, cfg)
expected_starts = cfg.expected_starts(:);
expected_ends = expected_starts + cfg.window_length;

assert(height(window_table) == numel(expected_starts), ...
    '应有 %d 个时间窗，但识别到 %d 个。', ...
    numel(expected_starts), height(window_table));

for i = 1:numel(expected_starts)
    assert(abs(window_table.start_time(i) - expected_starts(i)) ...
        <= cfg.window_tolerance, ...
        '第 %d 个窗口起始时间异常：实际 %.10g，预期 %.10g。', ...
        i, window_table.start_time(i), expected_starts(i));
    assert(abs(window_table.end_time(i) - expected_ends(i)) ...
        <= cfg.window_tolerance, ...
        '第 %d 个窗口终止时间异常：实际 %.10g，预期 %.10g。', ...
        i, window_table.end_time(i), expected_ends(i));
end
end


function mask = inside_to_mask(inside, n_pos)
if islogical(inside)
    mask = inside(:);
    assert(numel(mask) == n_pos, 'logical inside 长度与 pos 不一致。');
else
    mask = false(n_pos, 1);
    inside = inside(:);
    assert(all(inside >= 1 & inside <= n_pos), 'inside 索引超出范围。');
    mask(inside) = true;
end
end


function validate_source_cells(conf, cong, template_grid, expected_n, input_file)
assert(iscell(conf) && iscell(cong), ...
    'all_source_conf/cong 必须是 cell：%s', input_file);
assert(numel(conf) == expected_n && numel(cong) == expected_n, ...
    '被试数错误：%s；conf=%d，cong=%d，预期=%d。', ...
    input_file, numel(conf), numel(cong), expected_n);

n_pos = size(template_grid.pos, 1);
for i = 1:expected_n
    assert(~isempty(conf{i}) && ~isempty(cong{i}), ...
        'subject %d source 为空：%s', i, input_file);
    assert(isfield(conf{i}, 'avg') && isfield(conf{i}.avg, 'pow') && ...
           isfield(cong{i}, 'avg') && isfield(cong{i}.avg, 'pow'), ...
        'subject %d 缺少 avg.pow：%s', i, input_file);
    assert(numel(conf{i}.avg.pow) == n_pos && ...
           numel(cong{i}.avg.pow) == n_pos, ...
        'subject %d 的 pow 点数与模板不一致：%s', i, input_file);
    assert(isfield(conf{i}, 'dim') && isfield(cong{i}, 'dim') && ...
           isequal(double(conf{i}.dim(:)), double(template_grid.dim(:))) && ...
           isequal(double(cong{i}.dim(:)), double(template_grid.dim(:))), ...
        'subject %d 的 dim 与模板不一致：%s', i, input_file);
end
end


function source_out = make_stat_source(source_in, parameter, template_grid, inside)
source_out = source_in;
if isfield(source_out, 'cfg')
    source_out = rmfield(source_out, 'cfg');
end
source_out.pos = template_grid.pos;
source_out.dim = template_grid.dim;
source_out.inside = inside;
source_out.unit = template_grid.unit;
if isfield(template_grid, 'coordsys')
    source_out.coordsys = template_grid.coordsys;
end
source_out.avg.pow = parameter(:);
end


function source_out = make_interpolation_source(template_grid, inside)
source_out = [];
source_out.pos = template_grid.pos;
source_out.dim = template_grid.dim;
source_out.inside = inside;
source_out.unit = template_grid.unit;
if isfield(template_grid, 'coordsys')
    source_out.coordsys = template_grid.coordsys;
end
end


function write_brainnet_nifti(template_mri, data, parameter_name, ...
    filename_no_ext, datatype, overwrite)
output_file = [filename_no_ext '.nii'];
if isfile(output_file) && ~overwrite
    fprintf('Skip existing NIfTI: %s\n', output_file);
    return;
end

assert(isequal(size(data), double(template_mri.dim(:).')), ...
    'NIfTI 数据维度与模板 MRI 不一致：%s', filename_no_ext);

volume = [];
volume.dim = template_mri.dim;
volume.transform = template_mri.transform;
volume.unit = template_mri.unit;
if isfield(template_mri, 'coordsys')
    volume.coordsys = template_mri.coordsys;
end
volume.(parameter_name) = data;

cfg_write = [];
cfg_write.filename = filename_no_ext; % ft_volumewrite Requires no extension
cfg_write.filetype = 'nifti';
cfg_write.parameter = parameter_name;
cfg_write.datatype = datatype;
cfg_write.scaling = 'no';
ft_volumewrite(cfg_write, volume);

assert(isfile(output_file), 'NIfTI 写出失败：%s', output_file);
end


function row = summarize_window(stat_res, sig_mask, mean_logratio, ...
    final_inside, template_grid, atlas_grid, atlas, ...
    i_win, t_start, t_end, input_file, alpha)
F = double(stat_res.stat(:));
F_valid = F;
F_valid(~final_inside | ~isfinite(F_valid)) = -Inf;
[max_F_all, peak_all_idx] = max(F_valid);

if any(sig_mask)
    F_sig = F;
    F_sig(~sig_mask | ~isfinite(F_sig)) = -Inf;
    [max_F_sig, peak_sig_idx] = max(F_sig);
    peak_pos = template_grid.pos(peak_sig_idx, :);
    peak_effect = mean_logratio(peak_sig_idx);
    tissue_id = atlas_grid.tissue(peak_sig_idx);
    if tissue_id >= 1 && tissue_id <= numel(atlas.tissuelabel)
        peak_region = string(atlas.tissuelabel{tissue_id});
    else
        peak_region = "Unlabeled";
    end
else
    max_F_sig = NaN;
    peak_sig_idx = NaN;
    peak_pos = [NaN NaN NaN];
    peak_effect = NaN;
    peak_region = "No significant cluster";
end

n_sig_clusters = 0;
min_cluster_p = NaN;
if isfield(stat_res, 'posclusters') && ~isempty(stat_res.posclusters)
    cluster_p = [stat_res.posclusters.prob];
    n_sig_clusters = nnz(cluster_p <= alpha);
    if ~isempty(cluster_p)
        min_cluster_p = min(cluster_p);
    end
end

row = table(i_win, t_start, t_end, mean([t_start t_end]), ...
    string(input_file), nnz(final_inside), nnz(sig_mask), ...
    n_sig_clusters, min_cluster_p, max_F_all, peak_all_idx, ...
    max_F_sig, peak_sig_idx, peak_pos(1), peak_pos(2), peak_pos(3), ...
    peak_effect, peak_region, ...
    'VariableNames', {'window_index', 'start_time_s', 'end_time_s', ...
    'center_time_s', 'input_file', 'n_tested_source_points', ...
    'n_sig_source_points', 'n_sig_clusters', 'min_cluster_p', ...
    'max_F_all', 'peak_all_source_index', 'max_F_sig', ...
    'peak_sig_source_index', 'peak_x_mm', 'peak_y_mm', 'peak_z_mm', ...
    'peak_signed_logratio', 'peak_AAL_region'});
end
