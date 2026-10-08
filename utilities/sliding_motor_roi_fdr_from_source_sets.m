function slide_res = sliding_motor_roi_fdr_from_source_sets(source_right_by_win, source_left_by_win, cfg)
% SLIDING_MOTOR_ROI_FDR_FROM_SOURCE_SETS
%
% For multiple sliding time windows' right/left source power results perform motor ROI-FDR, and export BrainNet Usable NIfTI.
%
% Input: 
%   source_right_by_win : n_window x 1 cell
%       Each cell Contains all_source_right, That is n_subject x 1 cell, 
%       all_source_right{sub}.avg.pow For this time window right-response source power.
%
%   source_left_by_win : n_window x 1 cell
%       Each cell Contains all_source_left, That is n_subject x 1 cell.
%
%   cfg.window_start : n_window x 1, Units: seconds, For example -0.50:0.05:-0.10
%   cfg.window_end   : n_window x 1, Units: seconds, For example cfg.window_start + 0.20
%
% Required/Common cfg: 
%   cfg.output_dir
%   cfg.output_prefix
%   cfg.q = 0.05
%   cfg.roi_method = 'sphere' or 'aal'
%   cfg.left_handknob_mni_mm  = [-38 -24 56]
%   cfg.right_handknob_mni_mm = [ 38 -24 56]
%   cfg.sphere_radius_mm = 12
%   cfg.ft_path
%
% Dependencies: 
%   save_roi_fdr_source_map_for_brainnet.m
%
% Output: 
%   slide_res.window_table
%   slide_res.window(k).roi_fdr
%   Each window's *_F_masked.nii / *_signedT_masked.nii / *_cluster_label.nii and related files.

if nargin < 3
    cfg = [];
end

cfg = set_default(cfg, 'q', 0.05);
cfg = set_default(cfg, 'roi_method', 'sphere');
cfg = set_default(cfg, 'left_handknob_mni_mm', [-38, -24, 56]);
cfg = set_default(cfg, 'right_handknob_mni_mm', [38, -24, 56]);
cfg = set_default(cfg, 'sphere_radius_mm', 12);
cfg = set_default(cfg, 'output_dir', pwd);
cfg = set_default(cfg, 'output_prefix', 'sliding_beta_motor');
cfg = set_default(cfg, 'write_nifti', true);
cfg = set_default(cfg, 'write_interpolated_to_mri', true);
cfg = set_default(cfg, 'verbose', true);

n_window = numel(source_right_by_win);
assert(numel(source_left_by_win) == n_window, ...
    'source_right_by_win 和 source_left_by_win 的窗口数量不一致。');
assert(isfield(cfg, 'window_start') && isfield(cfg, 'window_end'), ...
    'cfg 必须包含 window_start 和 window_end。');
assert(numel(cfg.window_start) == n_window && numel(cfg.window_end) == n_window, ...
    'cfg.window_start/end 的长度必须等于窗口数量。');

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

slide_res = [];
slide_res.cfg = cfg;
slide_res.window = repmat(struct(), n_window, 1);

window_id = (1:n_window).';
window_start = cfg.window_start(:);
window_end = cfg.window_end(:);
window_center = (window_start + window_end) / 2;
n_fdr_all = nan(n_window, 1);
n_fdr_left = nan(n_window, 1);
n_fdr_right = nan(n_window, 1);
n_cluster = nan(n_window, 1);
crit_p = nan(n_window, 1);
out_mat_file = strings(n_window, 1);

for k = 1:n_window
    win_tag = sprintf('win_%+.3f_%+.3f', window_start(k), window_end(k));
    win_tag = strrep(win_tag, '+', 'p');
    win_tag = strrep(win_tag, '-', 'm');
    win_tag = strrep(win_tag, '.', 'p');

    cfg_k = cfg;
    cfg_k.output_prefix = sprintf('%s_%s', cfg.output_prefix, win_tag);

    if cfg.verbose
        fprintf('\n[%d/%d] ROI-FDR for window %.3f to %.3f s\n', ...
            k, n_window, window_start(k), window_end(k));
    end

    roi_fdr = save_roi_fdr_source_map_for_brainnet( ...
        source_right_by_win{k}, source_left_by_win{k}, cfg_k);

    slide_res.window(k).window_start = window_start(k);
    slide_res.window(k).window_end = window_end(k);
    slide_res.window(k).window_center = window_center(k);
    slide_res.window(k).roi_fdr = roi_fdr;

    n_fdr_all(k) = roi_fdr.summary.n_fdr_all;
    n_fdr_left(k) = roi_fdr.summary.n_fdr_left;
    n_fdr_right(k) = roi_fdr.summary.n_fdr_right;
    n_cluster(k) = roi_fdr.summary.n_clusters;
    crit_p(k) = roi_fdr.crit_p;
    if isfield(roi_fdr, 'output_mat_file')
        out_mat_file(k) = string(roi_fdr.output_mat_file);
    end
end

slide_res.window_table = table(window_id, window_start, window_end, window_center, ...
    n_fdr_all, n_fdr_left, n_fdr_right, n_cluster, crit_p, out_mat_file);

summary_file = fullfile(cfg.output_dir, [cfg.output_prefix, '_sliding_roiFDR_summary.mat']);
save(summary_file, 'slide_res', '-v7.3');
slide_res.summary_file = summary_file;

csv_file = fullfile(cfg.output_dir, [cfg.output_prefix, '_sliding_roiFDR_summary.csv']);
writetable(slide_res.window_table, csv_file);
slide_res.summary_csv = csv_file;

if cfg.verbose
    fprintf('\nSliding ROI-FDR summary saved:\n%s\n%s\n', summary_file, csv_file);
end

end


function cfg = set_default(cfg, field_name, default_value)
if ~isfield(cfg, field_name) || isempty(cfg.(field_name))
    cfg.(field_name) = default_value;
end
end

