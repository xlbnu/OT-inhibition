function roi_fdr = save_roi_fdr_source_map_for_brainnet(all_source_right, all_source_left, cfg)
% SAVE_ROI_FDR_SOURCE_MAP_FOR_BRAINNET
%
% Small-volume ROI-FDR source statistics and export for BrainNet/source plotting.
%
% This function computes voxel-wise statistics across the whole source grid:
%   diff = log(right-response source power) - log(left-response source power)
%   one-sample t-test against zero at each grid point
%   F = t^2
%
% But FDR correction is applied only inside a pre-defined motor ROI:
%   1) bilateral hand-knob sphere ROI, or
%   2) AAL motor ROI.
%
% It then saves:
%   - a .mat file with statistics, masks and cluster labels
%   - NIfTI files for BrainNet / MRIcron / FieldTrip plotting
%
% Usage:
%   cfg = [];
%   cfg.roi_method = 'sphere';
%   cfg.left_handknob_mni_mm  = [-38 -24 56];
%   cfg.right_handknob_mni_mm = [ 38 -24 56];
%   cfg.sphere_radius_mm = 12;
%   cfg.q = 0.05;
%   cfg.output_dir = 'F:\...\motor_ROI_FDR';
%   cfg.output_prefix = 'beta_peri_motor_roiFDR';
%   roi_fdr = save_roi_fdr_source_map_for_brainnet(all_source_conf, all_source_cong, cfg);
%
% Notes:
%   all_source_right = right-response source power
%   all_source_left  = left-response source power
%
% NIfTI output: *_F_masked.nii (F values at ROI-FDR significant points).
% Signed t, mean log differences, ROI/FDR masks, and cluster labels are
% retained in the MAT result, rather than exported as additional NIfTIs.
%
% Interpretation for beta ERD:
%   mean_diff > 0 means right-response beta power > left-response beta power.
%   Because beta ERD is power decrease, this often means stronger left-response
%   beta ERD at that source point.

if nargin < 3
    cfg = [];
end

cfg = set_default(cfg, 'ft_path', fileparts(which('ft_defaults')));
cfg = set_default(cfg, 'q', 0.05);
cfg = set_default(cfg, 'roi_method', 'sphere'); % 'sphere' or 'aal'
cfg = set_default(cfg, 'left_handknob_mni_mm', [-38, -24, 56]);
cfg = set_default(cfg, 'right_handknob_mni_mm', [38, -24, 56]);
cfg = set_default(cfg, 'sphere_radius_mm', 12);
cfg = set_default(cfg, 'aal_include_postcentral', false);
cfg = set_default(cfg, 'output_dir', pwd);
cfg = set_default(cfg, 'output_prefix', 'motor_roiFDR_logRight_minus_logLeft');
cfg = set_default(cfg, 'write_nifti', true);
cfg = set_default(cfg, 'write_interpolated_to_mri', true);
cfg = set_default(cfg, 'template_mri_file', fullfile(cfg.ft_path, 'template', 'anatomy', 'single_subj_T1_1mm.nii'));
cfg = set_default(cfg, 'verbose', true);

if ~exist(cfg.output_dir, 'dir')
    mkdir(cfg.output_dir);
end

num_subjects = numel(all_source_right);
if numel(all_source_left) ~= num_subjects
    error('Right and left source arrays have different participant counts.');
end

first_idx = find_first_nonempty(all_source_right);
template = all_source_right{first_idx};

if ~isfield(template, 'avg') || ~isfield(template.avg, 'pow')
    error('Source structures must contain .avg.pow.');
end
if ~isfield(template, 'pos') || ~isfield(template, 'dim')
    error('Source structures must contain .pos and .dim.');
end

num_vox = numel(template.avg.pow);
pow_diff = nan(num_subjects, num_vox);

for sub = 1:num_subjects
    if isempty(all_source_right{sub}) || isempty(all_source_left{sub})
        warning('Participant %d is empty and will be skipped.', sub);
        continue;
    end

    p_right = all_source_right{sub}.avg.pow(:);
    p_left  = all_source_left{sub}.avg.pow(:);

    if numel(p_right) ~= num_vox || numel(p_left) ~= num_vox
        error('Participant %d source grid size differs from the template.', sub);
    end

    pow_diff(sub, :) = log(max(p_right, eps)) - log(max(p_left, eps));
end

valid_subject = all(isfinite(pow_diff), 2);
pow_diff = pow_diff(valid_subject, :);
n_valid = size(pow_diff, 1);
if n_valid < 3
    error('At least three valid participants are required.');
end

if isfield(template, 'inside')
    inside = logical(template.inside(:));
else
    inside = true(num_vox, 1);
end
inside = inside & all(isfinite(pow_diff), 1).';

switch lower(cfg.roi_method)
    case 'sphere'
        [roi_mask, left_roi_mask, right_roi_mask, roi_info] = make_sphere_roi_masks(template, cfg);
    case 'aal'
        [roi_mask, left_roi_mask, right_roi_mask, roi_info] = make_aal_roi_masks(template, cfg);
    otherwise
        error('cfg.roi_method must be ''sphere'' or ''aal''.');
end

roi_mask = roi_mask(:) & inside;
left_roi_mask = left_roi_mask(:) & inside;
right_roi_mask = right_roi_mask(:) & inside;

mean_diff = nan(num_vox, 1);
tval = nan(num_vox, 1);
fval = nan(num_vox, 1);
pval = nan(num_vox, 1);

mean_diff(inside) = mean(pow_diff(:, inside), 1).';

inside_idx = find(inside);
for ii = 1:numel(inside_idx)
    v = inside_idx(ii);
    x = pow_diff(:, v);
    [t_this, p_this] = one_sample_t_two_tailed(x);
    tval(v) = t_this;
    fval(v) = t_this.^2;
    pval(v) = p_this;
end

% Small-volume FDR only inside the motor ROI.
[h_roi, crit_p, q_roi] = simple_bh_fdr(pval(roi_mask), cfg.q);

qval_roi = nan(num_vox, 1);
roi_idx = find(roi_mask);
qval_roi(roi_idx) = q_roi;

fdr_mask = false(num_vox, 1);
fdr_mask(roi_idx) = h_roi;

fdr_left_mask = fdr_mask & left_roi_mask;
fdr_right_mask = fdr_mask & right_roi_mask;

% Directional masks.
fdr_pos_mask = fdr_mask & mean_diff > 0;
fdr_neg_mask = fdr_mask & mean_diff < 0;

% Cluster labels within ROI-FDR mask using 6-neighbor connectivity in the source grid.
cluster_label = label_grid_clusters(fdr_mask, template.dim);

% Masked maps for visualization.
F_masked = zeros(num_vox, 1);
F_masked(fdr_mask) = fval(fdr_mask);

signedT_masked = zeros(num_vox, 1);
signedT_masked(fdr_mask) = tval(fdr_mask);

meanDiff_masked = zeros(num_vox, 1);
meanDiff_masked(fdr_mask) = mean_diff(fdr_mask);

qval_masked = ones(num_vox, 1);
qval_masked(roi_mask) = qval_roi(roi_mask);

roi_mask_double = double(roi_mask);

source_out = rmfield_if_exists(template, {'avg'});
source_out.stat = fval;
source_out.tstat = tval;
source_out.prob = pval;
source_out.qval_roi = qval_roi;
source_out.mean_diff = mean_diff;
source_out.roi_mask = roi_mask_double;
source_out.fdr_mask = double(fdr_mask);
source_out.fdr_left_mask = double(fdr_left_mask);
source_out.fdr_right_mask = double(fdr_right_mask);
source_out.F_masked = F_masked;
source_out.signedT_masked = signedT_masked;
source_out.meanDiff_masked = meanDiff_masked;
source_out.qval_masked = qval_masked;
source_out.cluster_label = cluster_label;

roi_fdr = [];
roi_fdr.cfg = cfg;
roi_fdr.n_subjects = n_valid;
roi_fdr.valid_subject = valid_subject;
roi_fdr.roi_info = roi_info;
roi_fdr.df1 = 1;
roi_fdr.df2 = n_valid - 1;
roi_fdr.crit_p = crit_p;
roi_fdr.source = source_out;
roi_fdr.pow_diff = pow_diff;
roi_fdr.summary = [];
roi_fdr.summary.n_roi_voxels = sum(roi_mask);
roi_fdr.summary.n_left_roi_voxels = sum(left_roi_mask);
roi_fdr.summary.n_right_roi_voxels = sum(right_roi_mask);
roi_fdr.summary.n_fdr_all = sum(fdr_mask);
roi_fdr.summary.n_fdr_left = sum(fdr_left_mask);
roi_fdr.summary.n_fdr_right = sum(fdr_right_mask);
roi_fdr.summary.n_fdr_pos = sum(fdr_pos_mask);
roi_fdr.summary.n_fdr_neg = sum(fdr_neg_mask);
roi_fdr.summary.n_clusters = max(cluster_label);

cluster_table = make_cluster_table(cluster_label, source_out, template);
roi_fdr.cluster_table = cluster_table;

mat_file = fullfile(cfg.output_dir, [cfg.output_prefix, '_roiFDR.mat']);
save(mat_file, 'roi_fdr', '-v7.3');
roi_fdr.output_mat_file = mat_file;

if cfg.verbose
    fprintf('\n============================================================\n');
    fprintf('ROI small-volume FDR source map\n');
    fprintf('ROI method        : %s\n', cfg.roi_method);
    fprintf('N subjects        : %d\n', n_valid);
    fprintf('ROI voxels        : %d\n', roi_fdr.summary.n_roi_voxels);
    fprintf('FDR q             : %.4f\n', cfg.q);
    fprintf('ROI FDR crit p    : %.10f\n', crit_p);
    fprintf('FDR voxels all    : %d\n', roi_fdr.summary.n_fdr_all);
    fprintf('FDR voxels left   : %d\n', roi_fdr.summary.n_fdr_left);
    fprintf('FDR voxels right  : %d\n', roi_fdr.summary.n_fdr_right);
    fprintf('Clusters          : %d\n', roi_fdr.summary.n_clusters);
    fprintf('Saved MAT         : %s\n', mat_file);
    fprintf('============================================================\n\n');
end

if cfg.write_nifti
    nifti_files = write_nifti_maps(source_out, cfg);
    roi_fdr.output_nifti_files = nifti_files;
    save(mat_file, 'roi_fdr', '-v7.3');
end

end


function nifti_files = write_nifti_maps(source_out, cfg)

% params = {'F_masked', 'signedT_masked', 'meanDiff_masked', ...
%           'cluster_label', 'roi_mask', 'fdr_mask'};
params = {'F_masked'};

nifti_files = struct();

if cfg.write_interpolated_to_mri
    if ~exist(cfg.template_mri_file, 'file')
        error('Template MRI not found: %s', cfg.template_mri_file);
    end

    template_mri = ft_read_mri(cfg.template_mri_file);

    interp_cfg = [];
    interp_cfg.parameter = params;
    interp_cfg.interpmethod = 'nearest';
    source_interp = ft_sourceinterpolate(interp_cfg, source_out, template_mri);

    for i = 1:numel(params)
        param = params{i};
        out_name = fullfile(cfg.output_dir, [cfg.output_prefix, '_', param, '.nii']);
        write_cfg = [];
        write_cfg.filename = strip_nii_ext(out_name);
        write_cfg.filetype = 'nifti';
        write_cfg.parameter = param;
        ft_sourcewrite(write_cfg, source_interp);
        nifti_files.(param) = out_name;
    end
else
    % Write on original source grid. This requires source_out to contain
    % correct dim/transform information. Interpolated MRI export is safer.
    for i = 1:numel(params)
        param = params{i};
        out_name = fullfile(cfg.output_dir, [cfg.output_prefix, '_', param, '.nii']);
        write_cfg = [];
        write_cfg.filename = strip_nii_ext(out_name);
        write_cfg.filetype = 'nifti';
        write_cfg.parameter = param;
        ft_sourcewrite(write_cfg, source_out);
        nifti_files.(param) = out_name;
    end
end

end


function cluster_table = make_cluster_table(cluster_label, source_out, template)

n_clusters = max(cluster_label);
if n_clusters == 0
    cluster_table = table();
    return;
end

pos = template.pos;
if max(abs(pos(:))) < 30
    pos_mm = pos * 10;
else
    pos_mm = pos;
end

cluster_id = (1:n_clusters).';
n_voxels = nan(n_clusters, 1);
peak_F = nan(n_clusters, 1);
peak_signedT = nan(n_clusters, 1);
peak_meanDiff = nan(n_clusters, 1);
peak_x = nan(n_clusters, 1);
peak_y = nan(n_clusters, 1);
peak_z = nan(n_clusters, 1);
hemi = strings(n_clusters, 1);

for c = 1:n_clusters
    idx = find(cluster_label == c);
    n_voxels(c) = numel(idx);
    [peak_F(c), rel] = max(source_out.F_masked(idx));
    peak_idx = idx(rel);
    peak_signedT(c) = source_out.signedT_masked(peak_idx);
    peak_meanDiff(c) = source_out.meanDiff_masked(peak_idx);
    peak_x(c) = pos_mm(peak_idx, 1);
    peak_y(c) = pos_mm(peak_idx, 2);
    peak_z(c) = pos_mm(peak_idx, 3);

    if mean(pos_mm(idx, 1)) < 0
        hemi(c) = "left";
    elseif mean(pos_mm(idx, 1)) > 0
        hemi(c) = "right";
    else
        hemi(c) = "midline";
    end
end

cluster_table = table(cluster_id, hemi, n_voxels, peak_F, peak_signedT, ...
    peak_meanDiff, peak_x, peak_y, peak_z);

end


function label_vec = label_grid_clusters(mask_vec, dim)

mask_vol = reshape(logical(mask_vec), dim);
label_vol = zeros(dim);
current_label = 0;

sx = dim(1);
sy = dim(2);
sz = dim(3);

for x = 1:sx
    for y = 1:sy
        for z = 1:sz
            if mask_vol(x, y, z) && label_vol(x, y, z) == 0
                current_label = current_label + 1;
                queue = zeros(nnz(mask_vol), 3);
                q_start = 1;
                q_end = 1;
                queue(q_end, :) = [x, y, z];
                label_vol(x, y, z) = current_label;

                while q_start <= q_end
                    p = queue(q_start, :);
                    q_start = q_start + 1;

                    neigh = [ ...
                        p(1)-1, p(2),   p(3); ...
                        p(1)+1, p(2),   p(3); ...
                        p(1),   p(2)-1, p(3); ...
                        p(1),   p(2)+1, p(3); ...
                        p(1),   p(2),   p(3)-1; ...
                        p(1),   p(2),   p(3)+1];

                    for ni = 1:6
                        nx = neigh(ni, 1);
                        ny = neigh(ni, 2);
                        nz = neigh(ni, 3);
                        if nx >= 1 && nx <= sx && ny >= 1 && ny <= sy && nz >= 1 && nz <= sz
                            if mask_vol(nx, ny, nz) && label_vol(nx, ny, nz) == 0
                                label_vol(nx, ny, nz) = current_label;
                                q_end = q_end + 1;
                                queue(q_end, :) = [nx, ny, nz];
                            end
                        end
                    end
                end
            end
        end
    end
end

label_vec = label_vol(:);

end


function [roi_mask, left_mask, right_mask, roi_info] = make_sphere_roi_masks(template, cfg)

pos = template.pos;
if max(abs(pos(:))) < 30
    pos_mm = pos * 10;
else
    pos_mm = pos;
end

if isfield(template, 'inside')
    inside = logical(template.inside(:));
else
    inside = true(size(pos_mm, 1), 1);
end

d_left = sqrt(sum((pos_mm - cfg.left_handknob_mni_mm).^2, 2));
d_right = sqrt(sum((pos_mm - cfg.right_handknob_mni_mm).^2, 2));

left_mask = d_left <= cfg.sphere_radius_mm & inside;
right_mask = d_right <= cfg.sphere_radius_mm & inside;
roi_mask = left_mask | right_mask;

roi_info = [];
roi_info.method = 'sphere';
roi_info.left_handknob_mni_mm = cfg.left_handknob_mni_mm;
roi_info.right_handknob_mni_mm = cfg.right_handknob_mni_mm;
roi_info.sphere_radius_mm = cfg.sphere_radius_mm;

end


function [roi_mask, left_mask, right_mask, roi_info] = make_aal_roi_masks(template, cfg)

atlas_file = fullfile(cfg.ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
atlas = ft_read_atlas(atlas_file);

interp_cfg = [];
interp_cfg.interpmethod = 'nearest';
interp_cfg.parameter = 'tissue';
atlas_grid = ft_sourceinterpolate(interp_cfg, atlas, template);

left_labels = {'Precentral_L'};
right_labels = {'Precentral_R'};
if cfg.aal_include_postcentral
    left_labels{end+1} = 'Postcentral_L';
    right_labels{end+1} = 'Postcentral_R';
end

left_ids = find_label_ids(atlas.tissuelabel, left_labels);
right_ids = find_label_ids(atlas.tissuelabel, right_labels);

left_mask = false(size(atlas_grid.tissue(:)));
right_mask = false(size(atlas_grid.tissue(:)));

for i = 1:numel(left_ids)
    left_mask = left_mask | atlas_grid.tissue(:) == left_ids(i);
end
for i = 1:numel(right_ids)
    right_mask = right_mask | atlas_grid.tissue(:) == right_ids(i);
end

if isfield(template, 'inside')
    inside = logical(template.inside(:));
else
    inside = true(size(left_mask));
end

left_mask = left_mask & inside;
right_mask = right_mask & inside;
roi_mask = left_mask | right_mask;

roi_info = [];
roi_info.method = 'aal';
roi_info.atlas_file = atlas_file;
roi_info.left_labels = left_labels;
roi_info.right_labels = right_labels;
roi_info.left_tissue_ids = left_ids;
roi_info.right_tissue_ids = right_ids;

end


function ids = find_label_ids(labels, target_labels)

ids = [];
for i = 1:numel(target_labels)
    target = target_labels{i};
    idx = find(strcmp(labels, target));
    if isempty(idx)
        idx = find(contains(labels, target, 'IgnoreCase', true));
    end
    if isempty(idx)
        error('AAL label not found: %s', target);
    end
    ids = [ids; idx(:)]; %#ok<AGROW>
end
ids = unique(ids);

end


function [tval, pval] = one_sample_t_two_tailed(x)

x = x(:);
x = x(isfinite(x));
n = numel(x);
if n < 3
    tval = nan;
    pval = nan;
    return;
end

m = mean(x);
s = std(x, 0);
if s == 0
    if m == 0
        tval = 0;
        pval = 1;
    else
        tval = sign(m) * inf;
        pval = 0;
    end
    return;
end

tval = m / (s / sqrt(n));
pval = 2 * tcdf(-abs(tval), n - 1);

end


function [h, crit_p, adj_p] = simple_bh_fdr(pvals, q)

p = pvals(:);
valid = isfinite(p);
pv = p(valid);

h = false(size(p));
adj_p = nan(size(p));
crit_p = 0;

if isempty(pv)
    return;
end

[ps, sort_idx] = sort(pv, 'ascend');
m = numel(ps);

thresh = ((1:m)' / m) * q;
below = ps <= thresh;

h_valid = false(size(pv));
if any(below)
    k = find(below, 1, 'last');
    crit_p = ps(k);
    h_valid = pv <= crit_p;
end

adj_sorted = ps .* m ./ (1:m)';
adj_sorted = flipud(cummin(flipud(adj_sorted)));
adj_sorted(adj_sorted > 1) = 1;

adj_valid = nan(size(pv));
adj_valid(sort_idx) = adj_sorted;

h(valid) = h_valid;
adj_p(valid) = adj_valid;

end


function out = rmfield_if_exists(in, fields)

out = in;
for i = 1:numel(fields)
    if isfield(out, fields{i})
        out = rmfield(out, fields{i});
    end
end

end


function s = strip_nii_ext(filename)

[p, n, e] = fileparts(filename);
if strcmpi(e, '.nii')
    s = fullfile(p, n);
else
    s = filename;
end

end


function cfg = set_default(cfg, field_name, default_value)

if ~isfield(cfg, field_name) || isempty(cfg.(field_name))
    cfg.(field_name) = default_value;
end

end


function idx = find_first_nonempty(x)

idx = [];
for i = 1:numel(x)
    if ~isempty(x{i})
        idx = i;
        return;
    end
end
error('All input source cells are empty.');

end

