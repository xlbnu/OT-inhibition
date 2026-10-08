
%% OXTR grid index
ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(ft_path);
ft_defaults;
fmrs_nii = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Group_Prob_Map.nii';
oxtr_nii ='E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\MEG_mask\article\oxtr_mask\oxtr_50.nii';
cfg = [];
cfg.threshold = 0.5;
cfg.interpmethod = 'nearest';
cfg.verbose = true;
[idx_oxtr, oxtr_reference, oxtr_info] = convert_nii_roi_to_reference(oxtr_nii,fmrs_nii,cfg);

%% fMRS grid index
fmrs_nii = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Group_Prob_Map.nii';
oxtr_nii ='E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Group_Prob_Map.nii';
cfg = [];
cfg.threshold = 29;
cfg.interpmethod = 'nearest';
cfg.verbose = true;
[idx_fmrs, fmrs_reference, fmrs_info] = convert_nii_roi_to_reference(oxtr_nii,fmrs_nii,cfg);

%% MEG grid index: gamma30-100hz
load(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)\stat_cluster.mat'],'stat_res');
ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
template_grid_file = fullfile(ft_path, 'template', 'sourcemodel', 'standard_sourcemodel3d5mm.mat');
load(template_grid_file);
disp('▶ 正在加载 AAL 图谱并与网格对齐...');
atlas = ft_read_atlas(fullfile(ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));
template_grid = sourcemodel; % load(template_grid_file) The loaded variable name is sourcemodel
cfg_atlas = [];
cfg_atlas.interpmethod = 'nearest';
cfg_atlas.parameter    = 'tissue';
atlas_grid = ft_sourceinterpolate(cfg_atlas, atlas, template_grid);
cfg_cluster = [];
cfg_cluster.alpha = 0.01;
cfg_cluster.cluster_type = 'pos';
idx_cluster = getFirstSignificantCluster(stat_res, atlas, atlas_grid, template_grid, cfg_cluster);
fmrs_nii = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Group_Prob_Map.nii';
cfg = [];
% The final mapping target MEG source ROI
cfg.idx_roi = idx_cluster;
% If stat_res No unit, Current standard 5 mm grid uses cm
cfg.source_unit = 'cm';
cfg.source_coordsys = 'mni';
cfg.verbose = true;
[idx_meg_tfr, meg_roi_tfr, meg_info_tfr] =convert_meg_stat_roi_to_reference(stat_res,fmrs_nii,cfg);

%% MEG grid index: sigmaV
load(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\stat_cluster.mat'],'stat_res');
ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
template_grid_file = fullfile(ft_path, 'template', 'sourcemodel', 'standard_sourcemodel3d5mm.mat');
load(template_grid_file);
disp('▶ 正在加载 AAL 图谱并与网格对齐...');
atlas = ft_read_atlas(fullfile(ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));
template_grid = sourcemodel; % load(template_grid_file) The loaded variable name is sourcemodel
cfg_atlas = [];
cfg_atlas.interpmethod = 'nearest';
cfg_atlas.parameter    = 'tissue';
atlas_grid = ft_sourceinterpolate(cfg_atlas, atlas, template_grid);
cfg_cluster = [];
cfg_cluster.alpha = 0.05;
cfg_cluster.cluster_type = 'pos';
idx_cluster = getFirstSignificantCluster(stat_res, atlas, atlas_grid, template_grid, cfg_cluster);

fmrs_nii = 'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Group_Prob_Map.nii';

cfg = [];
% The final mapping target MEG source ROI
cfg.idx_roi = idx_cluster;
% If stat_res No unit, Current standard 5 mm grid uses cm
cfg.source_unit = 'cm';
cfg.source_coordsys = 'mni';
cfg.verbose = true;
[idx_meg_efr, meg_roi_efr, meg_info_efr] =convert_meg_stat_roi_to_reference(stat_res,fmrs_nii,cfg);

%% conjunction 
idx_fmrs_oxtr_overlap =intersect(idx_fmrs,idx_oxtr);% oxtr ∩ fMRS
idx_tfr_fmrs_overlap = intersect(idx_meg_tfr,idx_fmrs);% meg_gamma-band source ∩ fMRS
idx_efr_fmrs_overlap = intersect(idx_meg_efr,idx_fmrs);% meg_coefficient(ΣV) source ∩ fMRS
idx_three_overlap1 = intersect(intersect(idx_meg_tfr,idx_oxtr),idx_fmrs); % meg_gamma-band source ∩ oxtr & fMRS
idx_three_overlap2 = intersect(intersect(idx_meg_efr,idx_oxtr),idx_fmrs); % meg_coefficient(ΣV) source ∩ oxtr & fMRS

idx_tfr_oxtr = intersect(idx_meg_tfr,idx_oxtr);% oxtr ∩ meg_gamma-band source
idx_efr_oxtr = intersect(idx_meg_efr,idx_oxtr);% oxtr ∩ meg_coefficient(ΣV) source

oxtr_tfr1=unique([idx_tfr_fmrs_overlap;idx_fmrs_oxtr_overlap]);% (fmrs ∩ oxtr) ∪(fmrs ∩ meg_gamma-band source)
oxtr_tfr2=unique([idx_efr_fmrs_overlap;idx_fmrs_oxtr_overlap]);% (fmrs ∩ oxtr) ∪(fmrs ∩ meg_coefficient(ΣV) source)

idx_meg=unique([idx_tfr_oxtr;idx_efr_oxtr]);% oxtr ∩ (meg_coefficient(ΣV) source ∪ meg_gamma-band source)
idx_meg_fmrs_overlap=intersect(idx_meg,idx_fmrs);% fmrs ∩ (oxtr ∩ (meg_coefficient(ΣV) source ∪ meg_gamma-band source))

%%
save('E:\xianliang\matlab_m\social_decision_m\data\OT_all\MEG\figure_data\MEG_OXTR_fMRS_overlap.mat','idx_fmrs','idx_oxtr','idx_meg_tfr','idx_meg_efr');













%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx_positive, roi_reference, info] = ...
    convert_nii_roi_to_reference(input_nii, reference_nii, cfg)
% CONVERT_NII_ROI_TO_REFERENCE Resample a NIfTI ROI to a reference grid.
%
% [idx_positive, roi_reference, info] = ...
%     convert_nii_roi_to_reference(input_nii, reference_nii, cfg)
%
% Inputs
%   input_nii     : NIfTI ROI/statistical image to convert.
%   reference_nii : NIfTI defining the required output dimension, affine,
%                   orientation, field of view, and voxel size.
%   cfg.threshold : Return finite voxels with value > threshold (default 0).
%   cfg.interpmethod : 'nearest' (default) or 'linear'. Use 'nearest' for
%                      binary masks, cluster labels, and thresholded ROIs.
%   cfg.affine_tolerance_mm : Geometry-comparison tolerance (default 1e-6).
%   cfg.verbose   : true/false (default true).
%
% Outputs
%   idx_positive  : Linear indices in the REFERENCE image grid.
%   roi_reference : FieldTrip MRI structure in exactly the reference grid;
%                   ROI/statistical values are in roi_reference.anatomy.
%   info          : Geometry and conversion QC, including voxel coordinates
%                   and physical XYZ coordinates for idx_positive.
%
% Important
%   Matching voxel size alone is insufficient. Direct indexing is used only
%   when dimension, array size, and the full voxel-to-world affine all match.

if nargin < 3 || isempty(cfg)
    cfg = struct();
end

if ~isfield(cfg, 'threshold') || isempty(cfg.threshold)
    cfg.threshold = 0;
end
if ~isfield(cfg, 'interpmethod') || isempty(cfg.interpmethod)
    cfg.interpmethod = 'nearest';
end
if ~isfield(cfg, 'affine_tolerance_mm') || ...
        isempty(cfg.affine_tolerance_mm)
    cfg.affine_tolerance_mm = 1e-6;
end
if ~isfield(cfg, 'verbose') || isempty(cfg.verbose)
    cfg.verbose = true;
end

input_nii = char(string(input_nii));
reference_nii = char(string(reference_nii));

assert(isfile(input_nii), 'Input NIfTI does not exist: %s', input_nii);
assert(isfile(reference_nii), ...
    'Reference NIfTI does not exist: %s', reference_nii);
assert(isscalar(cfg.threshold) && isfinite(cfg.threshold), ...
    'cfg.threshold must be one finite scalar.');
assert(ismember(lower(string(cfg.interpmethod)), ["nearest", "linear"]), ...
    'cfg.interpmethod must be ''nearest'' or ''linear''.');

input_mri = ft_read_mri(input_nii);
reference_mri = ft_read_mri(reference_nii);

assert(isfield(input_mri, 'anatomy') && ndims(input_mri.anatomy) == 3, ...
    'Input NIfTI must contain one 3-D anatomy volume.');
assert(isfield(reference_mri, 'anatomy') && ...
       ndims(reference_mri.anatomy) == 3, ...
    'Reference NIfTI must contain one 3-D anatomy volume.');

% Compare all geometry in millimetres.
input_mri = ft_convert_units(input_mri, 'mm');
reference_mri = ft_convert_units(reference_mri, 'mm');

input_affine = double(input_mri.transform);
reference_affine = double(reference_mri.transform);

input_dim = double(input_mri.dim(:)');
reference_dim = double(reference_mri.dim(:)');

input_voxel_size = vecnorm(input_affine(1:3,1:3), 2, 1);
reference_voxel_size = vecnorm(reference_affine(1:3,1:3), 2, 1);

same_dim = isequal(input_dim, reference_dim);
same_array_size = isequal(size(input_mri.anatomy), ...
                          size(reference_mri.anatomy));
same_affine = max(abs(input_affine-reference_affine), [], 'all') <= ...
              cfg.affine_tolerance_mm;
same_geometry = same_dim && same_array_size && same_affine;

if same_geometry
    % Already in the exact reference grid: do not interpolate.
    roi_reference = input_mri;
    was_resampled = false;
else
    interp_cfg = [];
    interp_cfg.parameter = 'anatomy';
    interp_cfg.interpmethod = char(lower(string(cfg.interpmethod)));

    roi_reference = ft_sourceinterpolate( ...
        interp_cfg, input_mri, reference_mri);
    was_resampled = true;
end

% Enforce the reference geometry explicitly in the returned object.
assert(isequal(size(roi_reference.anatomy), ...
               size(reference_mri.anatomy)), ...
    'Resampled output array size does not match the reference image.');
assert(isequal(double(roi_reference.dim(:)'), reference_dim), ...
    'Resampled output dimension does not match the reference image.');

output_affine_error = max(abs( ...
    double(roi_reference.transform)-reference_affine), [], 'all');
assert(output_affine_error <= cfg.affine_tolerance_mm, ...
    'Resampled output affine does not match the reference image.');

values = double(roi_reference.anatomy);
positive_mask = isfinite(values) & values > cfg.threshold;
idx_positive = find(positive_mask);

% Convert returned linear indices to voxel and physical coordinates.
[voxel_i, voxel_j, voxel_k] = ind2sub(reference_dim, idx_positive);
voxel_ijk1 = [voxel_i, voxel_j, voxel_k, ...
              ones(numel(idx_positive),1)];
xyz = (reference_affine * voxel_ijk1')';

info = struct();
info.input_file = input_nii;
info.reference_file = reference_nii;
info.was_resampled = was_resampled;
info.interpmethod = char(lower(string(cfg.interpmethod)));
info.threshold = cfg.threshold;
info.input_dim = input_dim;
info.reference_dim = reference_dim;
info.input_voxel_size_mm = input_voxel_size;
info.reference_voxel_size_mm = reference_voxel_size;
info.input_affine = input_affine;
info.reference_affine = reference_affine;
info.same_dim = same_dim;
info.same_affine = same_affine;
info.same_geometry = same_geometry;
info.output_affine_error_mm = output_affine_error;
info.n_positive = numel(idx_positive);
info.positive_volume_mm3 = numel(idx_positive) * ...
    abs(det(reference_affine(1:3,1:3)));
info.voxel_ijk = [voxel_i, voxel_j, voxel_k];
info.xyz_mm = xyz(:,1:3);

if cfg.verbose
    fprintf('\n===================================================\n');
    fprintf('NIfTI ROI -> reference-grid conversion\n');
    fprintf('Input       : %s\n', input_nii);
    fprintf('Reference   : %s\n', reference_nii);
    fprintf('Input dim   : %s\n', mat2str(input_dim));
    fprintf('Target dim  : %s\n', mat2str(reference_dim));
    fprintf('Input voxel : %s mm\n', mat2str(input_voxel_size,4));
    fprintf('Target voxel: %s mm\n', mat2str(reference_voxel_size,4));
    fprintf('Resampled   : %s\n', string(was_resampled));
    fprintf('Threshold   : value > %.6g\n', cfg.threshold);
    fprintf('Positive voxels: %d\n', info.n_positive);
    fprintf('Positive volume: %.4f cm^3\n', ...
        info.positive_volume_mm3/1000);
    fprintf('===================================================\n\n');
end

end



function [idx_reference, roi_reference, info] = ...
    convert_meg_stat_roi_to_reference(stat_res, reference_nii, cfg)
% CONVERT_MEG_STAT_ROI_TO_REFERENCE Map a source-grid ROI to a NIfTI grid.
%
% Recommended use:
%   cfg.idx_roi = idx_cluster;  % final cluster/atlas-constrained ROI
%   [idx_reference, roi_reference, info] = ...
%       convert_meg_stat_roi_to_reference(stat_res, reference_nii, cfg);
%
% If cfg.idx_roi is empty, the function uses stat_res.(cfg.mask_field).
%
% Inputs
%   stat_res       : FieldTrip source/statistics structure with pos and dim.
%   reference_nii  : NIfTI defining the target space.
%   cfg.idx_roi     : Source-grid indices to map (recommended).
%   cfg.mask_field  : Fallback mask field (default 'mask').
%   cfg.source_unit : Used only if stat_res.unit is absent (default 'cm').
%   cfg.source_coordsys : Used if stat_res.coordsys is absent (default 'mni').
%   cfg.verbose     : true/false (default true).
%
% Outputs
%   idx_reference  : Linear indices in the reference NIfTI grid.
%   roi_reference  : FieldTrip volume; roi_reference.meg_roi is the mapped ROI.
%   info           : Geometry and conversion QC.

if nargin < 3 || isempty(cfg)
    cfg = struct();
end
if ~isfield(cfg, 'idx_roi')
    cfg.idx_roi = [];
end
if ~isfield(cfg, 'mask_field') || isempty(cfg.mask_field)
    cfg.mask_field = 'mask';
end
if ~isfield(cfg, 'source_unit') || isempty(cfg.source_unit)
    cfg.source_unit = 'cm';
end
if ~isfield(cfg, 'source_coordsys') || isempty(cfg.source_coordsys)
    cfg.source_coordsys = 'mni';
end
if ~isfield(cfg, 'verbose') || isempty(cfg.verbose)
    cfg.verbose = true;
end

reference_nii = char(string(reference_nii));
assert(isfile(reference_nii), ...
    'Reference NIfTI does not exist: %s', reference_nii);
assert(isfield(stat_res, 'pos') && isfield(stat_res, 'dim'), ...
    'stat_res must contain pos and dim.');

source_pos = double(stat_res.pos);
if size(source_pos,1) == 3 && size(source_pos,2) ~= 3
    source_pos = source_pos.';
end
assert(size(source_pos,2) == 3, ...
    'stat_res.pos must be Ngrid-by-3 or 3-by-Ngrid.');

n_grid = size(source_pos,1);
source_dim = double(stat_res.dim(:)');
assert(prod(source_dim) == n_grid, ...
    'prod(stat_res.dim) does not equal the number of source positions.');

if isempty(cfg.idx_roi)
    mask_field = char(string(cfg.mask_field));
    if ~isfield(stat_res, mask_field)
        if strcmp(mask_field, 'mask') && isfield(stat_res, 'fdr_mask')
            mask_field = 'fdr_mask';
        else
            error('stat_res does not contain mask field: %s', mask_field);
        end
    end
    mask_value = double(stat_res.(mask_field)(:));
    assert(numel(mask_value) == n_grid, ...
        'The selected stat_res mask does not match the source grid.');
    idx_roi = find(isfinite(mask_value) & mask_value ~= 0);
    roi_definition = sprintf('stat_res.%s ~= 0', mask_field);
else
    idx_roi = unique(double(cfg.idx_roi(:)));
    roi_definition = 'cfg.idx_roi';
end

assert(~isempty(idx_roi), 'The MEG source ROI is empty.');
assert(all(isfinite(idx_roi) & idx_roi == round(idx_roi) & ...
           idx_roi >= 1 & idx_roi <= n_grid), ...
    'cfg.idx_roi contains invalid source-grid indices.');

source_roi = struct();
source_roi.pos = source_pos;
source_roi.dim = source_dim;

if isfield(stat_res, 'inside') && ~isempty(stat_res.inside)
    source_roi.inside = stat_res.inside;
    inside_mask = false(n_grid,1);
    if islogical(stat_res.inside)
        inside_mask = stat_res.inside(:);
    else
        inside_mask(double(stat_res.inside(:))) = true;
    end
else
    inside_mask = true(n_grid,1);
    if isfield(stat_res, 'stat') && numel(stat_res.stat) == n_grid
        inside_mask = isfinite(double(stat_res.stat(:)));
    end
    source_roi.inside = inside_mask;
end

assert(all(inside_mask(idx_roi)), ...
    'At least one requested ROI index is outside stat_res.inside.');

if isfield(stat_res, 'unit') && ~isempty(stat_res.unit)
    source_roi.unit = stat_res.unit;
else
    source_roi.unit = char(string(cfg.source_unit));
end
if isfield(stat_res, 'coordsys') && ~isempty(stat_res.coordsys)
    source_coordsys = stat_res.coordsys;
else
    source_coordsys = char(string(cfg.source_coordsys));
end

source_roi.meg_roi = zeros(n_grid,1);
source_roi.meg_roi(~inside_mask) = NaN;
source_roi.meg_roi(idx_roi) = 1;

source_roi = ft_convert_units(source_roi, 'mm');
reference_mri = ft_read_mri(reference_nii);
reference_mri = ft_convert_units(reference_mri, 'mm');

interp_cfg = [];
interp_cfg.parameter = 'meg_roi';
interp_cfg.interpmethod = 'nearest';
roi_reference = ft_sourceinterpolate( ...
    interp_cfg, source_roi, reference_mri);

assert(isequal(size(roi_reference.meg_roi), ...
               size(reference_mri.anatomy)), ...
    'Mapped MEG ROI does not match the reference array size.');
assert(max(abs(double(roi_reference.transform) - ...
               double(reference_mri.transform)), [], 'all') < 1e-6, ...
    'Mapped MEG ROI affine does not match the reference affine.');

reference_mask = isfinite(roi_reference.meg_roi) & ...
                 roi_reference.meg_roi > 0.5;
idx_reference = find(reference_mask);

[voxel_i, voxel_j, voxel_k] = ind2sub( ...
    reference_mri.dim, idx_reference);
xyz = (double(reference_mri.transform) * ...
       [voxel_i, voxel_j, voxel_k, ...
        ones(numel(idx_reference),1)]')';

reference_voxel_volume_mm3 = abs(det( ...
    double(reference_mri.transform(1:3,1:3))));

info = struct();
info.reference_file = reference_nii;
info.roi_definition = roi_definition;
info.idx_source = idx_roi;
info.n_source_grid = n_grid;
info.n_source_roi = numel(idx_roi);
info.source_dim = source_dim;
info.source_unit = source_roi.unit;
info.source_coordsys = source_coordsys;
info.reference_dim = double(reference_mri.dim(:)');
info.reference_voxel_size_mm = vecnorm( ...
    double(reference_mri.transform(1:3,1:3)),2,1);
info.n_reference_roi = numel(idx_reference);
info.reference_volume_mm3 = ...
    numel(idx_reference) * reference_voxel_volume_mm3;
info.voxel_ijk = [voxel_i, voxel_j, voxel_k];
info.xyz_mm = xyz(:,1:3);

if cfg.verbose
    fprintf('\n===================================================\n');
    fprintf('MEG source ROI -> NIfTI reference space\n');
    fprintf('ROI definition       : %s\n', roi_definition);
    fprintf('Source grid dim      : %s\n', mat2str(source_dim));
    fprintf('Source ROI points    : %d\n', numel(idx_roi));
    fprintf('Reference dim        : %s\n', ...
        mat2str(info.reference_dim));
    fprintf('Reference voxel size : %s mm\n', ...
        mat2str(info.reference_voxel_size_mm,4));
    fprintf('Reference ROI voxels : %d\n', numel(idx_reference));
    fprintf('Reference ROI volume : %.4f cm^3\n', ...
        info.reference_volume_mm3/1000);
    fprintf('===================================================\n\n');
end

end



function [idx_cluster, cluster_info] = getFirstSignificantCluster(stat_res, atlas, atlas_grid, template_grid, cfg)

% =========================================================================
% Get the first significant component in the statistical results cluster
%
% Function: 
%   1. Automatically find the first significant cluster
%   2. Return this cluster all grid points idx_cluster
%   3. Output peak voxel 's MNI coordinates and anatomical labels
%
% Default: 
%   cfg.alpha = 0.05;
%   cfg.cluster_type = 'pos';
%
% Applicable to: 
%   FieldTrip cluster permutation results stat_res
% =========================================================================

disp('▶ 正在提取统计结果中的第一个显著 Cluster...');

if nargin < 5 || isempty(cfg)
    cfg = [];
end

if ~isfield(cfg, 'alpha')
    cfg.alpha = 0.05;
end

if ~isfield(cfg, 'cluster_type')
    cfg.cluster_type = 'pos';   % 'pos' or 'neg'
end

alpha = cfg.alpha;

% ---------------------------------------------------------
% Step one: Read positive/negative cluster information
% ---------------------------------------------------------
switch lower(cfg.cluster_type)

    case 'pos'
        assert(isfield(stat_res, 'posclusters'), ...
            'stat_res 中没有 posclusters 字段。');

        assert(isfield(stat_res, 'posclusterslabelmat'), ...
            'stat_res 中没有 posclusterslabelmat 字段。');

        cluster_prob = [stat_res.posclusters.prob];
        labelmat = stat_res.posclusterslabelmat;

    case 'neg'
        assert(isfield(stat_res, 'negclusters'), ...
            'stat_res 中没有 negclusters 字段。');

        assert(isfield(stat_res, 'negclusterslabelmat'), ...
            'stat_res 中没有 negclusterslabelmat 字段。');

        cluster_prob = [stat_res.negclusters.prob];
        labelmat = stat_res.negclusterslabelmat;

    otherwise
        error('cfg.cluster_type 必须是 ''pos'' 或 ''neg''。');
end

% ---------------------------------------------------------
% Step two: Find the first significant cluster
% ---------------------------------------------------------
sig_cluster_ids = find(cluster_prob < alpha);

if isempty(sig_cluster_ids)
    error('未找到显著 cluster：没有 cluster 的 p < %.4f。', alpha);
end

target_cluster_id = sig_cluster_ids(1);
target_cluster_p = cluster_prob(target_cluster_id);

% ---------------------------------------------------------
% Step three: Extract this cluster all grid points
% ---------------------------------------------------------
idx_cluster = find(labelmat(:) == target_cluster_id);

if isempty(idx_cluster)
    error('找到显著 cluster ID=%d，但 labelmat 中没有对应点。', target_cluster_id);
end

% ---------------------------------------------------------
% Step four: Find this cluster within peak voxel
% ---------------------------------------------------------
stat_vec = stat_res.stat(:);
cluster_stat_vals = stat_vec(idx_cluster);

switch lower(cfg.cluster_type)
    case 'pos'
        [peak_stat, peak_relative_idx] = max(cluster_stat_vals);
    case 'neg'
        [peak_stat, peak_relative_idx] = min(cluster_stat_vals);
end

peak_voxel_idx = idx_cluster(peak_relative_idx);

% MNI coordinates, FieldTrip template_grid.pos Units are typically cm, Convert to mm
peak_pos_mni = template_grid.pos(peak_voxel_idx, :) * 10;

% Anatomical labels
peak_tissue_idx = atlas_grid.tissue(peak_voxel_idx);

if peak_tissue_idx > 0 && peak_tissue_idx <= numel(atlas.tissuelabel)
    peak_exact_region = atlas.tissuelabel{peak_tissue_idx};
else
    peak_exact_region = 'Unlabeled / White Matter';
end

% ---------------------------------------------------------
% Step five: Organize outputs
% ---------------------------------------------------------
cluster_info = struct();

cluster_info.alpha = alpha;
cluster_info.cluster_type = cfg.cluster_type;
cluster_info.target_cluster_id = target_cluster_id;
cluster_info.target_cluster_p = target_cluster_p;

cluster_info.idx_cluster = idx_cluster;
cluster_info.n_cluster_points = numel(idx_cluster);

cluster_info.peak_voxel_idx = peak_voxel_idx;
cluster_info.peak_pos_mni = peak_pos_mni;
cluster_info.peak_stat = peak_stat;
cluster_info.peak_region = peak_exact_region;

% ---------------------------------------------------------
% Step six: Output report
% ---------------------------------------------------------
fprintf('\n===================================================\n');
fprintf('🎯 【第一个显著 Cluster 定位报告】\n');
fprintf('   Cluster 类型       : %s\n', cfg.cluster_type);
fprintf('   显著性阈值 alpha   : %.4f\n', alpha);
fprintf('   第一个显著 Cluster : ID = %d\n', target_cluster_id);
fprintf('   Cluster p 值       : %.6f\n', target_cluster_p);
fprintf('   Cluster 网格点数   : %d\n', numel(idx_cluster));
fprintf('---------------------------------------------------\n');
fprintf('📌 【Peak Voxel】\n');
fprintf('   解剖学定位 : %s\n', peak_exact_region);
fprintf('   MNI 坐标   : [X: %3.0f, Y: %3.0f, Z: %3.0f] mm\n', ...
    peak_pos_mni(1), peak_pos_mni(2), peak_pos_mni(3));
fprintf('   统计值     : %.4f\n', peak_stat);
fprintf('===================================================\n\n');

end
