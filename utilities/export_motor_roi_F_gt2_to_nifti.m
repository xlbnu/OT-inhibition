function out = export_motor_roi_F_gt2_to_nifti( ...
    stat_res,mri_ref_file,output_dir,cfg)
% EXPORT_MOTOR_ROI_F_GT2_TO_NIFTI
%
% Restrict source-level F statistic to the precentral gyrus/postcentral gyrus, and retain
% ROI all within F > cfg.F_threshold connected components, Export as BrainNet readable NIfTI.
%
% Input: 
%   stat_res       FieldTrip source statistics structure, Must contain stat/pos/dim/inside
%   mri_ref_file   For output space and affine 's MNI MRI NIfTI
%   output_dir     Output directory
%   cfg            Optional settings: 
%       .ft_path       FieldTrip path
%       .aal_file      AAL NIfTI path
%       .roi_mode      'm1' or 'sensorimotor', Default sensorimotor
%       .F_threshold   Default 2
%       .connectivity  Default 26
%
% Output: 
%   out.F_nifti       F>threshold 's NIfTI
%   out.cluster_nifti Connected-component ID NIfTI
%   out.cluster_table Size of each component, Peak value and coordinates
%
% Note: F>2 is an uncorrected exploratory threshold, Does not represent cluster-level significance.

if nargin < 4 || isempty(cfg), cfg = struct(); end

cfg = default_local(cfg,'ft_path', ...
    'D:\software\matlab_toolbox\fieldtrip-20251218');
cfg = default_local(cfg,'aal_file',fullfile(cfg.ft_path, ...
    'template','atlas','aal','ROI_MNI_V4.nii'));
cfg = default_local(cfg,'roi_mode','sensorimotor');
cfg = default_local(cfg,'F_threshold',2);
cfg = default_local(cfg,'connectivity',26);

assert(isfield(stat_res,'stat') && isfield(stat_res,'pos') && ...
    isfield(stat_res,'dim') && isfield(stat_res,'inside'), ...
    'stat_res must contain stat, pos, dim and inside.');
assert(isfile(mri_ref_file),'Reference MRI does not exist.');
if ~exist(output_dir,'dir'), mkdir(output_dir); end

if exist('ft_defaults','file') ~= 2
    addpath(cfg.ft_path);
end
ft_defaults;

% Standardize source and atlas spatial units.
source_mm = ft_convert_units(stat_res,'mm');
atlas = ft_read_atlas(cfg.aal_file);
atlas_mm = ft_convert_units(atlas,'mm');

% Restrict AAL tissue Interpolate to source grid.
ci = [];
ci.interpmethod = 'nearest';
ci.parameter = 'tissue';
atlas_grid = ft_sourceinterpolate(ci,atlas_mm,source_mm);

roi_labels = get_roi_labels_local(atlas,cfg.roi_mode);
roi_mask = false(numel(source_mm.inside),1);
tissue = double(atlas_grid.tissue(:));

for k = 1:numel(roi_labels)
    tissue_id = find(strcmpi(atlas.tissuelabel,roi_labels{k}),1);
    assert(~isempty(tissue_id),'AAL label not found: %s',roi_labels{k});
    roi_mask = roi_mask | tissue == tissue_id;
end

roi_mask = roi_mask & logical(source_mm.inside(:));

F = double(source_mm.stat(:));
assert(numel(F)==numel(roi_mask),'stat and source grid sizes differ.');

keep = roi_mask & isfinite(F) & F > cfg.F_threshold;
assert(any(keep),'No source grid has F > %.4g inside the requested ROI.', ...
    cfg.F_threshold);

% Create a three-dimensional binary volume, and label all spatially connected components.
keep_vol = reshape(keep,source_mm.dim);
if exist('bwconncomp','file') == 2
    cc = bwconncomp(keep_vol,cfg.connectivity);
    label_vol = zeros(source_mm.dim,'uint16');
    for k = 1:cc.NumObjects
        label_vol(cc.PixelIdxList{k}) = k;
    end
else
    error(['bwconncomp is unavailable. Install MATLAB Image Processing ', ...
        'Toolbox or replace this block with an equivalent 3-D labeling method.']);
end

cluster_id = label_vol(:);
n_cluster = max(cluster_id);
cluster_size = zeros(n_cluster,1);
peak_F = nan(n_cluster,1);
peak_idx = nan(n_cluster,1);
peak_pos_mm = nan(n_cluster,3);

for k = 1:n_cluster
    idx = find(cluster_id==k);
    cluster_size(k) = numel(idx);
    [peak_F(k),j] = max(F(idx));
    peak_idx(k) = idx(j);
    peak_pos_mm(k,:) = source_mm.pos(peak_idx(k),:);
end

cluster_table = table((1:n_cluster)',cluster_size,peak_F, ...
    peak_pos_mm(:,1),peak_pos_mm(:,2),peak_pos_mm(:,3), ...
    'VariableNames',{'ClusterID','N_grid','Peak_F', ...
    'Peak_X_mm','Peak_Y_mm','Peak_Z_mm'});

% Create source parameter structure.
source_out = [];
source_out.pos = source_mm.pos;
source_out.dim = source_mm.dim;
source_out.inside = logical(source_mm.inside(:));
source_out.unit = 'mm';
if isfield(source_mm,'coordsys')
    source_out.coordsys = source_mm.coordsys;
end
% BrainNet For an entire background of NaN 's NIfTI suprathreshold voxels may not display.
source_out.F_gt_threshold = zeros(size(F));
source_out.F_gt_threshold(keep) = F(keep);
source_out.cluster_id = double(cluster_id);

% NIfTI spatial matrix should use mm to represent.Reference MRI May be cm(The old output was
% 0.5 cm/voxel written as 0.5 mm/voxel, Coordinates reduced tenfold).
mri_ref = ft_read_mri(mri_ref_file);
mri_ref = ft_convert_units(mri_ref,'mm');
ci = [];
ci.interpmethod = 'nearest';
ci.parameter = {'F_gt_threshold','cluster_id'};
mri_out = ft_sourceinterpolate(ci,source_out,mri_ref);

% Use zero background for exported maps; Formal statistics use only the original stat_res grid.
mri_out.F_gt_threshold(~isfinite(mri_out.F_gt_threshold)) = 0;
mri_out.cluster_id(~isfinite(mri_out.cluster_id)) = 0;
assert(nnz(mri_out.F_gt_threshold > cfg.F_threshold)>0 && ...
    nnz(mri_out.cluster_id>0)>0, ...
    'ROI voxels were lost during interpolation to the reference MRI.');
assert(strcmpi(mri_out.unit,'mm'), ...
    'Interpolated output is not in millimetres.');

f_file = fullfile(output_dir, ...
    sprintf('motor_%s_Fgt%g.nii',cfg.roi_mode,cfg.F_threshold));
c_file = fullfile(output_dir, ...
    sprintf('motor_%s_Fgt%g_clusterLabel.nii',cfg.roi_mode,cfg.F_threshold));

cw = [];
cw.parameter = 'F_gt_threshold';
cw.filename = f_file;
cw.filetype = 'nifti';
ft_sourcewrite(cw,mri_out);

% Verify the on-disk NIfTI, not just the pre-interpolation source grid.
cw.parameter = 'cluster_id';
cw.filename = c_file;
ft_sourcewrite(cw,mri_out);

% After writing both files, reread from disk for verification.
saved_F = ft_read_mri(f_file);
saved_labels = ft_read_mri(c_file);
assert(nnz(isfinite(saved_F.anatomy) & ...
    saved_F.anatomy > cfg.F_threshold)>0 && ...
    nnz(isfinite(saved_labels.anatomy) & ...
    saved_labels.anatomy > 0)>0, ...
    'Written NIfTI files have no surviving ROI voxels.');

out = struct();
out.F_nifti = f_file;
out.cluster_nifti = c_file;
out.cluster_table = cluster_table;
out.roi_mode = cfg.roi_mode;
out.roi_labels = roi_labels;
out.F_threshold = cfg.F_threshold;
out.n_roi_grid = nnz(roi_mask);
out.n_thresholded_grid = nnz(keep);
out.n_cluster = n_cluster;
out.source_mask = roi_mask;
out.threshold_mask = keep;
out.n_output_voxels = nnz(saved_F.anatomy > cfg.F_threshold);
out.output_transform_mm = saved_F.transform;

fprintf('\nExported motor ROI F>%.4g map.\n',cfg.F_threshold);
fprintf('ROI mode: %s (%s)\n',cfg.roi_mode,strjoin(roi_labels,', '));
fprintf('ROI grid points: %d\n',nnz(roi_mask));
fprintf('Retained grid points: %d\n',nnz(keep));
fprintf('Connected clusters: %d\n',n_cluster);
fprintf('Output NIfTI voxels above threshold: %d\n',out.n_output_voxels);
fprintf('Output NIfTI voxel spacing (mm): [%g %g %g]\n', ...
    vecnorm(saved_F.transform(1:3,1:3),2,1));
disp(cluster_table);
fprintf('F map: %s\n',f_file);
fprintf('Cluster labels: %s\n',c_file);
end

function labels = get_roi_labels_local(atlas,mode)
labels = {'Precentral_L','Precentral_R'};
if strcmpi(mode,'sensorimotor')
    labels = [labels,{'Postcentral_L','Postcentral_R'}];
elseif ~strcmpi(mode,'m1')
    error('roi_mode must be m1 or sensorimotor.');
end

for k = 1:numel(labels)
    idx = find(strcmpi(atlas.tissuelabel,labels{k}),1);
    assert(~isempty(idx),'AAL label not found: %s',labels{k});
    labels{k} = atlas.tissuelabel{idx};
end
end

function cfg = default_local(cfg,name,value)
if ~isfield(cfg,name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end

