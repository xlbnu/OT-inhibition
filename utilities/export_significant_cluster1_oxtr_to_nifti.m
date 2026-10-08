function info = export_significant_cluster1_oxtr_to_nifti( ...
    stat_res, output_base, user_cfg)
%% Export significant positive cluster(s) x OXTR as a BrainNet NIfTI map
%
% Usage: 
%   % First significant positive cluster and OXTR intersection(Default)
%   info = export_significant_cluster1_oxtr_to_nifti( ...
%       stat_res, output_base, 'first');
%
%   % All significant positive clusters and OXTR intersection
%   info = export_significant_cluster1_oxtr_to_nifti( ...
%       stat_res, output_base, 'all');
%
%   % Use a structure to modify additional parameters
%   opt = struct();
%   opt.selection_mode = 'all';
%   opt.cluster_alpha = 0.05;
%   opt.oxtr_threshold = 0.5;
%   info = export_significant_cluster1_oxtr_to_nifti( ...
%       stat_res, output_base, opt);
%
% Output voxel values: 
%   selected significant cluster(s) AND OXTR -> Original stat_res.stat
%   Other voxels                                   -> 0
%
% If no significant cluster, or significant cluster and OXTR no intersection, this function still
% outputs complete spatial geometry, but all voxel values are 0 's NIfTI.BrainNet Loading this file
% will not display significant mapping, and batch plotting/animation processing will continue.

if nargin < 3 || isempty(user_cfg)
    user_cfg = struct();
elseif ischar(user_cfg) || (isstring(user_cfg) && isscalar(user_cfg))
    user_cfg = struct('selection_mode', char(user_cfg));
end
assert(isstruct(user_cfg), ...
    'The third input must be ''first'', ''all'', or a configuration struct.');

%% ============================== Configuration ===========================
cfg = struct();
cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg.selection_mode = 'first';       % 'first' or 'all'
cfg.cluster_alpha = 0.05;           % corrected cluster-level p threshold
cfg.oxtr_file = ...
    ['E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\' ...
     'MEG_data\MEG_mask\article\oxtr_mask\oxtr_50.nii'];
cfg.oxtr_parameter = 'anatomy';
cfg.oxtr_threshold = 0.5;
cfg.oxtr_interpmethod = 'nearest';  % binary mask Recommended nearest
cfg.template_mri_file = '';         % empty -> FieldTrip standard MRI

cfg = overwrite_cfg(cfg, user_cfg);
cfg.selection_mode = validatestring(lower(string(cfg.selection_mode)), ...
    {'first','all'}, mfilename, 'selection_mode');
cfg.selection_mode = char(cfg.selection_mode);

assert(isscalar(cfg.cluster_alpha) && isfinite(cfg.cluster_alpha) && ...
    cfg.cluster_alpha > 0 && cfg.cluster_alpha <= 1, ...
    'cluster_alpha must be in (0, 1].');
assert(isscalar(cfg.oxtr_threshold) && isfinite(cfg.oxtr_threshold), ...
    'oxtr_threshold must be a finite scalar.');

addpath(cfg.ft_path);
ft_defaults;

if isempty(cfg.template_mri_file)
    cfg.template_mri_file = fullfile(cfg.ft_path, 'template', 'anatomy', ...
        'single_subj_T1.nii');
end

assert(isstruct(stat_res) && isfield(stat_res,'stat') && ...
    isfield(stat_res,'posclusters') && ...
    isfield(stat_res,'posclusterslabelmat'), ...
    ['stat_res requires stat, posclusters, and ' ...
     'posclusterslabelmat.']);
assert(isfile(cfg.oxtr_file), 'Cannot find OXTR mask: %s', cfg.oxtr_file);
assert(isfile(cfg.template_mri_file), ...
    'Cannot find template MRI: %s', cfg.template_mri_file);

output_base = strip_nifti_extension(char(output_base));
output_dir = fileparts(output_base);
if ~isempty(output_dir) && ~exist(output_dir,'dir')
    mkdir(output_dir);
end

%% ================= Select significant positive cluster(s) ==============
cluster_p = double([stat_res.posclusters.prob]);
cluster_label = double(stat_res.posclusterslabelmat(:));

significant_cluster_ids = find(isfinite(cluster_p) & ...
    cluster_p <= cfg.cluster_alpha);

switch cfg.selection_mode
    case 'first'
        if isempty(significant_cluster_ids)
            selected_cluster_ids = [];
        else
            % "First"By FieldTrip posclusters defined by label order.
            selected_cluster_ids = significant_cluster_ids(1);
        end
    case 'all'
        selected_cluster_ids = significant_cluster_ids;
end

if isempty(selected_cluster_ids)
    selected_cluster_mask = false(size(cluster_label));
    warning(['No positive cluster survived cluster-level p <= %.4f. ' ...
        'An all-zero NIfTI will be written.'], cfg.cluster_alpha);
else
    selected_cluster_mask = ismember(cluster_label, selected_cluster_ids);
    missing_ids = selected_cluster_ids(~ismember( ...
        selected_cluster_ids, unique(cluster_label(cluster_label > 0))));
    if ~isempty(missing_ids)
        warning('Selected cluster label(s) absent from label matrix: %s', ...
            mat2str(missing_ids));
    end
end

fprintf('\nSelected positive cluster mode: %s\n', cfg.selection_mode);
fprintf('Cluster-level alpha          : %.4f\n', cfg.cluster_alpha);
if isempty(selected_cluster_ids)
    fprintf('Selected cluster IDs         : none\n');
else
    fprintf('Selected cluster IDs         : %s\n', ...
        mat2str(selected_cluster_ids));
    for i_cluster = selected_cluster_ids
        fprintf('  cluster %d: p=%.6f, grid points=%d\n', ...
            i_cluster, cluster_p(i_cluster), ...
            sum(cluster_label == i_cluster));
    end
end

%% ================= Interpolate OXTR mask to source grid =================
source_geometry = copy_source_geometry(stat_res);
oxtr_mri = ft_read_mri(cfg.oxtr_file);
cfg_oxtr = [];
cfg_oxtr.parameter = cfg.oxtr_parameter;
cfg_oxtr.interpmethod = cfg.oxtr_interpmethod;
oxtr_grid = ft_sourceinterpolate(cfg_oxtr, oxtr_mri, source_geometry);

assert(isfield(oxtr_grid, cfg.oxtr_parameter), ...
    'Interpolated OXTR structure lacks %s.', cfg.oxtr_parameter);

oxtr_value = double(oxtr_grid.(cfg.oxtr_parameter)(:));
assert(numel(oxtr_value) == numel(cluster_label), ...
    'OXTR/source-grid point counts differ.');
oxtr_mask = isfinite(oxtr_value) & oxtr_value > cfg.oxtr_threshold;

conjunction_mask = selected_cluster_mask & oxtr_mask;
has_intersection = any(conjunction_mask);

if ~has_intersection
    if isempty(selected_cluster_ids)
        warning('No selected significant cluster; writing an all-zero NIfTI.');
    else
        warning(['Selected cluster(s) %s have no intersection with the ' ...
            'OXTR mask; writing an all-zero NIfTI.'], ...
            mat2str(selected_cluster_ids));
    end
end

%% ========== Preserve F statistic in conjunction; zero elsewhere ========
stat_value = double(stat_res.stat(:));
assert(numel(stat_value) == numel(conjunction_mask), ...
    'stat and cluster-label point counts differ.');

export_value = zeros(size(stat_value));
if has_intersection
    export_value(conjunction_mask) = stat_value(conjunction_mask);
end
export_value(~isfinite(export_value)) = 0;

source_export = source_geometry;
source_export.stat = reshape(export_value, size(stat_res.stat));
source_export.statdimord = 'pos';

fprintf('\nSource-grid summary:\n');
fprintf('  all grid points                    : %d\n', numel(export_value));
fprintf('  selected significant cluster points: %d\n', ...
    sum(selected_cluster_mask));
fprintf('  OXTR grid points                   : %d\n', sum(oxtr_mask));
fprintf('  conjunction points retained       : %d\n', ...
    sum(conjunction_mask));
if has_intersection
    fprintf('  retained statistic range          : [%.4f, %.4f]\n', ...
        min(export_value(conjunction_mask)), ...
        max(export_value(conjunction_mask)));
else
    fprintf('  retained statistic range          : [0, 0]\n');
end

%% =================== Interpolate to standard MRI ========================
template_mri = ft_read_mri(cfg.template_mri_file);
template_mri = ft_convert_units(template_mri,'mm');

cfg_interp = [];
cfg_interp.parameter = 'stat';
cfg_interp.interpmethod = 'nearest';
stat_interp = ft_sourceinterpolate(cfg_interp, source_export, template_mri);

stat_interp.stat = double(stat_interp.stat);
stat_interp.stat(~isfinite(stat_interp.stat)) = 0;
if ~has_intersection
    % Force zeros again, Prevent nonfinite or boundary values when interpolating an all-zero source grid.
    stat_interp.stat(:) = 0;
end

%% ============================== Write NIfTI =============================
cfg_write = [];
cfg_write.filename = output_base;
cfg_write.filetype = 'nifti';
cfg_write.parameter = 'stat';
cfg_write.datatype = 'single';
cfg_write.scaling = 'no';
ft_volumewrite(cfg_write, stat_interp);

output_nii = [output_base '.nii'];
if ~isfile(output_nii) && isfile(output_base)
    output_nii = output_base;
end
assert(isfile(output_nii), 'NIfTI output was not created: %s', output_nii);

%% ============================== Return info =============================
info = struct();
info.selection_mode = cfg.selection_mode;
info.cluster_alpha = cfg.cluster_alpha;
info.significant_cluster_ids = significant_cluster_ids;
info.selected_cluster_ids = selected_cluster_ids;
info.selected_cluster_p = cluster_p(selected_cluster_ids);
info.n_selected_cluster_points = sum(selected_cluster_mask);
info.n_oxtr_points = sum(oxtr_mask);
info.n_conjunction_points = sum(conjunction_mask);
info.has_intersection = has_intersection;
info.is_blank_map = ~has_intersection;
info.output_nii = output_nii;

fprintf('\nNIfTI saved for BrainNet:\n  %s\n', output_nii);
if has_intersection
    fprintf(['Selected significant cluster(s) x OXTR retain the original ' ...
        'F statistic; all other voxels are zero.\n']);
else
    fprintf(['No selected-cluster x OXTR intersection was present; ' ...
        'the saved NIfTI contains zeros only.\n']);
end
end


function cfg = overwrite_cfg(cfg, user_cfg)
fields = fieldnames(user_cfg);
for i = 1:numel(fields)
    name = fields{i};
    assert(isfield(cfg, name), 'Unknown configuration field: %s', name);
    cfg.(name) = user_cfg.(name);
end
end


function source_geometry = copy_source_geometry(stat_res)
% Retain only spatial fields required for interpolation, Avoid FieldTrip attempting to interpret cluster statistical fields.
required_fields = {'pos','dim','inside'};
for i = 1:numel(required_fields)
    assert(isfield(stat_res, required_fields{i}), ...
        'stat_res lacks required geometry field: %s', required_fields{i});
end

source_geometry = struct();
geometry_fields = {'pos','dim','inside','unit','coordsys','transform'};
for i = 1:numel(geometry_fields)
    name = geometry_fields{i};
    if isfield(stat_res, name)
        source_geometry.(name) = stat_res.(name);
    end
end
end


function output_base = strip_nifti_extension(output_base)
% FieldTrip expects a base filename and appends .nii itself.
output_base = regexprep(char(output_base), '(?i)\.nii(\.gz)?$', '');
end
