function result = export_cluster1_from_stat_file( ...
    stat_file,output_dir,user_cfg)
%EXPORT_CLUSTER1_FROM_STAT_FILE
% From a single MAT statistical file, read stat_res, Output the first significant positive cluster 's NIfTI.
%
% Input: 
%   stat_file  : Included variable stat_res 's MAT file
%   output_dir : NIfTI Output folder
%   user_cfg   : Optional parameter structure
%       .ft_path            FieldTrip path
%       .alpha              cluster-level significance threshold, Default 0.05
%       .overwrite          Whether to overwrite existing files, Default true
%       .template_mri_file  Standard template used for interpolation MRI
%
% Output filename: 
%   <Input MAT filename>_cluster1.nii
%
% If no significant cluster, still create voxels all equal to 0 's NIfTI, Meanwhile: 
%   result.found    = false
%   result.is_blank = true

% Example: 
%   result = export_cluster1_from_stat_file( ...
%       'F:\stats\frame_01_stat.mat', ...
%       'F:\stats\cluster1_nifti');

if nargin < 3 || isempty(user_cfg)
    user_cfg = struct();
end

assert(isfile(stat_file),'Statistics file does not exist: %s',stat_file);
assert(ischar(output_dir) || isstring(output_dir), ...
    'output_dir must be a character vector or string scalar.');
output_dir = char(output_dir);

cfg = struct();
cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg.alpha = 0.05;
cfg.overwrite = true;

names = fieldnames(user_cfg);
for i = 1:numel(names)
    cfg.(names{i}) = user_cfg.(names{i});
end

if ~isfield(user_cfg,'template_mri_file') || ...
        isempty(user_cfg.template_mri_file)
    cfg.template_mri_file = fullfile(cfg.ft_path,'template','anatomy', ...
        'single_subj_T1_1mm.nii');
end

assert(isfolder(cfg.ft_path),'FieldTrip path does not exist: %s',cfg.ft_path);
assert(isfile(cfg.template_mri_file), ...
    'Template MRI does not exist: %s',cfg.template_mri_file);
assert(isnumeric(cfg.alpha) && isscalar(cfg.alpha) && ...
    isfinite(cfg.alpha) && cfg.alpha>0 && cfg.alpha<1, ...
    'cfg.alpha must be between 0 and 1.');

addpath(cfg.ft_path);
ft_defaults;

loaded = load(stat_file,'stat_res');
assert(isfield(loaded,'stat_res'), ...
    'Input MAT file does not contain the variable stat_res: %s',stat_file);
stat_res = loaded.stat_res;

assert(isfield(stat_res,'stat') && ~isempty(stat_res.stat), ...
    'stat_res lacks stat values.');
assert(isfield(stat_res,'pos') && ...
    size(stat_res.pos,1)==numel(stat_res.stat), ...
    'stat_res.pos does not match stat_res.stat.');
assert(isfield(stat_res,'dim'),'stat_res lacks dim.');

[~,input_name] = fileparts(stat_file);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
filename_no_ext = fullfile(output_dir,[input_name '_cluster1']);
output_file = [filename_no_ext '.nii'];

result = struct();
result.input_file = stat_file;
result.output_file = output_file;
result.found = false;
result.is_blank = true;
result.cluster_id = NaN;
result.cluster_p = NaN;
result.n_source_points = 0;
result.peak_F = NaN;
result.peak_source_index = NaN;
result.peak_pos = [NaN NaN NaN];

cluster_mask = false(numel(stat_res.stat),1);
cluster_id = NaN;
cluster_p_value = NaN;

if isfield(stat_res,'posclusters') && ~isempty(stat_res.posclusters) && ...
        isfield(stat_res,'posclusterslabelmat') && ...
        ~isempty(stat_res.posclusterslabelmat)
    cluster_p = double([stat_res.posclusters.prob]);
    significant_ids = find(isfinite(cluster_p) & cluster_p<=cfg.alpha);
    if ~isempty(significant_ids)
        cluster_id = significant_ids(1);
        cluster_p_value = cluster_p(cluster_id);
        cluster_mask = ...
            double(stat_res.posclusterslabelmat(:))==cluster_id;
        cluster_mask = cluster_mask & isfinite(double(stat_res.stat(:)));
    end
end

template_mri = ft_read_mri(cfg.template_mri_file);

if any(cluster_mask)
    result.found = true;
    result.is_blank = false;
    result.cluster_id = cluster_id;
    result.cluster_p = cluster_p_value;
    result.n_source_points = nnz(cluster_mask);

    if isfield(stat_res,'unit') && ~isempty(stat_res.unit)
        template_mri = ft_convert_units(template_mri,stat_res.unit);
    end

    F_source = zeros(numel(stat_res.stat),1);
    F_source(cluster_mask) = double(stat_res.stat(cluster_mask));

    source_F = make_source_geometry_local(stat_res);
    source_F.cluster1_F = F_source;

    source_mask = make_source_geometry_local(stat_res);
    source_mask.cluster1_mask = double(cluster_mask);

    cfg_interp = [];
    cfg_interp.parameter = 'cluster1_F';
    cfg_interp.interpmethod = 'linear';
    interp_F = ft_sourceinterpolate(cfg_interp,source_F,template_mri);

    cfg_interp = [];
    cfg_interp.parameter = 'cluster1_mask';
    cfg_interp.interpmethod = 'nearest';
    interp_mask = ft_sourceinterpolate(cfg_interp,source_mask,template_mri);

    output_data = double(interp_F.cluster1_F);
    mask_volume = double(interp_mask.cluster1_mask)>=0.5;
    output_data(~isfinite(output_data) | ~mask_volume) = 0;

    F_peak = double(stat_res.stat(:));
    F_peak(~cluster_mask) = -Inf;
    [result.peak_F,result.peak_source_index] = max(F_peak);
    result.peak_pos = double( ...
        stat_res.pos(result.peak_source_index,:));
else
    % No significant cluster when, Output matches the template MRI in size and is entirely 0 NIfTI.
    output_data = zeros(double(template_mri.dim(:).'),'single');
end

if isfile(output_file) && ~cfg.overwrite
    fprintf('Output exists and overwrite=false: %s\n',output_file);
    return;
end

volume = [];
volume.dim = template_mri.dim;
volume.transform = template_mri.transform;
volume.unit = template_mri.unit;
if isfield(template_mri,'coordsys')
    volume.coordsys = template_mri.coordsys;
end
volume.cluster1_F = output_data;

cfg_write = [];
cfg_write.filename = filename_no_ext;
cfg_write.filetype = 'nifti';
cfg_write.parameter = 'cluster1_F';
cfg_write.datatype = 'single';
cfg_write.scaling = 'no';
ft_volumewrite(cfg_write,volume);
assert(isfile(output_file),'NIfTI creation failed: %s',output_file);

if result.found
    fprintf(['Exported cluster ID %d (cluster p=%.6g, %d source points):\n' ...
        '%s\n'],result.cluster_id,result.cluster_p, ...
        result.n_source_points,output_file);
else
    fprintf(['No cluster with cluster-level p<=%.4g was found. ' ...
        'A zero-valued NIfTI was created:\n%s\n'], ...
        cfg.alpha,output_file);
end
end


function source = make_source_geometry_local(stat_res)
source = [];
source.pos = stat_res.pos;
source.dim = stat_res.dim;
if isfield(stat_res,'inside') && ~isempty(stat_res.inside)
    source.inside = stat_res.inside;
else
    source.inside = true(size(stat_res.pos,1),1);
end
if isfield(stat_res,'unit')
    source.unit = stat_res.unit;
end
if isfield(stat_res,'coordsys')
    source.coordsys = stat_res.coordsys;
end
end
