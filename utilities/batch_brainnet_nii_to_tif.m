function result = batch_brainnet_nii_to_tif(user_cfg)
%BATCH_BRAINNET_NII_TO_TIF Batch-render NIfTI maps with BrainNet Viewer.
%
% Before the first batch run:
%   1. Open BrainNet Viewer once and load the intended surface and one NIfTI.
%   2. Set view, threshold, colormap, colorbar and image size.
%   3. Save the BrainNet options as a MAT file (it must contain EC).
%
% Example:
%   cfg = struct();
%   cfg.mapping_dir = 'F:\my_analysis\NIfTI_F_clusterSig';
%   cfg.cfg_file = 'F:\my_analysis\BrainNet_cfg.mat';
%   cfg.output_dir = fullfile(cfg.mapping_dir,'BrainNet_TIF');
%   cfg.input_pattern = '*.nii';
%   cfg.color_min = 0;
%   cfg.color_max = [];       % [] = global maximum across all NIfTIs
%   cfg.save_tif_only = true; % Save TIFFs without CSV/MAT sidecar files
%   cfg.output_suffix = '_medial_left'; % Append before .tif; default: ''
%   result = batch_brainnet_nii_to_tif(cfg);
%
% save_tif_only defaults to false, preserving CSV/MAT output. When true,
% metadata remains available in result, but is not saved to disk. Existing
% sidecar files in output_dir are not deleted. Temporary rendering settings
% are cleaned up automatically, including when rendering fails.

if nargin<1 || isempty(user_cfg)
    user_cfg = struct();
end
assert(isstruct(user_cfg),'user_cfg must be a structure.');

%% ------------------------------ Settings -----------------------------
cfg = struct();
codeRoot = fileparts(fileparts(mfilename('fullpath')));
cfg.brainnet_dir = fullfile(codeRoot,'utilities','BrainNet-Viewer');
% Optional fallback for installations without BrainNet.fig in the source.
cfg.brainnet_resource_dir = '';
cfg.spm_dir = ''; % Resolve the configured SPM installation from the path.
cfg.surface_file = fullfile(cfg.brainnet_dir,'BrainMesh_ICBM152.nv');

% These three paths should normally be changed by the user.
cfg.mapping_dir = fullfile(codeRoot,'figure_data','MEG_results');
cfg.output_dir = fullfile(codeRoot,'figure_outputs','brainnet');
cfg.cfg_file = ''; % Supply the intended saved BrainNet display configuration.

cfg.input_pattern = '*.nii';
cfg.overwrite = true;
cfg.save_tif_only = false;
cfg.output_suffix = '';

% BrainNet's GUIDE canvas (Tag='NV_axes') can remain visible underneath
% the rendered brain and produce ticks along the left/bottom image edges.
cfg.remove_axes = true;

% Positive source/F maps. All images should use one common scale.
cfg.fix_global_color_scale = true;
cfg.color_min = 0;
cfg.color_max = [];          % [] = scan all NIfTIs for global maximum

% If true, an all-zero NIfTI is skipped. If false, BrainNet is still called.
cfg.skip_all_zero = false;

cfg = overwrite_cfg_local(cfg,user_cfg);
if isempty(cfg.spm_dir)
    cfg.spm_dir = fileparts(which('spm'));
end
validateattributes(cfg.save_tif_only,{'logical','numeric'}, ...
    {'scalar','binary'},mfilename,'cfg.save_tif_only');
assert((ischar(cfg.output_suffix) && (isempty(cfg.output_suffix) || isrow(cfg.output_suffix))) || ...
    (isstring(cfg.output_suffix) && isscalar(cfg.output_suffix) && ~ismissing(cfg.output_suffix)), ...
    'cfg.output_suffix must be a character vector or scalar string.');
cfg.output_suffix = char(cfg.output_suffix);
assert(isempty(regexp(cfg.output_suffix,'[<>:"/\\|?*\x00-\x1F]','once')), ...
    'cfg.output_suffix must not contain path separators or invalid filename characters.');

%% -------------------------- Validate software ------------------------
assert(isfolder(cfg.brainnet_dir), ...
    'BrainNet directory does not exist: %s',cfg.brainnet_dir);
fig_in_source = fullfile(cfg.brainnet_dir,'BrainNet.fig');
fig_in_resource = fullfile(cfg.brainnet_resource_dir,'BrainNet.fig');
assert(isfile(fig_in_source) || (~isempty(cfg.brainnet_resource_dir) && isfile(fig_in_resource)), ...
    ['BrainNet.fig is missing from both the BrainNet source and fallback ' ...
     'resource directories.']);
assert(isfolder(cfg.spm_dir), ...
    ['Configure SPM12 first with setup_figure_paths(fieldtripRoot,spmRoot), ' ...
     'or supply cfg.spm_dir. Resolved directory: %s'],cfg.spm_dir);
assert(isfile(cfg.surface_file), ...
    'BrainNet surface file does not exist: %s',cfg.surface_file);
assert(isfolder(cfg.mapping_dir), ...
    'NIfTI directory does not exist: %s',cfg.mapping_dir);
assert(isfile(cfg.cfg_file), ...
    'BrainNet configuration file does not exist: %s',cfg.cfg_file);

% Path order is intentional: use readable BrainNet.m/BrainNet_MapCfg.m from
% brainnet_dir. Add the fallback only when BrainNet.fig is absent there.
if ~isfile(fig_in_source)
    addpath(cfg.brainnet_resource_dir,'-end');
end
addpath(cfg.brainnet_dir,'-begin');
addpath(cfg.spm_dir,'-begin');

assert(exist('BrainNet','file')==2, ...
    'BrainNet.m is not available on the MATLAB path.');
assert(exist('BrainNet_MapCfg','file')==2, ...
    'BrainNet_MapCfg.m is not available on the MATLAB path.');
assert(exist('BrainNet.fig','file')==2, ...
    ['BrainNet.fig is not available on the MATLAB path. Add its resource ' ...
     'directory before calling BrainNet_MapCfg.']);

fprintf('\nBrainNet          : %s\n',which('BrainNet'));
fprintf('BrainNet_MapCfg   : %s\n',which('BrainNet_MapCfg'));
fprintf('BrainNet.fig      : %s\n',which('BrainNet.fig'));
fprintf('Surface           : %s\n',cfg.surface_file);
fprintf('BrainNet settings : %s\n',cfg.cfg_file);

if ~isfolder(cfg.output_dir)
    mkdir(cfg.output_dir);
end

%% ------------------------- Find NIfTI files ---------------------------
D = dir(fullfile(cfg.mapping_dir,cfg.input_pattern));
D = D(~[D.isdir]);
assert(~isempty(D),'No NIfTI file matches %s.', ...
    fullfile(cfg.mapping_dir,cfg.input_pattern));

D = sort_nifti_files_local(D);
n_file = numel(D);

fprintf('Input NIfTIs     : %d\n',n_file);
fprintf('TIF output       : %s\n\n',cfg.output_dir);

%% ----------------------- Scan one common scale ------------------------
file_min = nan(n_file,1);
file_max = nan(n_file,1);
is_all_zero = false(n_file,1);

for i = 1:n_file
    nii_file = fullfile(D(i).folder,D(i).name);
    V = double(niftiread(nii_file));
    finite_value = V(isfinite(V));
    assert(~isempty(finite_value),'NIfTI has no finite voxel: %s',nii_file);
    file_min(i) = min(finite_value);
    file_max(i) = max(finite_value);
    is_all_zero(i) = all(finite_value==0);
end

if isempty(cfg.color_max)
    global_max = max(file_max,[],'omitnan');
else
    global_max = double(cfg.color_max);
end
assert(isscalar(global_max) && isfinite(global_max), ...
    'color_max/global maximum must be a finite scalar.');
if global_max<=cfg.color_min
    warning(['Global maximum %.6g is not larger than color_min %.6g. ' ...
        'Using color_min+1 for the BrainNet display range.'], ...
        global_max,cfg.color_min);
    global_max = cfg.color_min+1;
end

fprintf('Observed data range: [%.6g, %.6g]\n', ...
    min(file_min,[],'omitnan'),max(file_max,[],'omitnan'));
fprintf('BrainNet positive scale: [%.6g, %.6g]\n\n', ...
    cfg.color_min,global_max);

%% -------------------- Prepare fixed BrainNet config ------------------
render_cfg_file = cfg.cfg_file;
temporary_render_cfg = false;
if cfg.fix_global_color_scale
    cfg_content = load(cfg.cfg_file);
    assert(isfield(cfg_content,'EC'), ...
        'BrainNet configuration does not contain EC: %s',cfg.cfg_file);
    assert(isfield(cfg_content.EC,'vol'), ...
        'BrainNet EC structure does not contain EC.vol.');

    % BrainNet volume color limits:
    %   pn/px = positive minimum/maximum;
    %   nx/nn = negative minimum/maximum. The current workflow is positive.
    cfg_content.EC.vol.pn = cfg.color_min;
    cfg_content.EC.vol.px = global_max;
    cfg_content.EC.vol.nx = 0;
    cfg_content.EC.vol.nn = 0;

    if cfg.save_tif_only
        render_cfg_file = [tempname '.mat'];
        temporary_render_cfg = true;
        cleanup_render_cfg = onCleanup(@() delete_temporary_cfg_local(render_cfg_file)); %#ok<NASGU>
    else
        render_cfg_file = fullfile(cfg.output_dir, ...
            'BrainNet_cfg_fixed_global_scale.mat');
    end
    save(render_cfg_file,'-struct','cfg_content');
end

%% --------------------------- Render TIFs ------------------------------
tif_file = cell(n_file,1);
status = strings(n_file,1);

for i = 1:n_file
    nii_file = fullfile(D(i).folder,D(i).name);
    base_name = remove_nifti_extension_local(D(i).name);
    tif_file{i} = fullfile(cfg.output_dir,[base_name cfg.output_suffix '.tif']);

    fprintf('[%03d/%03d] %s\n',i,n_file,D(i).name);

    if is_all_zero(i) && cfg.skip_all_zero
        status(i) = "skipped_all_zero";
        fprintf('  skipped: all-zero map\n');
        continue;
    end
    if isfile(tif_file{i}) && ~cfg.overwrite
        status(i) = "kept_existing";
        fprintf('  kept existing: %s\n',tif_file{i});
        continue;
    end

    % Draw first without an output filename. BrainNet_MapCfg otherwise
    % prints internally before the residual GUIDE axes can be hidden.
    h = BrainNet_MapCfg(cfg.surface_file,nii_file,render_cfg_file);
    drawnow;

    if cfg.remove_axes
        remove_brainnet_frame_axes_local(h);
        drawnow;
    end

    print_brainnet_tif_local(h,tif_file{i},render_cfg_file);

    if ~isempty(h) && all(isgraphics(h))
        close(h);
    else
        close all;
    end

    assert(isfile(tif_file{i}), ...
        'BrainNet did not create the requested TIF: %s',tif_file{i});
    status(i) = "rendered";
end

%% ----------------------------- Manifest -------------------------------
input_file = arrayfun(@(x) fullfile(x.folder,x.name),D, ...
    'UniformOutput',false);

% dir() can return either a row or a column structure array depending on
% how the file list was produced.  table() interprets rows as observations,
% so force every manifest variable to the same n_file-by-1 shape.
input_file   = reshape(input_file,[],1);
tif_file     = reshape(tif_file,[],1);
file_min     = reshape(file_min,[],1);
file_max     = reshape(file_max,[],1);
is_all_zero  = reshape(is_all_zero,[],1);
status       = reshape(status,[],1);

manifest_height = [numel(input_file),numel(tif_file),numel(file_min), ...
    numel(file_max),numel(is_all_zero),numel(status)];
assert(all(manifest_height==n_file), ...
    'Manifest fields have inconsistent lengths: %s', ...
    mat2str(manifest_height));

manifest = table(input_file,tif_file,file_min,file_max,is_all_zero,status, ...
    'VariableNames',{'input_file','tif_file','file_min','file_max', ...
                     'is_all_zero','status'});
manifest_file = '';
if ~cfg.save_tif_only
    manifest_file = fullfile(cfg.output_dir,'BrainNet_TIF_manifest.csv');
    writetable(manifest,manifest_file);
end

result = struct();
result.config = cfg;
result.render_cfg_file = render_cfg_file;
if temporary_render_cfg
    % Do not return a path to a file deleted on function exit.
    result.render_cfg_file = '';
end
result.global_color_min = cfg.color_min;
result.global_color_max = global_max;
result.manifest = manifest;
result.manifest_file = manifest_file;
result.tif_files = tif_file;

if ~cfg.save_tif_only
    save(fullfile(cfg.output_dir,'BrainNet_batch_settings.mat'), ...
        'result','-v7.3');
end

fprintf('\nFinished. Rendered/available TIFs: %d/%d\n', ...
    sum(status=="rendered" | status=="kept_existing"),n_file);
if ~cfg.save_tif_only
    fprintf('Manifest: %s\n',manifest_file);
end
end


function delete_temporary_cfg_local(filename)
if isfile(filename)
    delete(filename);
end
end


function cfg = overwrite_cfg_local(cfg,user_cfg)
names = fieldnames(user_cfg);
for i = 1:numel(names)
    cfg.(names{i}) = user_cfg.(names{i});
end
end


function D = sort_nifti_files_local(D)
% Prefer numeric frame order when every filename contains frameXX.
frame_index = nan(numel(D),1);
for i = 1:numel(D)
    token = regexp(D(i).name,'frame[_-]?(\d+)', ...
        'tokens','once','ignorecase');
    if ~isempty(token)
        frame_index(i) = str2double(token{1});
    end
end
if all(isfinite(frame_index)) && numel(unique(frame_index))==numel(D)
    [~,order] = sort(frame_index);
else
    [~,order] = sort(lower(string({D.name})));
end
D = D(order);
end


function base_name = remove_nifti_extension_local(filename)
if endsWith(filename,'.nii.gz','IgnoreCase',true)
    base_name = extractBefore(filename,strlength(filename)-6);
    base_name = char(base_name);
else
    [~,base_name] = fileparts(filename);
end
end


function remove_brainnet_frame_axes_local(hfig)
% Hide only BrainNet's empty GUIDE canvas. Do not hide the actual brain
% axes or colorbar, both of which are needed in the exported image.
assert(~isempty(hfig) && all(isgraphics(hfig,'figure')), ...
    'BrainNet did not return a valid figure handle.');

frame_axes = findall(hfig,'Type','axes','Tag','NV_axes');
if isempty(frame_axes)
    warning(['BrainNet GUIDE axes (Tag=NV_axes) were not found. ' ...
        'No axes were removed to avoid hiding the brain or colorbar.']);
    return;
end

set(frame_axes, ...
    'Visible','off', ...
    'Box','off', ...
    'XTick',[], ...
    'YTick',[], ...
    'ZTick',[], ...
    'XColor','none', ...
    'YColor','none', ...
    'ZColor','none');
end


function print_brainnet_tif_local(hfig,tif_file,cfg_file)
% Reproduce BrainNet_MapCfg's TIFF export after the GUIDE frame is hidden.
S = load(cfg_file,'EC');
assert(isfield(S,'EC') && isfield(S.EC,'img'), ...
    'BrainNet configuration does not contain EC.img: %s',cfg_file);

required = {'width','height','dpi'};
for i = 1:numel(required)
    assert(isfield(S.EC.img,required{i}), ...
        'BrainNet configuration is missing EC.img.%s.',required{i});
end

width = double(S.EC.img.width);
height = double(S.EC.img.height);
dpi = double(S.EC.img.dpi);
assert(all(isfinite([width,height,dpi])) && all([width,height,dpi]>0), ...
    'BrainNet image width, height and dpi must be positive finite values.');

set(hfig,'PaperPositionMode','manual');
set(hfig,'PaperUnits','inch');
set(hfig,'PaperPosition',[1 1 width/dpi height/dpi]);
print(hfig,tif_file,'-dtiff',['-r',num2str(dpi)]);
end
