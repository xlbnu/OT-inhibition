function result = brainnet_sliding_window_gif(analysisName,cfg)
% BRAINNET_SLIDING_WINDOW_GIF Render S7 source NIfTIs directly to one GIF.
%   setup_figure_paths([], '/path/to/spm12')
%   result = brainnet_sliding_window_gif('decision_gamma_band');
%   Supported names: decision_gamma_band, decision_beta_band,
%   motor_gamma_band, motor_beta_band. Only the final GIF is retained.
%   Each group's source/BrainNet_cfg.mat supplies its view and colormap.
%   Default positive limits: [5 15], [5 25], [1 5], [5 15], respectively.
%   cfg.zlim overrides limits; cfg.surface_file selects another surface.
%   Default surfaces: decision gamma left, decision beta right, motors whole.
%   cfg.frame_rate defaults to 1. Header text uses window midpoints and
%   frequency ranges (gamma [30 100], beta [13 30] Hz).
%   cfg.image_size is [width height] in pixels (default [700 700]).
%   Fonts follow the scalp GIF: 8 pt at final_height_cm=5. Optional
%   final_width_cm overrides height calibration. crop_fraction is
%   [left top width height] relative to the complete BrainNet rendering.
%   Set save_frames=true to retain labelled PNGs in frames_dir. Set
%   write_gif=false to export frames before a separate GIF assembly step.
if nargin<2, cfg = struct(); end
analysisName = char(analysisName);
names = {'decision_gamma_band','decision_beta_band','motor_gamma_band','motor_beta_band'};
index = find(strcmp(analysisName,names));
assert(~isempty(index),'Unknown analysis name.');
codeRoot = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(codeRoot,'figure_data','MEG_results','source_dynamic_gif',analysisName,'source');
limits = [5 15;5 25;1 5;5 15];
defaults = struct('input_dir',inputDir,'cfg_file','', ...
    'output_file',fullfile(codeRoot,'figure_outputs','figureS7',[analysisName '_source.gif']), ...
    'brainnet_dir',fullfile(codeRoot,'utilities','BrainNet-Viewer'),'spm_dir','', ...
    'surface_file','','zlim',limits(index,:),'frame_rate',1, ...
    'image_size',[700 700],'render_width_px',1000, ...
    'crop_fraction',[100/2000 0 1670/2000 1400/1500], ...
    'final_height_cm',5,'final_width_cm',[],'final_font_size_pt',8, ...
    'line_spacing_ratio',0.95,'text_top_px',3,'text_bottom_margin_px',2, ...
    'save_frames',false,'frames_dir','','write_gif',true);
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(cfg,fields{k}), cfg.(fields{k}) = defaults.(fields{k}); end
end
if isempty(cfg.cfg_file), cfg.cfg_file = fullfile(cfg.input_dir,'BrainNet_cfg.mat'); end
if isempty(cfg.surface_file)
    surfaces = {'BrainMesh_ICBM152Left.nv','BrainMesh_ICBM152Right.nv', ...
        'BrainMesh_ICBM152.nv','BrainMesh_ICBM152.nv'};
    cfg.surface_file = fullfile(cfg.brainnet_dir,surfaces{index});
end
if isempty(cfg.spm_dir), cfg.spm_dir = fileparts(which('spm')); end
assert(isfile(fullfile(cfg.spm_dir,'spm.m')),'Configure SPM12 or supply cfg.spm_dir.');
assert(isfile(fullfile(cfg.brainnet_dir,'BrainNet_MapCfg.m')) && ...
    isfile(fullfile(cfg.brainnet_dir,'BrainNet.fig')),'BrainNet files are missing.');
assert(isfile(cfg.surface_file) && isfile(cfg.cfg_file),'Surface or BrainNet configuration is missing.');
addpath(cfg.brainnet_dir,'-begin'); addpath(cfg.spm_dir,'-begin');
validateattributes(cfg.zlim,{'numeric'},{'vector','numel',2,'finite','nonnegative'});
assert(cfg.zlim(2)>cfg.zlim(1),'zlim must be increasing.');
validateattributes(cfg.frame_rate,{'numeric'},{'scalar','finite','positive','<=',100});
validateattributes(cfg.image_size,{'numeric'},{'vector','numel',2,'integer','>=',200});
validateattributes(cfg.render_width_px,{'numeric'},{'scalar','integer','>=',200});
validateattributes(cfg.crop_fraction,{'numeric'},{'vector','numel',4,'finite','nonnegative'});
assert(all(cfg.crop_fraction(3:4)>0) && all(cfg.crop_fraction(1:2)+cfg.crop_fraction(3:4)<=1), ...
    'crop_fraction must describe a rectangle inside the rendered image.');
for field = {'final_height_cm','final_font_size_pt','line_spacing_ratio'}
    validateattributes(cfg.(field{1}),{'numeric'},{'scalar','finite','positive'});
end
for field = {'text_top_px','text_bottom_margin_px'}
    validateattributes(cfg.(field{1}),{'numeric'},{'scalar','finite','nonnegative'});
end
if ~isempty(cfg.final_width_cm)
    validateattributes(cfg.final_width_cm,{'numeric'},{'scalar','finite','positive'});
end
[outputDir,~,ext] = fileparts(cfg.output_file);
assert(strcmpi(ext,'.gif'),'output_file must end in .gif.');
validateattributes(cfg.save_frames,{'logical','numeric'},{'scalar','binary'});
validateattributes(cfg.write_gif,{'logical','numeric'},{'scalar','binary'});
assert(cfg.save_frames || cfg.write_gif,'Enable save_frames or write_gif.');
if isempty(cfg.frames_dir)
    cfg.frames_dir = fullfile(outputDir,'frames',analysisName,'source');
end
if cfg.save_frames && ~isfolder(cfg.frames_dir), mkdir(cfg.frames_dir); end
files = dir(fullfile(cfg.input_dir,'frame_*.nii'));
assert(~isempty(files),'No frame NIfTIs in %s.',cfg.input_dir);
windows = zeros(numel(files),2);
for k = 1:numel(files)
    tokens = regexp(files(k).name,'^frame_\d+_t_([mp]?\d+p\d+)_([mp]?\d+p\d+)_.*\.nii$', ...
        'tokens','once');
    assert(~isempty(tokens),'Unrecognized time-window filename: %s',files(k).name);
    for j = 1:2
        token = tokens{j};
        if startsWith(token,'p'), token = token(2:end); end
        windows(k,j) = str2double(strrep(strrep(token,'m','-'),'p','.'));
    end
end
[~,order] = sort(windows(:,1)); windows = windows(order,:); files = files(order);
assert(all(diff(windows(:,1))>0) && all(windows(:,2)>windows(:,1)), ...
    'Time windows must be unique and ordered with positive durations.');
if contains(analysisName,'gamma'), frequency = [30 100]; else, frequency = [13 30]; end
% Modify a temporary configuration only; original MAT/NIfTI files are untouched.
settings = load(cfg.cfg_file);
assert(isfield(settings,'EC') && isfield(settings.EC,'vol') && isfield(settings.EC,'img'), ...
    'BrainNet configuration must contain EC.vol and EC.img.');
settings.EC.vol.pn = cfg.zlim(1); settings.EC.vol.px = cfg.zlim(2);
settings.EC.vol.display = 2;
settings.EC.vol.nn = 0; settings.EC.vol.nx = 0;
renderHeight = round(cfg.render_width_px*settings.EC.img.height/settings.EC.img.width);
temporaryCfg = [tempname '.mat'];
save(temporaryCfg,'-struct','settings');
tempCleanup = onCleanup(@() delete(temporaryCfg));
fig = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 cfg.image_size], ...
    'InvertHardcopy','off');
figureCleanup = onCleanup(@() close_if_valid(fig));
dpi = 100;
paperSize = cfg.image_size/dpi;
set(fig,'PaperUnits','inches','PaperPosition',[0 0 paperSize], ...
    'PaperSize',paperSize,'PaperPositionMode','manual');
if isempty(cfg.final_width_cm)
    fontPoints = cfg.final_font_size_pt*paperSize(2)/(cfg.final_height_cm/2.54);
else
    fontPoints = cfg.final_font_size_pt*paperSize(1)/(cfg.final_width_cm/2.54);
end
fontPixels = fontPoints*dpi/72; lineHeight = 1.25*fontPixels;
lineSpacing = cfg.line_spacing_ratio*fontPixels;
headerHeight = ceil(cfg.text_top_px+lineSpacing+lineHeight+cfg.text_bottom_margin_px);
assert(headerHeight<0.5*cfg.image_size(2),'Header is too large for the requested image.');
frames = cell(numel(files),1);
frameFiles = cell(numel(files),1);
for k = 1:numel(files)
    h = BrainNet_MapCfg(cfg.surface_file,fullfile(files(k).folder,files(k).name),temporaryCfg);
    renderCleanup = onCleanup(@() close_if_valid(h));
    set(h,'Visible','off');
    canvas = findall(h,'Type','axes','Tag','NV_axes');
    set(canvas,'Visible','off','Box','off','XTick',[],'YTick',[],'ZTick',[]);
    set(h,'PaperUnits','inches','PaperPosition',[0 0 cfg.render_width_px/dpi renderHeight/dpi], ...
        'PaperPositionMode','manual','InvertHardcopy','off');
    drawnow;
    rgb = print(h,'-RGBImage','-r100');
    clear renderCleanup;
    left = floor(cfg.crop_fraction(1)*size(rgb,2))+1;
    top = floor(cfg.crop_fraction(2)*size(rgb,1))+1;
    right = min(size(rgb,2),floor(sum(cfg.crop_fraction([1 3]))*size(rgb,2)));
    bottom = min(size(rgb,1),floor(sum(cfg.crop_fraction([2 4]))*size(rgb,1)));
    rgb = rgb(top:bottom,left:right,:);
    % Fit the cropped rendering below the header without changing aspect ratio.
    available = [cfg.image_size(1),cfg.image_size(2)-headerHeight];
    scale = min(available./[size(rgb,2) size(rgb,1)]);
    bodySize = scale*[size(rgb,2) size(rgb,1)];
    bodyLeft = (cfg.image_size(1)-bodySize(1))/2;
    bodyBottom = cfg.image_size(2)-headerHeight-bodySize(2);
    clf(fig);
    ax = axes('Parent',fig,'Units','pixels','Position',[bodyLeft bodyBottom bodySize]);
    image(ax,rgb); axis(ax,'image'); axis(ax,'off');
    center = mean(windows(k,:)); if abs(center)<1e-10, center = 0; end
    timeText = regexprep(sprintf('%.2f',center),'\.?0+$','');
    if isempty(timeText), timeText = '0'; end
    lines = {sprintf('Time: %ss',timeText),sprintf('Frequency window: [%g %g]Hz',frequency)};
    for line = 1:2
        y = cfg.text_top_px+(line-1)*lineSpacing;
        annotation(fig,'textbox',[0.01 1-(y+lineHeight)/cfg.image_size(2) ...
            0.98 lineHeight/cfg.image_size(2)],'String',lines{line}, ...
            'EdgeColor','none','Margin',0,'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontUnits','points','FontSize',fontPoints,'FontName','Arial','Interpreter','none');
    end
    drawnow; frames{k} = print(fig,'-RGBImage','-r100');
    if cfg.save_frames
        [~,base] = fileparts(files(k).name);
        frameFiles{k} = fullfile(cfg.frames_dir,[base '.png']);
        imwrite(frames{k},frameFiles{k});
    end
    fprintf('%s source: frame %d/%d, [%.2f, %.2f] s\n',analysisName,k,numel(files),windows(k,:));
end
if cfg.write_gif
samples = cellfun(@(rgb) rgb(round(linspace(1,size(rgb,1),128)), ...
    round(linspace(1,size(rgb,2),128)),:),frames,'UniformOutput',false);
[~,palette] = rgb2ind(cat(1,samples{:}),256,'nodither');
if ~isempty(outputDir) && ~isfolder(outputDir), mkdir(outputDir); end
for k = 1:numel(frames)
    indexed = rgb2ind(frames{k},palette,'nodither');
    if k==1
        imwrite(indexed,palette,cfg.output_file,'gif','LoopCount',Inf,'DelayTime',1/cfg.frame_rate);
    else
        imwrite(indexed,palette,cfg.output_file,'gif','WriteMode','append','DelayTime',1/cfg.frame_rate);
    end
end
end
result = struct('output_file',cfg.output_file,'time_windows',windows, ...
    'time_centers',mean(windows,2),'frame_count',numel(files),'zlim',cfg.zlim, ...
    'frequency_window_hz',frequency,'cfg_file',cfg.cfg_file,'surface_file',cfg.surface_file, ...
    'frame_files',{frameFiles},'gif_written',logical(cfg.write_gif));
end

function close_if_valid(h)
if isgraphics(h,'figure'), close(h); end
end
