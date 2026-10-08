function result = topoplot_sliding_window_gif(analysisName, cfg)
% TOPOPLOT_SLIDING_WINDOW_GIF Render S7 scalp maps directly to a GIF.
%   result = topoplot_sliding_window_gif('decision_gamma_band');
%   Supported names: decision_gamma_band, decision_beta_band,
%   motor_gamma_band, motor_beta_band. Configure FieldTrip first.
%   Optional cfg fields: input_dir, output_file, frame_rate (default 1),
%   zlim (analysis-specific defaults below), image_size ([700 750], width/height).
%   Default t-value limits: decision gamma [-3 3], decision beta [-7 7],
%   motor gamma [-6 6], motor beta [-3 6]. cfg.zlim overrides these limits.
%   Labels follow the original crop/annotation proportions: 8 pt at a final
%   height of 5 cm, with compact line spacing and a small gap above the map.
%   Optional final_width_cm overrides height-based font calibration.
%   Optional top_padding_px, text_top_px, text_bottom_margin_px and
%   line_spacing_ratio adjust the header; sizes refer to exported pixels.
%   Default: only the GIF is saved. Set save_frames=true to retain labelled
%   PNGs in frames_dir. Set write_gif=false for frame export without a GIF.
%   Frames are sorted by window start time parsed
%   from their filenames, and checked against stored window metadata.
%   Motor inputs use top-level ot_freq_right/ot_freq_left/pl_freq_right/pl_freq_left.
%   Optional t_win metadata is checked against the filename window.
%   Maps reproduce the supplied script: averaged OT/PL log contrasts,
%   followed by a participant-level one-sample t-test. Uncorrected P<0.05
%   controls opacity only; this is not a cluster-permutation result.
if nargin < 2, cfg = struct(); end
analysisName = char(analysisName);
names = {'decision_gamma_band','decision_beta_band', ...
    'motor_gamma_band','motor_beta_band'};
assert(ismember(analysisName,names),'Unknown analysis name.');
assert(exist('ft_topoplotER','file') == 2,'Initialize FieldTrip first.');
assert(exist('ttest','file') == 2,'Statistics and Machine Learning Toolbox is required.');
codeRoot = fileparts(fileparts(mfilename('fullpath')));
defaultLimits = [-3 3; -7 7; -6 6; -3 6];
analysisIndex = find(strcmp(analysisName,names));
defaults = struct('input_dir',fullfile(codeRoot,'figure_data','MEG_results', ...
    'source_dynamic_gif',analysisName,'topoplot'), ...
    'output_file',fullfile(codeRoot,'figure_outputs','figureS7', ...
    [analysisName '_topoplot.gif']),'frame_rate',1,'zlim',defaultLimits(analysisIndex,:), ...
    'image_size',[700 750],'final_height_cm',5,'final_width_cm',[], ...
    'final_font_size_pt',8,'line_spacing_ratio',0.95, ...
    'save_frames',false,'frames_dir','','write_gif',true);
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(cfg,fields{k}), cfg.(fields{k}) = defaults.(fields{k}); end
end
validateattributes(cfg.frame_rate,{'numeric'},{'scalar','finite','positive','<=',100});
validateattributes(cfg.zlim,{'numeric'},{'vector','numel',2,'finite'});
assert(cfg.zlim(2)>cfg.zlim(1),'zlim must be increasing.');
validateattributes(cfg.image_size,{'numeric'},{'vector','numel',2,'integer','>=',200});
validateattributes(cfg.final_height_cm,{'numeric'},{'scalar','finite','positive'});
if ~isempty(cfg.final_width_cm)
    validateattributes(cfg.final_width_cm,{'numeric'},{'scalar','finite','positive'});
end
validateattributes(cfg.final_font_size_pt,{'numeric'},{'scalar','finite','positive'});
validateattributes(cfg.line_spacing_ratio,{'numeric'},{'scalar','finite','positive'});
% Match the original 4200 x 4500 crop, including its 500-pixel top padding.
headerDefaults = struct('top_padding_px',round(cfg.image_size(2)*500/4500), ...
    'text_top_px',round(cfg.image_size(2)*18/4500), ...
    'text_bottom_margin_px',max(1,round(cfg.image_size(2)*8/4500)));
fields = fieldnames(headerDefaults);
for k = 1:numel(fields)
    if ~isfield(cfg,fields{k}), cfg.(fields{k}) = headerDefaults.(fields{k}); end
    validateattributes(cfg.(fields{k}),{'numeric'},{'scalar','finite','nonnegative'});
end
[outputDir,~,extension] = fileparts(cfg.output_file);
assert(strcmpi(extension,'.gif'),'output_file must end in .gif.');
validateattributes(cfg.save_frames,{'logical','numeric'},{'scalar','binary'});
validateattributes(cfg.write_gif,{'logical','numeric'},{'scalar','binary'});
assert(cfg.save_frames || cfg.write_gif,'Enable save_frames or write_gif.');
if isempty(cfg.frames_dir)
    cfg.frames_dir = fullfile(outputDir,'frames',analysisName,'topoplot');
end
if cfg.save_frames && ~isfolder(cfg.frames_dir), mkdir(cfg.frames_dir); end
files = dir(fullfile(cfg.input_dir,'frame_*_topo.mat'));
assert(~isempty(files),'No scalp-map MAT files found in %s.',cfg.input_dir);
windows = zeros(numel(files),2);
for k = 1:numel(files)
    token = regexp(files(k).name, ...
        '^frame_\d+_t_(m?\d+p\d+)_(m?\d+p\d+)_topo\.mat$', 'tokens','once');
    assert(~isempty(token),'Unrecognized window filename: %s',files(k).name);
    windows(k,:) = cellfun(@(s) str2double(strrep(strrep(s,'m','-'),'p','.')),token);
end
[~,order] = sort(windows(:,1));
windows = windows(order,:); files = files(order);
assert(all(diff(windows(:,1))>0) && all(windows(:,2)>windows(:,1)), ...
    'Window starts must be unique and window ends must follow starts.');
if contains(analysisName,'gamma'), band = '[60 80]Hz'; else, band = '[13 30]Hz'; end
fig = figure('Visible','off','Color','w','Units','pixels', ...
    'Position',[100 100 cfg.image_size],'InvertHardcopy','off');
cleanup = onCleanup(@() close(fig));
% Fix export dimensions independently of display DPI. At the target PPT width,
% the embedded text has the requested physical point size after scaling.
exportDpi = 100;
paperSize = cfg.image_size/exportDpi;
set(fig,'PaperUnits','inches','PaperPosition',[0 0 paperSize], ...
    'PaperSize',paperSize,'PaperPositionMode','manual');
if isempty(cfg.final_width_cm)
    labelFontSize = cfg.final_font_size_pt*paperSize(2)/(cfg.final_height_cm/2.54);
else
    labelFontSize = cfg.final_font_size_pt*paperSize(1)/(cfg.final_width_cm/2.54);
end
fontPixels = labelFontSize*exportDpi/72;
lineSpacing = cfg.line_spacing_ratio*fontPixels;
lineHeight = 1.25*fontPixels;
headerHeight = max(cfg.top_padding_px,ceil(cfg.text_top_px+lineSpacing+ ...
    lineHeight+cfg.text_bottom_margin_px));
bottomMargin = 0.015*cfg.image_size(2);
assert(headerHeight+bottomMargin<0.5*cfg.image_size(2), ...
    'Header is too large; reduce font size, spacing, or padding.');
frames = cell(numel(files),1);
frameFiles = cell(numel(files),1);
layout = [];
for k = 1:numel(files)
    path = fullfile(files(k).folder,files(k).name);
    variables = whos('-file',path);
    if any(strcmp({variables.name},'t_win'))
        metadata = load(path,'t_win');
        assert(numel(metadata.t_win)==2 && ...
            max(abs(metadata.t_win(:)'-windows(k,:)))<1e-6, ...
            'Stored time window mismatch.');
    end
    if startsWith(analysisName,'decision')
        s = load(path,'ot_freq_conf','ot_freq_cong','pl_freq_conf','pl_freq_cong');
        groups = {s.ot_freq_conf,s.ot_freq_cong,s.pl_freq_conf,s.pl_freq_cong};
    else
        s = load(path,'ot_freq_right','ot_freq_left','pl_freq_right','pl_freq_left');
        groups = {s.ot_freq_right,s.ot_freq_left,s.pl_freq_right,s.pl_freq_left};
    end
    reference = groups{1}{1};
    nSubjects = numel(groups{1}); nChannels = numel(reference.label);
    values = zeros(nSubjects,nChannels,4);
    for g = 1:4
        assert(numel(groups{g})==nSubjects,'Participant counts differ between conditions.');
        for subject = 1:nSubjects
            item = groups{g}{subject};
            [found,channelOrder] = ismember(reference.label,item.label);
            assert(all(found) && numel(unique(item.label))==nChannels && ...
                isequal(item.freq,reference.freq), ...
                'Channel sets or frequency bins differ between inputs.');
            power = item.powspctrm;
            assert(numel(power)==nChannels && all(isfinite(power(:))) && all(power(:)>0), ...
                'Expected one finite, positive power value per channel.');
            % Align by channel name: some supplied participants use a different order.
            values(subject,:,g) = log(power(channelOrder))';
        end
    end
    contrast = ((values(:,:,1)-values(:,:,2))+(values(:,:,3)-values(:,:,4)))/2;
    [significant,~,~,stat] = ttest(contrast,0,'Dim',1,'Alpha',0.05);
    data = struct('label',{reference.label},'time',0,'dimord','chan_time', ...
        'avg',stat.tstat(:),'sig_mask',0.18+0.82*double(significant(:)));
    assert(all(isfinite(data.avg)),'Nonfinite t-statistic in %s.',files(k).name);
    if isempty(layout)
        layoutCfg = struct('layout','neuromag306cmb.lay','channel',{data.label});
        layout = ft_prepare_layout(layoutCfg,data);
        labels = data.label;
    else
        assert(isequal(labels,data.label),'Channel labels change across frames.');
    end
    clf(fig);
    ax = axes('Parent',fig,'Position',[0.05 bottomMargin/cfg.image_size(2) 0.90 ...
        (cfg.image_size(2)-headerHeight-bottomMargin)/cfg.image_size(2)]);
    plotCfg = struct('layout',layout,'parameter','avg','zlim',cfg.zlim, ...
        'marker','off','comment','no','colorbar','no','interactive','no', ...
        'maskparameter','sig_mask','figure',ax);
    ft_topoplotER(plotCfg,data);
    colormap(fig,adaptive_cmap(cfg.zlim));
    centerTime = mean(windows(k,:));
    if abs(centerTime)<1e-10, centerTime = 0; end
    timeText = regexprep(sprintf('%.2f',centerTime),'\.?0+$','');
    if isempty(timeText), timeText = '0'; end
    lines = {sprintf('Time: %ss',timeText),sprintf('Frequency window: %s',band)};
    for line = 1:2
        top = cfg.text_top_px+(line-1)*lineSpacing;
        annotation(fig,'textbox',[0.01 1-(top+lineHeight)/cfg.image_size(2) ...
            0.98 lineHeight/cfg.image_size(2)],'String',lines{line}, ...
            'EdgeColor','none','HorizontalAlignment','center','Margin',0, ...
            'VerticalAlignment','top','FontUnits','points','FontSize',labelFontSize, ...
            'FontName','Arial','Interpreter','none');
    end
    drawnow;
    frames{k} = print(fig,'-RGBImage',sprintf('-r%d',exportDpi));
    if cfg.save_frames
        [~,base] = fileparts(files(k).name);
        frameFiles{k} = fullfile(cfg.frames_dir,[base '.png']);
        imwrite(frames{k},frameFiles{k});
    end
    fprintf('%s: frame %d/%d, [%.2f, %.2f] s\n',analysisName,k,numel(files),windows(k,:));
end
if cfg.write_gif
% Build one palette from all frames to prevent frame-specific color shifts.
samples = cellfun(@(rgb) rgb(round(linspace(1,size(rgb,1),100)), ...
    round(linspace(1,size(rgb,2),100)),:),frames,'UniformOutput',false);
[~,palette] = rgb2ind(cat(1,samples{:}),256,'nodither');
if ~isempty(outputDir) && ~isfolder(outputDir), mkdir(outputDir); end
for k = 1:numel(frames)
    indexed = rgb2ind(frames{k},palette,'nodither');
    if k == 1
        imwrite(indexed,palette,cfg.output_file,'gif','LoopCount',Inf,'DelayTime',1/cfg.frame_rate);
    else
        imwrite(indexed,palette,cfg.output_file,'gif','WriteMode','append','DelayTime',1/cfg.frame_rate);
    end
end
end
result = struct('output_file',cfg.output_file,'time_windows',windows, ...
    'input_files',{fullfile({files.folder},{files.name})},'frame_count',numel(files), ...
    'zlim',cfg.zlim,'map_statistic','t','time_centers',mean(windows,2), ...
    'final_width_cm',cfg.final_width_cm,'final_height_cm',cfg.final_height_cm, ...
    'final_font_size_pt',cfg.final_font_size_pt,'header_height_px',headerHeight, ...
    'frame_files',{frameFiles},'gif_written',logical(cfg.write_gif));
end

function cmap = adaptive_cmap(limits)
colors = [59 76 192;103 136 238;154 187 255;221 221 221; ...
    245 158 114;217 88 71;180 4 38]/255;
if limits(1)<0 && limits(2)>0
    positions = [limits(1)*[1 .66 .33] 0 limits(2)*[.33 .66 1]];
elseif limits(1)>=0
    colors = colors(4:7,:); positions = linspace(limits(1),limits(2),4);
else
    colors = colors(1:4,:); positions = linspace(limits(1),limits(2),4);
end
cmap = interp1((positions-limits(1))/diff(limits),colors,linspace(0,1,256),'pchip');
cmap = max(0,min(1,cmap));
end
