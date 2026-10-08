%% Step 1. Configure this host once per MATLAB session
% Open this saved script in the MATLAB Editor and run sections in order.
% Keep the complete OT_code directory together when moving to another host.
% Run this section from the file, not by pasting mfilename into the Command Window.
scriptPath = mfilename('fullpath');
assert(~isempty(scriptPath),'Run initialization from the saved script.');
codeRoot = fileparts(fileparts(scriptPath));
addpath(codeRoot);

% Leave empty if already available on the MATLAB path. Otherwise enter this
% host's installation root containing ft_defaults.m or spm.m, respectively.
fieldtripRoot = '';
spmRoot = '';
setup_figure_paths(fieldtripRoot,spmRoot);

%% Step 2. Select a panel and configure display/export settings
% Change analysisName to reproduce another panel, then rerun Steps 2-6.
% Panel  Analysis               Scalp range  Source range  Source surface
% S7A    decision_gamma_band     [-3 3]       [5 15]        left hemisphere
% S7B    decision_beta_band      [-7 7]       [5 25]        right hemisphere
% S7C    motor_gamma_band        [-6 6]       [1 5]         whole brain
% S7D    motor_beta_band         [-3 6]       [5 15]        whole brain
% Gamma scalp data: 60-80 Hz. Gamma sources: 30-100 Hz. Beta: 13-30 Hz.
analysisName = 'decision_gamma_band';
analysisNames = {'decision_gamma_band','decision_beta_band', ...
    'motor_gamma_band','motor_beta_band'};
analysisIndex = find(strcmp(analysisName,analysisNames));
assert(~isempty(analysisIndex),'Select one of the four documented analysis names.');
scalpLimits = [-3 3;-7 7;-6 6;-3 6];
sourceLimits = [5 15;5 25;1 5;5 15];
surfaces = {'BrainMesh_ICBM152Left.nv','BrainMesh_ICBM152Right.nv', ...
    'BrainMesh_ICBM152.nv','BrainMesh_ICBM152.nv'};
inputRoot = fullfile(codeRoot,'figure_data','MEG_results','source_dynamic_gif',analysisName);
outputRoot = fullfile(codeRoot,'figure_outputs','figureS7');
frameRoot = fullfile(outputRoot,'frames',analysisName);

% Render PNGs first; assemble GIFs separately in Step 5. The wrapper defaults
% remain save_frames=false and write_gif=true for direct GIF-only calls.
commonCfg = struct('save_frames',true,'write_gif',false, ...
    'final_height_cm',5,'final_font_size_pt',8,'line_spacing_ratio',0.95);
frameRate = 1; % Frames per second; change this and rerun only Step 5.
previewFrames = false; % Set true to display the first PNGs in Step 6.

%% Step 3. Calculate scalp t maps and export 15 labelled PNG frames
% Requires FieldTrip and Statistics and Machine Learning Toolbox.
% Conflict-congruent (decision) or right-left (motor) log-power contrasts
% are averaged over OT/PL within participants, then tested against zero.
% Uncorrected P<0.05 controls opacity; color represents the t-statistic.
topoCfg = commonCfg;
topoCfg.input_dir = fullfile(inputRoot,'topoplot');
topoCfg.output_file = fullfile(outputRoot,[analysisName '_topoplot.gif']);
topoCfg.frames_dir = fullfile(frameRoot,'topoplot');
topoCfg.zlim = scalpLimits(analysisIndex,:);
topoCfg.image_size = [700 750];
topoResult = topoplot_sliding_window_gif(analysisName,topoCfg);
assert(~topoResult.gif_written,'This step should export PNG frames only.');
fprintf('Scalp PNG frames: %s\n',topoCfg.frames_dir);

%% Step 4. Render source maps and export 15 labelled PNG frames
% Requires the bundled BrainNet Viewer and this host's SPM12 installation.
% Read the saved view/colormap from this group's source/BrainNet_cfg.mat.
% The NIfTI values and their existing statistical selection are unchanged.
sourceCfg = commonCfg;
sourceCfg.input_dir = fullfile(inputRoot,'source');
sourceCfg.cfg_file = fullfile(sourceCfg.input_dir,'BrainNet_cfg.mat');
sourceCfg.surface_file = fullfile(codeRoot,'utilities','BrainNet-Viewer',surfaces{analysisIndex});
sourceCfg.output_file = fullfile(outputRoot,[analysisName '_source.gif']);
sourceCfg.frames_dir = fullfile(frameRoot,'source');
sourceCfg.zlim = sourceLimits(analysisIndex,:);
sourceCfg.image_size = [700 700];
sourceCfg.crop_fraction = [100/2000 0 1670/2000 1400/1500];
% To retain the complete saved view, use crop_fraction=[0 0 1 1].
sourceResult = brainnet_sliding_window_gif(analysisName,sourceCfg);
assert(~sourceResult.gif_written,'This step should export PNG frames only.');
fprintf('Source PNG frames: %s\n',sourceCfg.frames_dir);

%% Step 5. Assemble GIFs from the exported PNGs without rendering again
% Window-start times in filenames determine frame order. A shared palette
% prevents frame-specific color shifts. Existing GIFs of the same name are replaced.
% Use the exact returned file lists to avoid including stale PNGs in a folder.
topoGif = png_frames_to_gif(topoResult.frame_files,topoCfg.output_file,frameRate);
sourceGif = png_frames_to_gif(sourceResult.frame_files,sourceCfg.output_file,frameRate);

% In a later MATLAB session, after initializing codeRoot and the utilities
% path, a frame directory can be supplied instead of an in-memory file list:
% png_frames_to_gif(fullfile(codeRoot,'figure_outputs','figureS7','frames', ...
%     'decision_gamma_band','topoplot'), ...
%     fullfile(codeRoot,'figure_outputs','figureS7','decision_gamma_band_topoplot.gif'),2);

%% Step 6. Check frame counts, time order, and optional previews
assert(topoGif.frame_count==15 && sourceGif.frame_count==15,'Expected 15 frames per GIF.');
assert(max(abs(topoGif.time_windows(:)-sourceGif.time_windows(:)))<1e-6, ...
    'Scalp and source windows do not match.');
assert(all(diff(topoGif.time_windows(:,1))>0),'Frames must be chronological.');
topoInfo = imfinfo(topoGif.output_file);
sourceInfo = imfinfo(sourceGif.output_file);
assert(numel(topoInfo)==15 && numel(sourceInfo)==15,'Unexpected GIF frame count.');
fprintf('Verified %s: 15 matching windows, from [-0.80,-0.60] to [-0.10,0.10] s.\n',analysisName);
fprintf('Scalp GIF: %s\nSource GIF: %s\n',topoGif.output_file,sourceGif.output_file);
if previewFrames
    figure('Name',[analysisName ' scalp: first frame'],'Color','w');
    image(imread(topoGif.frame_files{1})); axis image off;
    figure('Name',[analysisName ' source: first frame'],'Color','w');
    image(imread(sourceGif.frame_files{1})); axis image off;
end