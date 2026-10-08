# Shared statistical and plotting functions

Add this directory to the MATLAB path before running analyses or figures. It contains helpers shared across scripts; local functions already embedded in scripts do not need to be extracted.

| File | Purpose |
|---|---|
| get2_2StatisticsValue.m | 2×2 repeated-measures results |
| getThreeFactorsInteraction_stat.m | Three-factor interaction tests |
| corr_bootstrap.m | Paired-resampling correlation, confidence intervals, and two-sided P values |
| get2CorrInfor.m | Calls corr_bootstrap and reports correlation information |
| mediation_version2.m | Standardized single-mediator model and bootstrap |
| plot_column_significance.m | Column-wise one-sample t-tests and temporal significance markers |
| errorCIbound.m | Regression lines and confidence bands |
| errorbound.m | Shaded error bounds; FSLboost, Apache 2.0 (LICENSE_FSLboost_Apache-2.0.txt) |
| figureTmp1_nonGroup.m | Distributions without cohort grouping; used in Fig.S7K |
| figureX1.m | Figure windows and axes |
| figureTmp1_OT.m / figureTmp1_OT_SinglePair.m | Four-condition/single-pair distributions, paired data, and annotations |
| RemoveBoxplotlines.m | Boxplot element handling |
| gcaf1.m / setFixedPlotArea.m / setAxesForPPT.m / cutBackColor.m | Fonts, plot sizes, and backgrounds |
| export_motor_roi_F_gt2_to_nifti.m | Export exploratory motor-ROI source components with F>2 |
| export_significant_cluster1_oxtr_to_nifti.m | Export significant source clusters intersected with the OXTR map using a 0.5 threshold |
| export_cluster1_from_stat_file.m | Source cluster-statistics file to NIfTI; requires FieldTrip |
| topoplot_sliding_window_gif.m | S7A-D scalp topographies directly to chronological GIFs |
| brainnet_sliding_window_gif.m | S7A-D source NIfTIs directly to chronological BrainNet GIFs |
| png_frames_to_gif.m | Assemble retained RGB PNG frames into a chronologically ordered GIF |
| batch_brainnet_nii_to_tif.m | Static BrainNet NIfTI-to-TIFF rendering; SPM12 required |
| run_sliding_source_F_cluster_to_brainnet.m | Window-wise source cluster-permutation F tests and NIfTI export; configure external source inputs |
| sliding_motor_roi_fdr_from_source_sets.m | Loop over motor source windows and aggregate ROI-FDR results |
| save_roi_fdr_source_map_for_brainnet.m | Motor ROI-FDR statistics, MAT results, and F_masked NIfTI export |

corr_bootstrap defaults to 10000 resamples; get2CorrInfor explicitly uses 100000 and Seed=42. The number of resamples for mediation_version2 is supplied by the caller; Fig.4/S8 use 100000. plot_column_significance defaults to Correction='none'; correction is applied only when fdr is explicitly specified and is not equivalent to a cluster permutation test.

Statistical functions mainly depend on Statistics and Machine Learning Toolbox; NIfTI export uses FieldTrip.

save_roi_fdr_source_map_for_brainnet.m computes motor ROI-FDR statistics and exports F_masked NIfTI maps; additional statistics and masks are stored in its MAT result. Configure FieldTrip on the MATLAB path and provide the source sets and ROI inputs.

## Sliding-window scalp GIFs

Initialize FieldTrip using setup_figure_paths, then call `topoplot_sliding_window_gif('decision_gamma_band')`. The other supported names are `decision_beta_band`, `motor_gamma_band`, and `motor_beta_band`. Inputs are read from figure_data/MEG_results/source_dynamic_gif/<name>/topoplot; NIfTI files in the sibling source directories are not used. The nested motor-gamma frequency structure is handled automatically.

Only GIFs are written, to figure_outputs/figureS7/<name>_topoplot.gif; rerunning replaces the same GIF. Fifteen windows are played from [-0.80,-0.60] s to [-0.10,0.10] s, sorted by the times in MAT filenames and checked against stored window metadata where available. Each frame displays two header lines, for example "Time: -0.7s" and "Frequency window: [60 80]Hz". Time is the midpoint of the analysis window. The default playback rate is one frame per second, with continuous looping. Optional cfg fields include input_dir, output_file, frame_rate, zlim, image_size (default [700 750], width/height in pixels), final_height_cm (5), and final_font_size_pt (8). Font settings are calibrated as follows: font size is calibrated to 8 pt at a final image height of 5 cm, line spacing is 0.95 times the pixel font size, and the small top/bottom text margins scale from the original 18/8 pixels. Thus, at 4 cm image height the default text is approximately 6.4 pt. The top padding is expanded only as needed to fit the two lines, keeping the header close to the scalp map. Export dimensions are fixed independently of screen DPI. Set top_padding_px, text_top_px, text_bottom_margin_px, or line_spacing_ratio to adjust spacing. Supplying final_width_cm overrides the height-based font calibration; for example final_width_cm=4 and final_font_size_pt=8 restores 8 pt at 4 cm image width. Keep the complete image uncropped and lock its aspect ratio in PowerPoint.

Channel powers are aligned by channel label before calculation. Participant cell order must match across the four condition arrays. Each map displays t-values, not raw power, from a one-sample t-test of the participant-wise average OT/PL log-power contrast: conflict minus congruent for decision, right minus left for motor. Uncorrected P<0.05 controls channel opacity; it does not indicate cluster-corrected significance. Within each GIF, all frames share the same t-value limits and one palette. Default limits are decision gamma [-3,3], decision beta [-7,7], motor gamma [-6,6], and motor beta [-3,6]. Set cfg.zlim to override the default. These are display limits, not statistical thresholds; values outside them saturate at the end colors without being excluded. Gamma labels indicate the supplied 60-80 Hz topography data, not the 30-100 Hz source maps. No spectral re-estimation is performed. MATLAB, FieldTrip, and Statistics and Machine Learning Toolbox are required. The returned result includes zlim and map_statistic ('t').

```matlab
names = {'decision_gamma_band','decision_beta_band', ...
         'motor_gamma_band','motor_beta_band'};
cfg = struct('frame_rate',2); % Two frames per second; omit cfg for defaults.
for k = 1:numel(names)
    result = topoplot_sliding_window_gif(names{k},cfg);
end
```

## BrainNet TIFF export

For source GIFs instead of TIFFs, use `brainnet_sliding_window_gif`; see the next section. The TIFF helper remains available for static panels.

Use cfg.output_suffix to distinguish views of the same NIfTI in one output directory, for example '_medial_left' and '_lateral'. The TIFF name is the NIfTI basename plus this suffix and '.tif'. The default empty suffix preserves the original filename. Fig.S6K uses distinct suffixes for its two views. This suffix affects TIFF filenames; use separate output directories if CSV/MAT sidecars from multiple runs must also be retained.

The bundled BrainNet-Viewer directory supplies BrainNet.m, BrainNet_MapCfg.m, BrainNet.fig, and hemisphere surfaces. Configure an external SPM12 root with the second setup_figure_paths argument, or supply cfg.spm_dir explicitly. By default the renderer resolves SPM from the MATLAB path and BrainNet from this package, without author-specific installation/cache paths. Supply cfg.cfg_file for the intended saved display configuration.

For batch_brainnet_nii_to_tif, set `cfg.save_tif_only = true` to save only TIFF images. CSV manifests and MAT settings are still returned in memory through result but are not written to disk; temporary fixed-scale configurations are automatically removed. Existing output files are not deleted. The default, false, preserves TIFF, CSV, and MAT output. This option does not change the display settings or image format.

## Sliding-window source GIFs

Run from the package root after configuring SPM12:

```matlab
setup_figure_paths([], '/your/local/path/to/spm12')
names = {'decision_gamma_band','decision_beta_band', ...
         'motor_gamma_band','motor_beta_band'};
for k = 1:numel(names)
    result = brainnet_sliding_window_gif(names{k});
end
```

Each call reads frame_*.nii and BrainNet_cfg.mat from figure_data/MEG_results/source_dynamic_gif/<name>/source. The saved view, map projection, and colormap are preserved. Default bundled surfaces are BrainMesh_ICBM152Left.nv for decision gamma, BrainMesh_ICBM152Right.nv for decision beta, and BrainMesh_ICBM152.nv for both motor groups; cfg.surface_file can override the choice. MATLAB, the bundled BrainNet-Viewer files, and SPM12 are required; FieldTrip is not required for this renderer. Configure SPM12 using setup_figure_paths or cfg.spm_dir on each host.

If SPM12 is already on the MATLAB path, run `setup_figure_paths` without an SPM argument; the renderer resolves its root using `fileparts(which('spm'))`. If SPM12 is installed but is not on the MATLAB path, supply that host's root directory (the folder containing spm.m) once per session using `setup_figure_paths([],spmRoot)`. Alternatively, pass `cfg.spm_dir=spmRoot` to the renderer. Installing or copying SPM12 alone does not make MATLAB discover its directory. Add the root rather than using genpath, which can introduce conflicting toolbox subfolders. Check `which spm -all` and `which spm_vol -all` if multiple SPM versions are installed. SPM12 must include compiled components compatible with the host OS; the renderer does not scan drives for installations. Keep the complete OT_code folder together when moving it.

| Analysis | Surface | Positive display limits | Frequency label |
|---|---|---|---|
| decision_gamma_band | Left hemisphere | [5 15] | [30 100]Hz |
| decision_beta_band | Right hemisphere | [5 25] | [13 30]Hz |
| motor_gamma_band | Whole brain | [1 5] | [30 100]Hz |
| motor_beta_band | Whole brain | [5 15] | [13 30]Hz |

Frames are sorted by the actual window-start times encoded in filenames, accepting both p0p10 and 0p10 for positive values. The 15 windows run from [-0.80,-0.60] to [-0.10,0.10] s. The header shows the window midpoint (for example Time: -0.7s) and the full frequency range on the second line. Source gamma labels refer to 30-100 Hz, whereas the scalp gamma GIFs use the separate 60-80 Hz data. The function displays the supplied NIfTI values without recomputing source statistics or changing statistical masks. In particular, the motor-gamma F_gt2 maps retain their supplied selection and do not acquire corrected significance through rendering. Positive BrainNet display limits are fixed across frames, with values below the lower limit shown using the saved null color and high values saturated at the upper color. cfg.zlim overrides the defaults.

Only figure_outputs/figureS7/<name>_source.gif is retained; an existing GIF with that name is replaced. A temporary MAT configuration is automatically removed, including on errors. Original NIfTIs and saved configurations are not modified. Intermediate renderings and a shared 256-color GIF palette stay in memory. Empty maps still contribute a frame to preserve time order. The default frame rate is 1 frame/s with continuous looping.

Defaults use a 700 x 700 output canvas. The default crop corresponds to: 100 pixels from the left, width 1670, and height 1400 on a 2000 x 1500 rendering. The crop is fitted below a compact two-line header without changing its aspect ratio. Use cfg.crop_fraction=[0 0 1 1] to retain the full rendering, or change [left top width height] fractions if the saved view needs another crop. No automatic per-frame crop is applied, so the brain remains in the same position across frames.

Font calibration follows the scalp GIF: final_height_cm=5 and final_font_size_pt=8 (approximately 6.4 pt at 4 cm image height). Optional final_width_cm overrides height calibration. Other options include cfg_file, input_dir, output_file, brainnet_dir, surface_file, spm_dir, frame_rate, zlim, image_size, render_width_px, line_spacing_ratio, text_top_px, and text_bottom_margin_px. Example: `cfg=struct('zlim',[5 20],'frame_rate',2); result=brainnet_sliding_window_gif('decision_gamma_band',cfg);`.

## Staged example with retained PNG frames

Open MEG/MEG_source_Sliding_Window_gif.m and run its six sections in order. Step 1 configures the current host, Step 2 selects S7A-D, Steps 3/4 export scalp/source frames, Step 5 assembles GIFs, and Step 6 checks counts and matching time windows. The script defaults to decision_gamma_band. It uses the same display limits, surfaces, header text, and statistical calculations as the direct wrappers.

Both wrappers accept save_frames (default false), frames_dir (default under figure_outputs/figureS7/frames/<name>/<type>), and write_gif (default true). Set save_frames=true and write_gif=false to render labelled RGB PNGs first; the returned frame_files list identifies the files written by that call, and gif_written is false. These PNGs are the complete cropped/labelled frames used for the GIF, not unlabelled raw maps. No existing PNGs or unrelated outputs are deleted; files with the same names are overwritten. The direct default calls still retain only GIFs. A previous GIF is left untouched when write_gif=false.

```matlab
cfg = struct('save_frames',true,'write_gif',false);
r = topoplot_sliding_window_gif('decision_gamma_band',cfg);
g = png_frames_to_gif(r.frame_files,r.output_file,1);
% Change playback to 2 frames/s without rendering again:
g = png_frames_to_gif(r.frame_files,r.output_file,2);
```

The same pattern applies to brainnet_sliding_window_gif. png_frames_to_gif accepts either the returned file list or a directory containing frame_*.png. It sorts by window-start times in filenames, validates matching RGB dimensions, and uses one shared palette. Prefer the returned file list within a session to avoid including stale files in the directory. It does not alter PNGs and writes only the requested GIF. The function supports the RGB PNGs exported by these wrappers; arbitrary indexed/transparent images require conversion first.
