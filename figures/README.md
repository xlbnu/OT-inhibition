# Article figure reproduction

The scripts generate MATLAB panels for Fig.1–4 and Fig.S1–S8 from the supplied ../figure_data inputs. They remain research scripts intended for manual execution by section. Refitting behavioral models, preprocessing raw MEG, or fitting fMRS concentrations is not required for these plots. Statistical tables follow Methods; no separate table-export workflow is included.

## Quick start on another computer

Copy the entire package, keeping figures, utilities, and figure_data under the same root. The root and parent directories can be renamed. Use a location with write permission. Use MATLAB R2024b with Statistics and Machine Learning Toolbox. Complete Fig.3/S6/S7 scripts additionally require FieldTrip fieldtrip-20251218 and SPM12; Fig.S5 requires SPM12. Fig.S3 reads saved bmc results and does not require VBA for replotting; VBA V3.0 remains a dependency for recomputing model comparisons.

Set MATLAB's Current Folder to the package root, then run in the Command Window:

```matlab
setup_figure_paths
report = check_figure_inputs;
figure1_code
```

For Fig.3/S6/S7, configure the installed FieldTrip and SPM12 roots first:

```matlab
setup_figure_paths('/your/local/path/to/fieldtrip-20251218','/your/local/path/to/spm12')
figure3_code
```

Replace the example tool path with your actual installation directory. setup_figure_paths adds figures and utilities, and initializes FieldTrip when supplied or already on the MATLAB path. It does not change the current directory or run analyses. Repeat setup in a new MATLAB session.

check_figure_inputs reads figure_inputs.json for direct MAT loads and figure_assets.json for BrainNet/S7 assets, and verifies the requested top-level variables without running statistics. MATLAB supports both ordinary and v7.3 MAT files. Update the manifests when changing script inputs. A successful check does not validate nested fields, participant matching, or complete panel rendering.

Alternatively, open a saved figure script and run it in the MATLAB Editor. Each script locates the package root using mfilename('fullpath'), adds utilities, and creates its output directory. Do not paste its path-discovery lines into the Command Window: no script is executing there, so mfilename is empty. Configure FieldTrip before running scripts that require it.

| Script | Main contents and prerequisite results |
|---|---|
| figure1_code.m | Pooled behavior, synergy, curves, and GLMM |
| figure2_code.m | Full-model parameters, prediction coefficients, and correlations |
| figure3_code.m | MEG-cohort behavior, model trajectories, mean-field results, sensor/source signals and spectra |
| figure4_code.m | fMRS-cohort behavior, parameters, GABA, and mediation |
| figureS1_code.m | Cohort-specific behavioral results |
| figureS2_code.m | Perceptual-stage and congruent-trial control analyses |
| figureS3_code.m | Model comparison, recovery, and predictive validation |
| figureS4_code.m | Two-direction inhibition model |
| figureS5_code.m | ROI cross-tests and source signals |
| figureS6_code.m | Mean-field beta and spectra from different source ROIs |
| figureS7_code.m | Decision- and motor-related source/spectral controls |
| figureS8_code.m | Stage, trial-type, metabolite, and mediation controls |

## Section order and outputs

Run initialization, data-loading, and shared-variable sections before the required panels. Fig.3H/I require the trajectory computations in G; S5C/D require B; S7K requires the preceding decision/motor TFR computations. Sections within a group share workspace variables. Running a complete script from an empty workspace avoids reliance on a different figure's variables. Local functions remain in their script files.

Outputs are written to ../figure_outputs/figure1, figure2, ..., figureS8. Directories are created automatically; rerunning a panel overwrites outputs with the same names. f_file retains a trailing platform-specific separator for existing export calls. SVG filenames are determined by those calls. Scripts export panels; manuscript layout and schematics require separate assembly.

## Brain surface panels

Fig.3F/T render automatically with the bundled utilities/BrainNet-Viewer files and configured SPM12. F uses decision_sigmaV_source_OXTRconjunction_.nii with BrainMesh_ICBM152Right.nv and BrainNet_cfg_medial_right.mat; T uses decision_gamma_band_source_OXTRconjunction_.nii with BrainMesh_ICBM152Left.nv and BrainNet_cfg_medial_left.mat. Maps and configurations are in figure_data/MEG_results. Saved thresholds, color limits, and image settings are preserved; only TIFFs are saved in figure_outputs/figure3. The two outputs retain their source-map filenames.

Keep the BrainNet source, BrainNet.fig, hemisphere surfaces, NIfTI maps, and saved EC configurations when migrating the package. SPM12 is an external dependency: configure it with the second setup_figure_paths argument or add its root to the MATLAB path. check_figure_inputs checks the MAT-load and BrainNet/S7 asset manifests.

S5A exports its static brain map automatically. S6F/K/P also export brain maps, with distinct suffixes for the two S6K views. S7A-D call topoplot_sliding_window_gif for scalp maps and brainnet_sliding_window_gif for source surfaces within the figure script; see utilities/README.md. Call either function with decision_gamma_band, decision_beta_band, motor_gamma_band, and motor_beta_band for S7A, B, C, and D, respectively. The animation calls save only GIFs in figure_outputs/figureS7, with 15 chronologically ordered, midpoint-time/frequency-labelled frames per animation; later S7 sections export static SVG panels. The source renderer reads each group's source/BrainNet_cfg.mat and requires SPM12; scalp maps require FieldTrip. Gamma frequency labels are 30-100 Hz for sources and 60-80 Hz for scalp maps.

## Randomness and shared helpers

For a staged S7A-D example with intermediate PNG export and separate GIF assembly, run MEG/MEG_source_Sliding_Window_gif.m section by section. PNGs are saved under figure_outputs/figureS7/frames/<analysisName>/topoplot and source. The direct wrapper calls retain their GIF-only defaults.

Mean-field plots load one saved random pairing of high- and low-inhibition groups: Fig.3 reads pair_data.mat and S6 reads pair_data_low.mat, both included in figure_data/MEG_results. Their indices are used in the corresponding paired calculations. Do not regenerate them for replotting.

Scatter jitter, mediation bootstrap, and frequency-domain permutation tests also use random draws. Use saved pairing indices and the seeds defined in the scripts for repeat analyses. Correlation and trajectory bootstrap functions use their specified seeds.

errorbound.m (FSLboost; Apache 2.0 license in utilities) and figureTmp1_nonGroup.m are included in utilities. errorCIbound.m is a separate regression confidence-band function. FieldTrip's neuromag306cmb.lay is used for sensor topographies.

## S6P-T ROI and reconstruction

decision_ROIsigmaV_sourceTFR.mat contains time-frequency results reconstructed from source grids in the SigmaV source cluster intersected with the OXTR ROI. Reconstruction uses the gamma-band source-reconstruction procedure and gamma-derived spatial filters, with the ROI grids replaced by the SigmaV/OXTR conjunction. The SigmaV source-localization filters are not used for this reconstruction. The S6P-T sections load this file from figure_data/MEG_results.
