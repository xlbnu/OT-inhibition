# MEG analysis

Dependencies: MATLAB R2024b, FieldTrip fieldtrip-20251218, and statistical and parallel computing toolboxes. ICA uses runica bundled with FieldTrip. SPM12 is used for related structural processing. Add the FieldTrip root directory and call ft_defaults, then add the project's MEG and utilities directories.

export_motor_roi_F_gt2_to_nifti.m additionally requires Image Processing Toolbox for bwconncomp. This dependency is for recomputing connected-component NIfTI exports, not for rendering the already supplied source maps.

## Manual execution order

| Stage | File | Inputs and purpose |
|---|---|---|
| 1 | MEG_Preprocessing_workflow.m | Raw MEG, run-level behavior, and trigger information; filtering, trial alignment, repair, and ICA |
| 2 | MEG_headmodel_leadfield.m | T1, fiducial registration, and sensor geometry; generate head models, leadfields, and head-movement information |
| 3 | MEG_SensorLevel_GLM.m | Cleaned MEG and trial behavior; sensor regression and task/ITI coefficients |
| 3 | MEG_sensorLevel_TimeFrequency.m | Cleaned MEG and behavior for 19 participants; conflict/congruent and left/right TFR |
| 4a | MEG_sigmaV_source_localization.m | Behavior, MEG, leadfields, and head movement; ΣV source localization |
| 4b | MEG_gamma_band_source_localization.m | Corresponding data and CSD; gamma source localization |
| 5 | MEG_paired_source_cluster_permutation_F_test.m | Generated source maps; spatial cluster statistics and export |
| 6 | MEG_OXTR_fMRS_overlap.m | Source statistics, OXTR, and fMRS voxel maps; spatial correspondence and overlap |
| 7a | MEG_sigmaV_source_Reconstruction.m | Defined ROIs and matching filters; source signal reconstruction |
| 7b | MEG_gamma_band_source_Reconstruction_TFR.m | Corresponding ROIs and gamma filters; reconstruction and source TFR |
| 8 | MEG_sigmaV_source_Reconstruction_GLM.m | roi_amp and matching trial behavior; time-resolved source GLM |

4a/4b and 7a/7b are branches rather than sequential operations in the same workspace. Complete manual operations such as registration and visual ICA confirmation according to the scripts. Input and output paths appear in the top-level configurations and individual sections; set them before running each script. Results are used in Fig.3 and S5–S7.

## Control analyses and visualization

- MEG_sourceReconstruction_crossTest.m: cross-test signals and spectra using source ROIs; run with each branch's configuration.
- MEG_decision_frequency_source_localization_mutli_window.m: decision-related multi-window source localization.
- MEG_motor_frequency_source_localization_mutli_window.m: motor-related multi-window source localization.
- MEG_multi_slideWindow_test.m: read generated sliding-window results and run decision cluster and motor ROI/FDR statistics.
- MEG_source_Sliding_Window_gif.m: section-based S7A-D example using the scalp/source GIF wrappers. Run Steps 1-6 in the MATLAB Editor: configure this host's toolboxes, select analysisName, export scalp PNGs, export source PNGs, assemble the two GIFs, and verify the outputs. Defaults select S7A (decision_gamma_band). Change analysisName in Step 2 for S7B-D. Intermediate labelled RGB PNGs are retained in figure_outputs/figureS7/frames/<analysisName>/topoplot and source; final GIFs remain in figure_outputs/figureS7. Change frameRate and rerun Step 5 to adjust playback without repeating rendering. The functions retain GIF-only defaults when called directly. See utilities/README.md for options and dependencies.
- export_significant_cluster1_oxtr_to_nifti.m and export_motor_roi_F_gt2_to_nifti.m are export helpers in utilities. run_sliding_source_F_cluster_to_brainnet.m, sliding_motor_roi_fdr_from_source_sets.m, save_roi_fdr_source_map_for_brainnet.m, and export_cluster1_from_stat_file.m are in utilities; add utilities before running the sliding-window statistics scripts.

Raw MEG, structural MRI, leadfields, and intermediate source sets are external inputs. Configure their paths before upstream analysis and perform interactive registration/ICA at the relevant steps. Prepared source GIFs can be rendered directly from the supplied NIfTIs. See [the complete file index](../SCRIPT_INDEX.md) for script purposes and input requirements.

## Cross-test reconstruction

The second branch of MEG_sourceReconstruction_crossTest.m applies the gamma-band source-reconstruction procedure and gamma-derived spatial filters to the SigmaV/OXTR conjunction grids. Its reconstructed time-frequency results are supplied as figure_data/MEG_results/decision_ROIsigmaV_sourceTFR.mat for S6P-T. The first branch performs the complementary cross-test: gamma/OXTR grids with SigmaV-informed reconstruction and trajectory GLM, results are supplied as figure_data/MEG_results/trajectory_sourceEFR_GLM_beta_ROIgamma.mat for S5A-D.
