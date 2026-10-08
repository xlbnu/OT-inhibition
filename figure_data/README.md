# Figure input data

These prepared behavioral inputs and analysis outputs supply the MATLAB scripts in ../figures. Keep filenames and relative subdirectories unchanged. They are not a complete raw-acquisition dataset; figure reproduction does not require rerunning upstream analyses.

| Directory | Contents |
|---|---|
| behavioral_results | Accuracy/confidence summaries, evidence curves, cohort-specific and pooled GLMM results, prepared behavior |
| model_results | Full/two-direction parameters, saved model comparison including bmc, predictions, recovery, model trajectories |
| MEG_results | Sensor/source summaries, source TFRs, trajectory coefficients, mean-field spectra, pairing indices, NIfTI maps |
| MRS_results | GABA/Glx concentration structures and MEG/OXTR/fMRS overlap indices |

Run check_figure_inputs from the package root to check files and explicitly loaded top-level variables. The file-to-script manifest ../figures/figure_inputs.json lists direct MAT loads. ../figures/figure_assets.json additionally lists BrainNet code/resources, hemisphere surfaces, static NIfTIs/configurations, and the 15 scalp/source frames plus BrainNet_cfg.mat for each of the four S7 groups. These checks do not validate every nested field, NIfTI geometry, participant ordering, or scientific selection.

Four-condition behavior is ordered as ots_data, pls_data, otn_data, pln_data: OT-social, PL-social, OT-asocial, PL-asocial. The pooled sample contains 85 participants (37 behavioral, 19 MEG, 29 fMRS); MEG participants occupy 38:56 in the pooled ordering. Maintain participant correspondence when replacing data.

Full-model parameter order: [v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda,tau,t0]. Two-direction parameter order: [v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda1,lambda2,tau,t0]. Plotting uses optParams; required panel-specific log transforms are applied in scripts.

MEG_results/pair_data.mat and pair_data_low.mat contain one-based order_idx1/order_idx2 for 19 selected levels in each inhibition group. Retain saved indices for reproduction.

MEG_results/source_dynamic_gif/<analysisName>/topoplot contains prepared scalp power MAT files; source contains prepared source NIfTIs and its BrainNet_cfg.mat. Scalp gamma inputs represent 60-80 Hz, whereas source gamma inputs represent 30-100 Hz. Supplied maps retain their original statistical selections; display ranges do not recompute significance. See ../utilities/README.md for direct GIF rendering and intermediate-frame options.

Some supplied files support upstream or extended analyses. Panel-specific contrasts, coefficient indices, and time/frequency selections are defined in the corresponding figure scripts. Recomputed behavioral and prediction outputs are saved separately under analysis_outputs.

## SigmaV/OXTR ROI time-frequency data

decision_ROIsigmaV_sourceTFR.mat contains time-frequency results reconstructed from source grids in the SigmaV source cluster intersected with the OXTR ROI. Reconstruction uses the gamma-band source-reconstruction procedure and gamma-derived spatial filters, with the ROI grids replaced by the SigmaV/OXTR conjunction. The SigmaV source-localization filters are not used for this reconstruction.
