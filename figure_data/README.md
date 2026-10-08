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

## Behavioral model inputs

behavioral_results/bhv_data_85.mat contains four-condition trial data for all 85 retained participants. bhv_data_MEG_19.mat contains the corresponding MEG participants in the same order as rows 38:56 of the complete dataset. These files support behavioral GLMs and computational-model fitting.

model_results/model_fit_parameter_85_NOLambda.mat and model_fit_parameter_MEG19_NOLambda.mat contain independently refitted no-inhibition model parameters. Their four minb2 arrays contain optParams_FullArray with nine slots and lambda=0 in slot 7. Full-model parameters use the corresponding model_fit_parameter_85.mat and model_fit_parameter_MEG19.mat files and optParams fields. New analyses write to analysis_outputs, preserving the supplied figure inputs.

## Sliding-window scalp input format

Both motor_gamma_band/topoplot and motor_beta_band/topoplot use four top-level cell arrays: ot_freq_right, ot_freq_left, pl_freq_right, and pl_freq_left. Each cell contains one participant's FieldTrip channel-power structure. Motor-gamma files also contain subject_ids (participant order) and t_win (window boundaries). These are the only six top-level variables in motor-gamma scalp files; CSDs, fitting/export configurations, and trial-Fourier flags are not required for GIF reproduction. Decision inputs retain ot_freq_conf, ot_freq_cong, pl_freq_conf, and pl_freq_cong.

The topoplot_sliding_window_gif function uses one shared reading branch for both motor bands. Filenames determine chronological order, and t_win is checked where supplied. Power values, channel order, participant order, and the participant-wise t-statistic calculation are preserved.

## Participant IDs and cohort rows

participant_index.csv defines the public participant IDs and the corresponding rows in the complete and cohort-specific datasets. IDs follow the 85-row retained-sample order: sub-001 to sub-037 for the behavioral cohort, sub-038 to sub-056 for MEG, and sub-057 to sub-085 for fMRS. Single-cohort files retain these global IDs. The four behavioral conditions contain one own-participant subsName entry per record; partner names, expOrder, and pair_num are omitted. bhv_data_bhv_37.mat and bhv_data_MRS_29.mat provide the additional cohort subsets. Numerical parameter, coefficient, and concentration arrays retain their supplied participant order.
