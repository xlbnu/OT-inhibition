# fMRS spectral processing and concentration aggregation

This directory contains spectral processing, run-level concentration aggregation, and related analysis code. Spectral processing uses GANNET3.5.2 and SPM12; concentration aggregation reads existing fit results.

## Concentration analysis order

| Step | File | Inputs and outputs |
|---|---|---|
| 1 | neurotransmitter_stat.m | Loads four-condition behavior from OT_MRS_data.mat and AllData from Aggregated_MRS_Results.mat, then extracts and aggregates run-level concentrations |
| 2 | figures/figure4_code.m and figureS8_code.m | Load concentrations, behavioral GLMM results, and model parameters; run result-plotting and mediation-analysis sections |

Concentration aggregation uses MATLAB R2024b; statistics and plotting use utilities and Statistics and Machine Learning Toolbox.

## Run-level concentration aggregation

AllData is organized by subject/session/run/context. neurotransmitter_stat.m reads Conflict_Phase2_Fit, Consistent_Phase2_Fit, Conflict_Phase1_Fit, and Consistent_Phase1_Fit, matches the four conditions using mrs_num and expOrder, and generates gaba/glx/glu structures with participant × 2 runs fields. gaba and glx are saved in analysis_outputs/MRS_results/Neurotransmitter_data.mat relative to the package root. The prepared figure_data/MRS_results/Neurotransmitter_data.mat remains the input for figure reproduction.

cf/cg denote phase2 conflict/congruent, and cf1/cg1 denote phase1. The concentration field is ConcIU_AlphaTissCorr. Set behavioral-data and fit-result load paths to their actual locations.

Prepared Fig.4/S8 reproduction reads the supplied concentration structures and does not require a local GANNET installation. The spectral-processing branch requires the complete external GANNET environment and acquisition inputs.

## Spectral processing code

| File | Purpose |
|---|---|
| kh06_config.m | Set raw-data and behavioral-file paths |
| run_KH06.m | Batch processing by subject/run and four subsets |
| GannetLoad_MultiSite.m | Custom GANNET loader that loads frame ranges synchronized to behavior |
| LICENSE_Gannet_BSD-3-Clause.txt | License for the related third-party derivative code |

run_KH06 reads a behavioral file containing sub_data, TWIX data, a water reference, and T1 images. It constructs phase1/phase2 frame ranges using valid and conflict flags, then calls GannetLoad_MultiSite, GannetFit, GannetCoRegister, GannetSegment, and GannetQuantify. TR=1.6 s, with four frames per trial by default. Configure GANNET and SPM12 in the spectral processing environment.

See [the complete file index](../SCRIPT_INDEX.md) for script purposes and input requirements.
