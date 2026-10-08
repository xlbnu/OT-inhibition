# Two-k model trajectories for MEG

This module uses independently fitted parameters from 85 participants for a two-stage empirical-Bayes approximation, then generates corresponding model trajectories for 19 MEG participants. Dependencies: MATLAB R2024b, Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox.

## Inputs and parameters

Parameter order: `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda,tau,t0]`. The full model uses model_minb1_2k_85.mat; the no-lambda model uses independently refitted model_minb2_2k_85.mat, with slot 7 set to zero. Full-model parameter values are in optParams; the no-lambda export branch reads optParams_FullArray.

Behavioral inputs bhv_data_85.mat and bhv_data_MEG_19.mat must contain the four-condition variables ots_data/pls_data/otn_data/pln_data. MEG subject IDs should correspond to 38:56 in the 85-participant data. check_2k_trajectory_inputs also reads MEG19 parameter files for the full/no-lambda models. Configure the driver scripts for the external behavior, fitted-parameter, and trajectory inputs used by this module.

## Manual execution order

1. Add this directory to the path, check script input paths, and run check_2k_trajectory_inputs.
2. Run run_hierarchical_eb_85_to_meg19_2k.m to generate the hierarchy and MEG trajectories.
3. Run run_export_trial_accumulator_glm_data_2k.m to generate full-model trialAccumulator. The latest hierarchy is selected by default; specify a fixed version with trialAccumulatorOverride.hierarchyFile.
4. Run run_export_no_lambda_trial_accumulator_data_2k.m to generate trialAccumulatorNoLambda. The default is parameterMode='hierarchical85'; noLambdaOverride.parameterMode='direct19' uses independent MEG19 point estimates, repeatedly using those point estimates in this mode.
5. In model_trajectory_GLM.m, explicitly select the full/no-lambda exported files and MEG behavioral inputs. Run the difference and stay/change regression sections and save model_DDM_trajectory_GLM.mat.

run_all_2k_local.m is an optional module entry point containing only input checks and steps 2–4, excluding the final GLM. Some driver scripts call clearvars and preserve only the corresponding Override variables; do not rely on path settings from an earlier workspace.

## Core functions and outputs

- fit_hierarchical_eb_parameters_85_2k: EB parameter modeling, with 200 draws by default.
- simulate_hierarchical_meg_posterior_trajectories_2k: MEG predictive trajectories.
- build_hierarchical_meg_trial_trajectory_bank_2k: trial-level trajectory bank.
- build_choice_aligned_trial_accumulator_data_2k: alignment to observed choice, normalization, and reliability flags.

Defaults are 30 parameter draws × 100 diffusion paths per trial, maxSteps=1000, nGrid=101, b0=300, and delta=0.01, with four process workers. Valid trials require legal choice/confidence values, conflict==1, and finite others' evidence and ischange. Social input is conf/16-1/32; asocial input is inforate-0.5. Outputs retain stay/change branches, B0/dynamic-boundary normalization, baseline/endpoint alignment, path counts, and reliability information. They are saved in timestamped subdirectories under outputs in the code directory and used for trajectory analyses in Fig.3.

See [the complete file index](../SCRIPT_INDEX.md) for script purposes and input requirements.
