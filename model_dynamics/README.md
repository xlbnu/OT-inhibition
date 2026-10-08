# Model trajectories for MEG

This module constructs a two-stage empirical-Bayes parameter approximation from independently fitted behavioral parameters and generates trial-level trajectories for the 19 MEG participants. It requires MATLAB R2024b, Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox.

## Supplied inputs

Paths are resolved from the saved script location. Keep the package directories together when moving the project.

| Directory under figure_data | File | Purpose |
|---|---|---|
| behavioral_results | bhv_data_85.mat | Four-condition trial data for all 85 retained participants |
| behavioral_results | bhv_data_MEG_19.mat | MEG participant trial data, corresponding to rows 38:56 of the 85-participant dataset |
| model_results | model_fit_parameter_85.mat | Full-model parameters for EB estimation |
| model_results | model_fit_parameter_MEG19.mat | Independent MEG full-model estimates, validated by the input checker |
| model_results | model_fit_parameter_85_NOLambda.mat | Independently refitted no-inhibition parameters for EB estimation |
| model_results | model_fit_parameter_MEG19_NOLambda.mat | Independent MEG no-inhibition estimates for optional direct19 trajectories |

Behavior files contain ots_data, pls_data, otn_data, and pln_data. The parameter order is `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda,tau,t0]`. Full-model arrays are named ots_minb1/pls_minb1/otn_minb1/pln_minb1 and contain optParams. No-inhibition arrays use the minb2 suffix and optParams_FullArray, with lambda fixed at zero.

## Execution order

After adding model_dynamics to the MATLAB path, run:

```matlab
check_trajectory_inputs;
run_hierarchical_eb_85_to_meg19;
run_export_trial_accumulator_glm_data;
run_export_no_lambda_trial_accumulator_data;
model_trajectory_GLM;
```

Each stage writes timestamped products under analysis_outputs/model_dynamics. The full-model export selects the latest saved hierarchy. The final GLM selects the latest full-model and hierarchical85 no-inhibition accumulator exports and saves model_DDM_trajectory_GLM.mat in that output directory. Missing upstream products trigger an instruction identifying the stage to run first. Prepared figure_data inputs are preserved.

For a fixed hierarchy, set trialAccumulatorOverride.hierarchyFile before running the full-model export. For fixed final-GLM inputs, set fullTrajectoryFile and noLambdaTrajectoryFile before running model_trajectory_GLM. These explicit files also permit selection of a direct19 no-inhibition export. Use matching configurations and participant order when combining products.

The no-inhibition export defaults to parameterMode='hierarchical85'. Set noLambdaOverride.parameterMode='direct19' to use the independently fitted 19-participant point estimates. Override structures also permit explicit parameter/behavior files and simulation settings; their fields are documented by the driver scripts. Drivers clear workspace variables except their corresponding Override structure. run_all_local is an optional entry point for the input check and three generation stages; it excludes the final GLM.

## Functions and default settings

| Function | Purpose |
|---|---|
| check_trajectory_inputs | Validate all supplied behavior and parameter files without fitting |
| fit_hierarchical_eb_parameters_85 | Build the EB parameter hierarchy |
| simulate_hierarchical_meg_posterior_trajectories | Generate MEG posterior-predictive trajectories |
| build_hierarchical_meg_trial_trajectory_bank | Generate trial-specific stay/change trajectories |
| build_choice_aligned_trial_accumulator_data | Align trajectories to observed choices, normalize signals, and retain reliability information |

Defaults are 200 posterior draws for EB estimation, 30 parameter draws × 100 paths per trial for trajectory generation, maxSteps=1000, nGrid=101, b0=300, delta=0.01, and four workers. Valid trials require conflict==1, legal choice/confidence values, and finite evidence and ischange. Social evidence is conf/16-1/32; asocial evidence is inforate-0.5. Outputs retain alignment, normalization, path counts, and reliability information.

See the [complete file index](../SCRIPT_INDEX.md) for all script purposes.
