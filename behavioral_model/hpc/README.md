# Repeated fitting and recovery entry points

These functions use task environment variables for condition-by-repeat computation. They support SLURM jobs and local MATLAB execution with the same variables set. Dependencies are MATLAB R2024b, Global Optimization Toolbox, and Parallel Computing Toolbox.

| Entry function | Input | Output directory under analysis_outputs |
|---|---|---|
| run_slurm_fit_repeated | figure_data/behavioral_results/bhv_data_85.mat | model_fitting/full/RUN_TAG |
| run_slurm_fit_2lambda | Same behavioral data | model_fitting/two_lambda/RUN_TAG |
| run_slurm_compare16 | Same behavioral data; four candidate models | model_comparison/RUN_TAG |
| run_slurm_recovery | Generated prediction_data_85.mat, or the supplied prediction file | model_recovery/RUN_TAG |

Paths are relative to the saved package location. Each task writes participant checkpoints and results_all.mat into condition/repeat folders; comparison additionally includes the model folder. Recovery first checks analysis_outputs/model_prediction_results/prediction_data_85.mat, then figure_data/model_results/prediction_data_85.mat. Keep the entire package structure together when moving it.

## Task configuration

Set SLURM_ARRAY_TASK_ID, SLURM_CPUS_PER_TASK, FIT_NREPEAT, FIT_NSUB, and RUN_TAG before calling an entry function. SLURM_JOB_ID names the worker storage directory. The assertions in each entry function define allowed values. Set FIT_NSUB=85 for a complete dataset and FIT_NREPEAT=30 for 30 fitting repetitions. Task counts are 4 × FIT_NREPEAT for full fitting, two-direction fitting, and recovery, and 16 × FIT_NREPEAT for comparison. Each task represents one condition and one repeat, plus one model for comparison.

For a local first-task example, after adding behavioral_model/hpc to the MATLAB path:

```matlab
setenv('SLURM_ARRAY_TASK_ID','1');
setenv('SLURM_CPUS_PER_TASK','5');
setenv('SLURM_JOB_ID','local');
setenv('FIT_NREPEAT','30');
setenv('FIT_NSUB','85');
setenv('RUN_TAG','full_30repeats');
run_slurm_fit_repeated;
```

This runs only task 1, not the complete experiment. To complete a run, invoke every task ID with consistent settings and the same RUN_TAG. Use separate RUN_TAG values for independent configurations. Worker counts depend on the task's available CPU and participant settings.

Random seeds and input/code hashes protect resumed runs. Use a new RUN_TAG after changing input files or model/entry-point code. After complete full/two-direction/comparison runs, aggregate results using model_best_fit with the matching bestFitAnalysis value. Recovery results remain in the task output files.

See the [complete file index](../../SCRIPT_INDEX.md).
