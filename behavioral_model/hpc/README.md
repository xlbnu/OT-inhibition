# SLURM model computation entry points

The four entry points perform repeated full-model fitting, two-direction fitting, four-model comparison, and parameter recovery. They require a SLURM-configured computing environment; local MATLAB analyses do not require these entry points.

| File | Input |
|---|---|
| run_slurm_fit_repeated.m | bhv_data_85.mat and main_model.m in the parent directory |
| run_slurm_fit_2lambda.m | bhv_data_85.mat and main_model_2lambda.m in the parent directory |
| run_slurm_compare16.m | bhv_data_85.mat and main_model_compare_one.m in the parent directory |
| run_slurm_recovery_2k.m | prediction_data_85.mat and model_recovery_2k.m in the parent directory |

Before submission, set SLURM_ARRAY_TASK_ID, SLURM_CPUS_PER_TASK, FIT_NREPEAT, FIT_NSUB, and RUN_TAG. SLURM_JOB_ID is used for the parallel temporary directory. Allowed ranges are defined by each entry point's assert statements. Entry points select conditions and repeats by task and write participant checkpoints and results_all.mat under outputs/RUN_TAG. Random seeds and input/code hashes are used to verify resumed runs; use a new RUN_TAG if inputs or code change.

Entry points add their parent behavioral_model directory to the MATLAB path and hash model files from that directory. Inputs are read from the hpc directory if present; otherwise they fall back to figure_data/behavioral_results/bhv_data_85.mat or figure_data/model_results/prediction_data_85.mat under the package root. Keep the package structure together on the cluster. Results are aggregated by model_best_fit.m. Changes to entry-point code also affect checkpoint signatures; use a new RUN_TAG after updating code.

See [the complete file index](../../SCRIPT_INDEX.md) for script purposes.
