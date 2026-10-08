# Behavioral computational models

This module fits behavioral models, aggregates repeated fits, generates predictions, fits recovery models, and tests parameter equivalence. Dependencies are MATLAB R2024b, Global Optimization Toolbox, Statistics and Machine Learning Toolbox, and Parallel Computing Toolbox. Add utilities to the MATLAB path for statistical and plotting helpers. Bayesian model comparison in downstream figure analyses uses VBA V3.0.

## Inputs

figure_data/behavioral_results/bhv_data_85.mat contains the four-condition arrays ots_data/pls_data/otn_data/pln_data for all 85 retained participants. No additional participant-selection index is required. Their order is 37 behavioral, 19 MEG, and 29 fMRS participants. Fitting functions take one condition's participant structure array, or a single participant structure.

| File | Purpose |
|---|---|
| main_model.m | Fit the nine-parameter full model |
| main_model_2lambda.m | Fit the ten-parameter model with separate inhibition in two directions |
| main_model_compare_one.m | Fit modelType=1:4: full, no-inhibition, no-urgency, and basic models |
| model_best_fit.m | Read completed repeated-fit results, select full/two-direction parameters, or aggregate comparison metrics |
| main_model_prediction.m | Load all 85 participants and full-model parameters; generate single-run simulated data and repeated prediction GLMMs |
| model_recovery.m | Fit recovery models to simulated participant data |
| model_parameter_equivalence_test.m | Run condition-within-parameter and parameter-within-condition equivalence analyses by section |

## Fitting and aggregation

See the [HPC instructions](hpc/README.md) for task configuration and output locations. The entry functions can also run in local MATLAB after their environment variables are configured. Computational cost depends on participants, conditions, repetitions, and worker count.

After fitting, set bestFitAnalysis to 'full', 'two_lambda', or 'comparison', then run model_best_fit. The script selects the latest RUN_TAG directory for that analysis. Set bestFitInputDir to a specific RUN_TAG directory to select a fixed run instead. Only the selected branch executes. Completed input files must contain all 85 participant fits. Full/two-direction selection retains repeats with at most one parameter on its specified fitting boundary, then selects the lowest BIC. Comparison aggregation uses the minimum BIC, AIC, and NLL across repetitions for each model and participant; AIC uses the model's free-parameter count.

Outputs are saved under analysis_outputs/model_results as model_fit_parameter_85.mat, twolambda_model_fit_parameter_85.mat, or model_compare_data.mat. The symmetric model has nine parameter slots `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda,tau,t0]`; the two-direction model has ten slots `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda1,lambda2,tau,t0]`. The comparison models have 9/8/7/6 free parameters, respectively.

## Prediction and recovery

main_model_prediction loads the supplied full-model parameters from figure_data/model_results/model_fit_parameter_85.mat by default. To use a newly aggregated fit, set predictionParameterFile to its complete filename before running. All four behavior and fit arrays must contain 85 aligned participants.

Predictions are written to analysis_outputs/model_prediction_results/prediction_data_85.mat. Repeated-prediction GLMM results are saved there as all_GLM_beta_pred.mat. The three coefficient columns are intercept, ΣV, and ΔV. The choice response is stay=1/change=0; pooled GLMMs include cohort terms. Predicted RT is saved but is not analyzed by these GLMMs. The figure scripts continue to use the supplied prepared results in figure_data/model_results.

model_recovery takes simulated data containing SelfConfidence, OtherConfidence, InitialChoice, PredictedChoice, PredictedConfidence, and PredictedRT. The recovery entry point reads newly generated prediction_data_85.mat when available, otherwise the supplied file in figure_data/model_results. It writes repeated fits under analysis_outputs/model_recovery.

## Equivalence analyses

Run model_parameter_equivalence_test by section after its initialization. The full-model section explicitly loads model_fit_parameter_85.mat; the two-direction section loads twolambda_model_fit_parameter_85.mat, both from figure_data/model_results. The between-parameter section selects the corresponding file through model_n_parameters=9 or 10. These inputs are the supplied publication parameters; change the loading path explicitly to test a newly fitted model.

Bounds use paired standardized differences, dz=mean(difference)/SD(difference). Comparisons among the four conditions within each parameter use [-0.6,0.6]. Comparisons between parameters within each condition use [-0.75,0.75], with parameter-specific normalization bounds, RequiredConditions=4, and Holm correction across parameter pairs. These are separate analyses.

See the [complete file index](../SCRIPT_INDEX.md) for all script purposes.
