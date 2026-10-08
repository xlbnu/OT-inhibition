# Behavioral computational models

This module contains fitting functions and postprocessing scripts that depend on workspace variables. Dependencies: MATLAB R2024b, Global Optimization, Statistics and Machine Learning, and Parallel Computing toolboxes. Model comparison requires VBA V3.0; plotting requires utilities.

## Files and execution order

| Order | File | Usage and prerequisite inputs |
|---|---|---|
| 1 | main_model.m | Function; takes a participant structure for one condition and returns optParams/nll/bic |
| 1 | main_model_2lambda.m | Function; takes the same type of behavioral input and fits inhibition in two directions |
| 1 | main_model_compare_one.m | Function; takes behavioral input and modelType=1:4 |
| 2 | model_best_fit.m | Read repeated-fit results by section, select parameters, and aggregate model comparisons; prepare the corresponding results separately for the full and two-direction branches |
| 3 | main_model_prediction.m | Preload four-condition behavior, ots_minb1/pls_minb1/otn_minb1/pln_minb1, and ot_valid into the workspace; generates predictions and GLMM results |
| 4 | model_recovery_2k.m | Function; takes simulated data containing SelfConfidence/OtherConfidence/InitialChoice/PredictedChoice/PredictedConfidence/PredictedRT |
| 5 | model_parameter_equivalence_test.m | Preload four participant × parameter matrices from the same model, ots_xx1/pls_xx1/otn_xx1/pln_xx1; run by section |

See the [HPC instructions](hpc/README.md) for cluster entry points for repeated fitting and recovery. Locally, call fitting functions for each condition and aggregate repeated-fit results.

## Parameters and outputs

The symmetric model has nine slots: `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda,tau,t0]`. The two-direction model has ten slots: `[v0,b1_ratio,b2_ratio,k1,k2,sigma,lambda1,lambda2,tau,t0]`; direction definitions follow the model equations. modelType specifies full, no-lambda, no-urgency, and basic models, respectively, with 9/8/7/6 free parameters.

model_best_fit generates workspace variables for the best parameters and mc_aic/mc_bic/mc_nll in model_compare_data.mat. Configure its external repeated-fit input directories before running sections. main_model_prediction.m saves prediction_data_85.mat and all_GLM_beta_pred_2k.mat in analysis_outputs/model_prediction_results relative to the package root. It does not overwrite the prepared figure inputs. The three prediction GLMM columns are the intercept, ΣV, and ΔV. The choice response is isstay, with stay=1 and change=0; cohort terms are included in the pooled prediction GLMM. Predicted RT is retained in simulation output but is not analyzed by the prediction GLMM. Simulation and mixed-model helpers are local functions within the script. Results are used in Fig.2, S3, S4, and model_dynamics. The figure scripts read figure_data/model_results/all_GLM_beta_pred.mat; newly computed predictions are saved separately in analysis_outputs/model_prediction_results.

## Manual execution

Load the corresponding inputs separately for the full model, two-direction model, and recovery analyses. Before prediction, load four-condition behavior, fitted parameters, and retained-sample indices. Run best-fit selection, model comparison, and equivalence tests by their respective sections.

The equivalence script contains separate full-model and two-direction-model sections. Load the matching participant-by-parameter matrices and run the selected section. Its SDMultiplier argument sets the equivalence bound. See [the complete file index](../SCRIPT_INDEX.md) for script purposes.

## Equivalence-test bounds

model_parameter_equivalence_test.m implements two analyses using paired standardized differences, dz = mean(difference)/SD(difference):

| Analysis | Standardized equivalence bounds | Script section |
|---|---|---|
| Compare the four conditions within each parameter | [-0.6, 0.6] | equivalence_four_conditions_sd; six paired condition comparisons per parameter |
| Compare different parameters within each condition | [-0.75, 0.75] | equivalence_parameter_levels_partial_conjunction_sd |

For comparisons between different parameters, values are first normalized using the parameter-specific lower and upper bounds supplied in ParameterMin and ParameterMax. The script tests each parameter pair in the four conditions and combines the condition-level results using its configured RequiredConditions value. The supplied configuration uses RequiredConditions=4 and Holm correction across parameter pairs. These bounds apply to different analyses; they are not interchangeable settings for the same test.
