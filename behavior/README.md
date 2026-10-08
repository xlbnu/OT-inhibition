# Behavioral data preparation and analysis

Dependencies: MATLAB R2024b, Statistics and Machine Learning Toolbox, and utilities. Run the data-reading and analysis sections manually.

## Execution order

| Step | File or section | Input | Output and purpose |
|---|---|---|---|
| 1 | Three cohort sections in Behavior_data_proprecess.m | Raw OT/PL behavioral records for each cohort; fMRS participant-ID mapping | Social/asocial structures for each cohort |
| 2 | merge three cohort experiments' data | Cohort structures for the four conditions | OT_data.mat; choice, confidence, and log-odds fields |
| 3 | Apply the manually retained sample defined in Methods | Merged data and manually retained indices | Final analysis inputs for 85 participants with matching IDs across conditions |
| 4 | Input section in Behavioral_analysis.m | bhv_data_85.mat, containing ots_data/pls_data/otn_data/pln_data | Check the 85 participants and the 37/19/29 cohort ordering |
| 5 | Descriptive statistics, curves, overall GLMM, and cohort-specific GLMM sections | The same retained sample | Behavioral results in analysis_outputs/behavioral_results |

Manual exclusions follow experimental records and Methods, including abnormal behavior, sleepiness during the experiment, withdrawal, and disordered data. These exclusions are not automatically determined within this script. Use the established retained sample when running the code; do not redefine exclusion criteria.

## Variables and results

Main inputs include answer1/answer2, conf1/conf2, conflict, ischange, rtime1/rtime2, coh, correctness, and others' evidence. Preprocessing calculates lo_s1/lo_o1/lo_s2, lo_plus/lo_minus, lo_cplus/lo_cminus, and related fields. Social evidence is converted from confidence, whereas asocial evidence is converted from inforate.

Behavioral_analysis.m reads figure_data/behavioral_results/bhv_data_85.mat relative to the package root and writes all_accauc.mat, MEG_acc_auc.mat, coherence_data.mat, sigmoid_data.mat, all_GLM_beta_withGroup_fixExpID.mat, singleExp_GLM_beta.mat, MEG_GLM_beta.mat, and MRS_GLM_beta.mat to analysis_outputs/behavioral_results. It does not overwrite the prepared figure inputs. GLMM coefficient columns are the intercept, ΣV, and ΔV; the choice model fits isstay (stay=1, change=0) directly. No sign reversal is needed for the stay direction. These results support Fig.1, S1, S2, and cohort plots. 

## Manual execution

The retained-sample bhv_data_85.mat is supplied. Start directly with Behavioral_analysis.m to reproduce behavioral analyses; no additional participant exclusions are required. Behavior_data_proprecess.m is only needed to rebuild these structures from acquisition records. Each section requires its preceding variables.

Configure the raw cohort input paths in Behavior_data_proprecess.m before data preparation. Prepared plotting inputs remain in figure_data; recomputed summaries are written to analysis_outputs/behavioral_results.
