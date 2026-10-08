# Oxytocin social decision analysis code

This directory contains code for behavioral analyses, computational models, MEG, fMRS, and article figures. Most files are research scripts that depend on workspace variables and are run manually by section. Each module README describes its inputs, execution order, and outputs. Prepared figure inputs and retained-participant behavioral trial data are supplied in figure_data; raw MEG and MRS acquisition data are not included. Statistical tables are provided in the appendix following the article Methods.

## Quick start: reproduce figures

Copy the entire package and set MATLAB's Current Folder to its root. With MATLAB R2024b and Statistics and Machine Learning Toolbox installed, run:

```matlab
setup_figure_paths
report = check_figure_inputs;
figure1_code
```

Select another figure by its script name, for example figure2_code or figureS4_code. For Fig.3/S6/S7, initialize your installed FieldTrip first:

```matlab
setup_figure_paths('/your/local/path/to/fieldtrip-20251218','/your/local/path/to/spm12')
figure3_code
```

Replace the example paths with local installation roots. Complete Fig.3/S6/S7 scripts require both FieldTrip and SPM12; Fig.S5 requires SPM12 for its static brain map. Other figure scripts need MATLAB and Statistics and Machine Learning Toolbox. Figure scripts locate figure_data relative to their own saved location and create figure_outputs/<figure name> automatically. They do not require editing data paths when the complete package is moved or renamed. Run script initialization before individual panel sections. Static brain panels in Fig.3/S5/S6 export TIFFs automatically; Fig.S7A-D generate scalp/source GIFs. See the [figure instructions](figures/README.md), [data notes](figure_data/README.md), and [complete per-file index](SCRIPT_INDEX.md).

## Directories and execution order

| Directory | Contents | Instructions |
|---|---|---|
| behavior | Behavioral data preparation and analysis | [README](behavior/README.md) |
| behavioral_model | Fitting, result aggregation, prediction, recovery, and equivalence testing | [README](behavioral_model/README.md) |
| behavioral_model/hpc | SLURM entry points | [README](behavioral_model/hpc/README.md) |
| model_dynamics | Empirical-Bayes parameter modeling and model trajectories for the MEG cohort | [README](model_dynamics/README.md) |
| mean-field | Mean-field network and spectral simulations | [README](mean-field/README.md) |
| MEG | Preprocessing, localization, reconstruction, statistics, and control analyses | [README](MEG/README.md) |
| fMRS | Spectral processing and concentration aggregation | [README](fMRS/README.md) |
| figures | Plotting Fig.1–4 and Fig.S1–S8 | [README](figures/README.md) |
| utilities | Shared statistical and plotting functions | [README](utilities/README.md) |

For behavioral analyses, use the supplied 85-participant trial data, then fit, aggregate, and validate the behavioral models. The 19-participant MEG behavioral file is supplied separately. Behavioral and model drivers locate these files relative to their saved location and write new results under analysis_outputs; keep the package directories together. Model trajectories, mean-field simulations, MEG, and fMRS are separate analysis modules. Figure reproduction uses the supplied prepared results directly.

## Software dependencies

| Software | Version or purpose |
|---|---|
| MATLAB | R2024b |
| Statistics and Machine Learning Toolbox | Mixed models and statistical tests |
| Global Optimization Toolbox | Model fitting with particleswarm and patternsearch |
| Parallel Computing Toolbox | parfor and parallel simulations |
| Signal Processing Toolbox | Signal and spectral processing |
| Image Processing Toolbox | Upstream motor F>2 connected-component export (bwconncomp); not needed to replot supplied NIfTIs |
| FieldTrip | fieldtrip-20251218; MEG analysis |
| SPM12 | Structural imaging and tissue processing |
| VBA toolbox | V3.0; Bayesian model comparison |
| runica | ICA implementation bundled with FieldTrip |
| GANNET | 3.5.2; spectral processing and concentration fitting |
| Python | 3.9.18; mean-field simulations |
| NumPy | 1.26.4 |
| SciPy | 1.13.1 |
| Matplotlib | 3.8.3 |
| Seaborn | 0.13.2 |
| specparam | 2.0.0rc6 |
| tqdm | 4.67.1 |
| joblib | 1.2.0 |
| Jupyter Notebook | 6.4.8; running gamma.ipynb |

For the Python environment, install the pinned packages from the package root with `python -m pip install -r requirements.txt`. See [mean-field instructions](mean-field/README.md) for notebook execution.

## Preparation for manual execution

Start in the code directory and add the analysis modules and shared functions to the path. Add the FieldTrip root directory and call ft_defaults.

```matlab
codeRoot = pwd;  % Current directory: this project's code folder
addpath(fullfile(codeRoot,'utilities'));
addpath(fullfile(codeRoot,'behavior'));
addpath(fullfile(codeRoot,'behavioral_model'));
addpath(fullfile(codeRoot,'model_dynamics'));
addpath(fullfile(codeRoot,'MEG'));
addpath(fullfile(codeRoot,'fMRS'));
% For MEG analyses, add your local FieldTrip root and call ft_defaults.
% Add VBA for model comparison and SPM12 for structural processing as needed.
```

For upstream analysis scripts, set load/save paths and output directories to match the actual data locations. Figure scripts use the bundled relative paths described above. Behavioral_analysis.m, main_model_prediction.m, and neurotransmitter_stat.m write recomputed outputs to analysis_outputs, preserving the prepared figure inputs. Other MEG/MRS upstream scripts require external raw/intermediate inputs; SCRIPT_INDEX.md lists their purposes and configuration requirements. Run data-loading and parameter-setting sections before analysis sections. For scripts containing clear or clearvars, load variables according to their workspace preservation rules.

## Samples and conditions

Participants were excluded manually based on experimental behavior and data completeness, including sleepiness, early withdrawal, and disordered data; criteria and exclusions are described in Methods. The analysis sample includes 37 behavioral, 19 MEG, and 29 fMRS participants, totaling 85. Public IDs and cohort row indices are listed in figure_data/participant_index.csv. Cohort subsets retain the IDs from the complete sample. Use consistent participant IDs and ordering across the four conditions.

The condition order is `ots_data, pls_data, otn_data, pln_data`, corresponding to OT-social, PL-social, OT-asocial, and PL-asocial. See each module's instructions for the column order of model parameters and GLMM coefficients.

## External inputs

Raw MEG and MRS recordings, structural MRI, and upstream intermediate analysis files are not included. Figure reproduction uses the supplied figure_data and does not require these external inputs.

For MEG/MRS analyses using separately obtained acquisition data, set each relevant script's input and output paths to the local data locations. Behavioral model analyses use the supplied trial data and package-relative paths. Preserve the subject/session/run folder hierarchy and use a writable output directory. Recomputed behavioral, prediction, and concentration summaries are saved separately under analysis_outputs. Configure toolbox roots for the current MATLAB session.
