# Mean-field network and spectral simulations

gamma.ipynb contains network parameters, trial simulations, Welch spectra, specparam fitting of periodic components, simulations at 40 inhibition levels, and plotting. Execute cells in dependency order; large parallel simulations may take considerable time.

## Environment

Use Python 3.9.18 and execute gamma.ipynb cell by cell in Jupyter Notebook 6.4.8.

| Software or package | Version | Purpose |
|---|---|---|
| Python | 3.9.18 | Runtime environment |
| NumPy | 1.26.4 | Arrays and network simulations |
| SciPy | 1.13.1 | Welch spectra, statistics, and MAT-file I/O |
| Matplotlib | 3.8.3 | Plotting |
| Seaborn | 0.13.2 | Statistical plotting |
| specparam | 2.0.0rc6 | SpectralModel spectral parameterization |
| tqdm | 4.67.1 | Progress display |
| joblib | 1.2.0 | Parallel simulations |
| Jupyter Notebook | 6.4.8 | Interactive notebook execution |

## Execution order and outputs

First run fixed-parameter and function definitions, then get_raw_power/get_peak_power and the aggregate-result function definitions. Run the two-condition demonstration cells and the formal 40-level simulation cells separately. The 40-level section uses joblib Parallel(n_jobs=8), with 100 valid trials per level. The save cell outputs 40_level_data.mat, containing data, inhibition_strength, f, t_wind_start, window_length, and window_step. Later cells read this file and calculate gamma power and differences between strong and weak inhibition. The prepared article plotting input is supplied as ../figure_data/MEG_results/40_level_data.mat; the MATLAB figure scripts locate it automatically. Running the notebook is not required to reproduce these prepared plots.

## Frequency configuration

The notebook uses freq_range=(10,150) for specparam fitting and peak_width_limits=[20,999]. The supplied 40_level_data.mat frequency bins span approximately 10.005-145.073 Hz. The article uses the 10-120 Hz subset for downstream display/analysis, or narrower target bands selected explicitly (e.g., 40-70 Hz for the notebook gamma summary). Selecting saved frequency bins or limiting plot axes does not refit specparam over 10-120 Hz.

## Random pairing

Strong- and weak-inhibition groups are selected from the high and low groups and paired randomly. Indices from one random pairing are supplied for figure reproduction. Fig.3 uses pair_data.mat, and Fig.S6 uses order_idx1/order_idx2 in pair_data_low.mat. Load the saved indices when replotting rather than generating a new pairing. Pairing files are prerequisite plotting inputs. If generating another pairing, save both groups' indices and their corresponding inhibition levels, distinguishing MATLAB's one-based indices from Python's zero-based indices.

See [the complete file index](../SCRIPT_INDEX.md) for script purposes and input requirements.
