# Difference Plots

The code in this directory is designed to create difference plots using UPP output from two different model runs (e.g., RRFS vs RRFS, or RRFS vs WRF).

## Contents

- `hcrsxn_diffs_input.yml`: Input YAML file. Controls which fields to examine, which times and model runs to use, and plotting parameters (e.g., filled contour intervals).
- `make_submit_diff_plots.sh`: Script to create several difference plot jobs for different time periods. Essentially a very simple way to parallelize the code.
- `plot_hcrsxn_diffs.py`: Main driver for difference plots. Relies on `plot_model_data` and `upp_postprocess` from [pyDA_utils](https://github.com/ShawnMurdzek-NOAA/pyDA_utils).
- `run_diff.sh`: Simple batch submission script used to run `plot_hcrsxn_diffs.py`.
