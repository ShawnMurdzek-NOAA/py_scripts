# Convolutions of 2D Difference Fields

**Purpose**: To compare large amounts of output from two forecasts to identify mesoscale regions where there are large differences in a specific field

### Approach

For each initialization time and forecast hour combination...

1. Compute the squared difference field betwen the two forecasts
2. Perform a 2D convolution using a user-defined kernel (i.e., a square matrix of ones)
3. Divide the resulting convolutions by the size of the kernel and take the square root to obtain root-mean-squared differences (RMSD) for each "patch" (small region with the same size as the kernel).
4. Save output as a CSV file, sorted by RMSDs from largest to smallest

Thus, the patches with the largest differences between the two forecasts will be at the top of the CSV. These should represent mesoscale events where the two forecasts greatly diverged.

### Assumptions

1. Both forecasts are on the same grid
2. Output is in Unified Post Processor (UPP) GRIB2 format

### Running the code

Program options are controlled using the `input_diff_2d_convolve.yml` file. To run, use the following command:

`python diff_2d_convolve.py input_diff_2d_convolve.yml`

To run several instances as batch jobs (using `run_diff.sh`)

1. Create a clean directory
2. Copy `make_submit_diff_2d_convolve.sh` and `input_diff_2d_convolve_TEMPLATE.yml` to the new directory
3. Edit `make_submit_diff_2d_convolve.sh` and `input_diff_2d_convolve_TEMPLATE.yml`
4. `bash make_submit_diff_2d_convolve.sh`
