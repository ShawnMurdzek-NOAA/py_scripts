"""
Plot the Reduction in the Number of UAS Obs Owing to Wind Speed and Icing Limits

Raw Obs (Used Here)
-------------------
Use output from log files to determine how many obs are lost by applying meteorological limits


Superobs (Not Used Here)
------------------------
Counting the number of UAS superobs can be done using ECMWF's eccodes library, which is available in
my pygrib environment. The following code snippet can be used to count all obs in a BUFR file:

fname = '/work/noaa/wrfruc/murdzek/nature_run_spring/obs/uas_obs_300km/syn_bufr/2022042922.rap.t22z.prepbufr.tm00'
out_str = subprocess.check_output(f"bufr_ls -w typicalDate=20220429 -p numberOfSubsets {fname}", shell=True).decode('utf-8')
nobs = np.sum(np.array([int(i) for i in out_str.split('\n')[2:-4]]))


shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import os
import glob
import subprocess
import copy
import matplotlib.pyplot as plt


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

# Template for UAS ob directories
# Include placeholders for {season} and {exp}
dir_tmpl = '/work/noaa/wrfruc/murdzek/nature_run_{season}/obs/uas_obs_35km_{exp}'

# Seasons
all_seasons = ['spring', 'winter']

# Experiments
all_exp = ['icing', 'wspd20', 'wspd30', 'wspd40', 'icing_wspd20', 'icing_wspd30', 'icing_wspd40']

# All possible wind speeds
all_wspd = [20, 30, 40]

# Output plot file name
out_fname = 'UAS_count_reduction.png'


#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

ob_counts = {}
for season in all_seasons:
    ob_counts[season] = {}
    for exp in all_exp:
        print(f"Counting obs for {season} {exp}")

        # Untar and unzip log files
        d = dir_tmpl.format(season=season, exp=exp)
        os.system(f"tar xzf {d}/log.tar.gz")
        os.system("mkdir tmp && mv ./log ./tmp/")

        # Loop over each file
        files = glob.glob(f"./tmp/log/*limit_uas.log")
        initial = 0
        final = 0
        for f in files:
            try:
                out = subprocess.check_output(f"grep '136 =' {f}", shell=True).decode('utf-8').split('\n')
            except subprocess.CalledProcessError:
                print(f"Cannot perform grep on {f}. Skipping")
                continue
            if len(out) != 3:
                raise IOError(f"file {f} has the wrong number of grep matches")
            initial = initial + float(out[0].split('=')[1])
            final = final + float(out[1].split('=')[1].split('(')[0])

        # Update ob counts
        ob_counts[season][exp] = {'initial':initial, 
                                  'final':final,
                                  'pct_reduce': 1e2 * (final - initial) / initial}

        # Delete log file directory
        os.system(f"rm -r ./tmp")

# Create plot
print('\nMaking plot')
fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True, figsize=(6, 8))
plt.subplots_adjust(left=0.12, bottom=0.12, right=0.98, top=0.93)
for season, ax in zip(all_seasons, axes):

    for tag, label, c in zip(['icing_', ''], ['icing limit', 'no icing limit'], ['b', 'r']):
        xdata = copy.deepcopy(all_wspd) + [all_wspd[-1] + 10]
        ydata = []
        for wspd in xdata[:-1]:
            ydata.append(ob_counts[season][f"{tag}wspd{wspd}"]["pct_reduce"])
        if tag == 'icing_':
            ydata.append(ob_counts[season]["icing"]["pct_reduce"])
        else:
            ydata.append(0)

        ax.plot(xdata, ydata, c=c, ls='-', marker='o', label=label)

    # Add plot annotations
    ax.legend(fontsize=14)
    ax.grid()
    xlabels = [str(n) for n in all_wspd] + ['None']
    ax.set_xticks(all_wspd + [all_wspd[-1] + 10], labels=xlabels)
    ax.set_title(season, size=20)
    ax.set_ylabel('Percent Reduction in Ob Counts (%)', size=12)
    if season == all_seasons[-1]:
        ax.set_xlabel('Wind Speed Limit (m/s)', size=14)

plt.savefig(out_fname)


"""
End plot_reduction_number_UAS_obs_wspd_icing.py
"""
