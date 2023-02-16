"""
Add Observation Errors to a BUFR CSV File

Also include some code to plot histograms of the difference between the perfect obs and obs with 
added errors. These plots provide a sanity check to ensure that the distribution of the differences
is approximately Gaussian with a standard deviation close to the prescribed errors.

shawn.s.murdzek@noaa.gov
Date Created: 15 February 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
import meteo_util as mu

import bufr


#---------------------------------------------------------------------------------------------------
# Input Parameters
#---------------------------------------------------------------------------------------------------

bufr_fname = '/mnt/lfs4/BMC/wrfruc/murdzek/sample_real_obs/obs_rap/202204291200.rap.prepbufr.csv'
errtable = '/mnt/lfs4/BMC/wrfruc/murdzek/sample_real_obs/errtable.rrfs'

out_fname = './test.csv'

# Option to check obs errors by plotting differences between obs w/ and w/out errors
plot_diff_hist = True


#---------------------------------------------------------------------------------------------------
# Add Observation Errors
#---------------------------------------------------------------------------------------------------

start = dt.datetime.now()
print('start time = %s' % start.strftime('%H:%M:%S'))

bufr_csv = bufr.bufrCSV(bufr_fname)
out_df = pd.read_csv('out_bufr_errs.csv')
#out_df = bufr.add_obs_err_uncorr(bufr_csv.df, errtable)

end = dt.datetime.now()
print('elapsed time = %s s' % (end - start).total_seconds()) 


#---------------------------------------------------------------------------------------------------
# Plot Histograms of Differences
#---------------------------------------------------------------------------------------------------

if plot_diff_hist:

    # Compute RH
    q = bufr_csv.df['QOB'] * 1e-6
    mix = q  / (1. - q)
    bufr_csv.df['RHOB'] = 100 * (mix / mu.equil_mix(bufr_csv.df['TOB'] + 273.15, 
                                                    bufr_csv.df['POB'] * 1e2))

    q = out_df['QOB'] * 1e-6
    mix = q  / (1. - q)
    out_df['RHOB'] = 100 * (mix / mu.equil_mix(out_df['TOB'] + 273.15, out_df['POB'] * 1e2))

    # Add RH units to metadata
    bufr_csv.meta['RHOB'] = {}
    bufr_csv.meta['RHOB']['units'] = '%'

    # Read in error tables
    etable = bufr.read_ob_errors(errtable)
    eprs = etable[100]['prs'].values
    for e in etable.keys():
        etable[e]['RHerr'] = etable[e]['RHerr'] * 10
    typ = out_df['TYP'].unique()

    for t in typ:

        if t < 200:
            obs = ['TOB', 'RHOB']
            errs = ['Terr', 'RHerr']
        else:
            obs = ['UOB', 'VOB']
            errs = ['UVerr', 'UVerr']

        # Check to see if errors vary with pressure
        vprof = False
        for e in errs:
            if len(np.unique(etable[t][e])) > 2:
                vprof = True
                continue

        # Configure plot
        if vprof:
            fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 8), sharey=True)
            plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.92, hspace=0.35)
        else:
            fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))
            plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.9)

        tind = np.where(out_df['TYP'] == t)[0]

        for j, (o, e) in enumerate(zip(obs, errs)):

            # Compute difference before and after adding obs errors
            diff = out_df.loc[tind, o].values - bufr_csv.df.loc[tind, o].values
            prs = out_df.loc[tind, 'POB'].values
            idiff = np.where(np.logical_not(np.isnan(diff)))[0]
            if len(idiff) > 0:
                diff = diff[idiff]
                prs = prs[idiff]
            else:
                continue

            if vprof:

                # Create histogram that depends on pressure as well as vertical profiles of standard
                # deviations
                cts = np.zeros([len(eprs)-1, 21])
                stdev = np.zeros(len(eprs)-1)
                num = np.zeros(len(eprs)-1)
                maxval = np.amax(np.abs(diff))
                bins = np.linspace(-maxval - (0.25*maxval), maxval + (0.25*maxval), 22)
                eprs_plt = (eprs[1:] + eprs[:-1]) / 2.
                for k in range(len(eprs)-1):
                    idx = np.where(np.logical_and(prs <= eprs[k], prs > eprs[k+1]))[0]
                    if len(idx) == 0:
                        continue  
                    cts[k, :] = np.histogram(diff[idx], bins=bins)[0]
                    cts[k, np.isclose(cts[k, :], 0)] = np.nan
                    stdev[k] = np.std(diff[idx])
                    num[k] = len(idx)

                # Top plot: Histogram
                ax1 = axes[0, j]     
                x2d, y2d = np.meshgrid((bins[1:] + bins[:-1]) / 2., eprs_plt)
                cax = ax1.contourf(x2d, y2d, cts)
                cbar = plt.colorbar(cax, ax=ax1)
                cbar.set_label('counts') 
                ax1.set_xlabel('%s (%s)' % (o, bufr_csv.meta[o]['units']), size=12)
                if j == 0:  
                    ax1.set_ylabel('pressure (%s)' % bufr_csv.meta['POB']['units'], size=12)
                ax1.set_ylim([eprs.max(), eprs.min()])
                ax1.grid()

                # Bottom plot: Vertical profile of standard deviations
                ax2 = axes[1, j]
                ax21 = ax2.twiny()
                ax2.plot(stdev, eprs_plt, 'b-', label='actual')
                ax2.plot(etable[t][e], eprs, 'r-', label='errtable')
                ax21.plot(num, eprs_plt, 'k--', label='counts')
                ax2.set_xlabel('std %s (%s)' % (o, bufr_csv.meta[o]['units']), size=12)
                ax21.set_xlabel('data counts', size=12)
                if j == 0:  
                    ax2.set_ylabel('pressure (%s)' % bufr_csv.meta['POB']['units'], size=12)
                ax2.set_ylim([eprs.max(), eprs.min()])
                ax2.legend()
                ax2.grid()

            else:

                # Plot regular histogram
                ax = axes[j]
                ax.hist(diff)
                ax.set_title('%s (std = %.3e, estd = %.3e)' % (o, np.std(diff), etable[t].loc[0, e]))
                ax.set_xlabel('%s (%s)' % (o, bufr_csv.meta[o]['units']), size=12)
                ax.set_ylabel('counts', size=12)
                ax.grid()

        plt.suptitle('Type = %d Difference Histograms' % t, size=16)
        plt.savefig('diff_hist_%d.png' % t)
        plt.close()


"""
End add_obs_errors.py 
"""
