"""
Add Observation Errors to a BUFR CSV File

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
out_df = bufr.add_obs_err_uncorr(bufr_csv.df, errtable)

end = dt.datetime.now()
print('elapsed time = %s s' % (end - start).total_seconds()) 


#---------------------------------------------------------------------------------------------------
# Plot Histograms of Differences
#---------------------------------------------------------------------------------------------------

if plot_diff_hist:

    etable = bufr.read_ob_errors(errtable)
    eprs = etable[100]['prs'].values
    obs = ['TOB', 'QOB', 'UOB', 'VOB', 'PRSS', 'PWO']
    errs = ['Terr', 'RHerr', 'UVerr', 'UVerr', 'PSerr', 'PWerr']
    typ = out_df['TYP'].unique()

    for t in typ:

        fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(10, 8))
        tind = np.where(out_df['TYP'] == t)[0]

        for j, (o, e) in enumerate(zip(obs, errs)):
            ax = axes[int(j/3), j%3]

            # Compute difference before and after adding obs errors
            diff = out_df.loc[tind, o].values - bufr_csv.df.loc[tind, o].values
            prs = out_df.loc[tind, 'POB'].values
            idiff = np.where(np.logical_not(np.isnan(diff)))[0]
            if len(idiff) > 0:
                diff = diff[idiff]
                prs = prs[idiff]
            else:
                continue

            if len(np.unique(etable[t][e])) <= 2:
                # Plot regular histogram
                ax.hist(diff)
                ax.set_title('%s (std = %.3e, estd = %.3e)' % (o, np.std(diff), etable[t].loc[0, e]))
                ax.set_xlabel('%s (%s)' % (o, bufr_csv.meta[o]['units']), size=12)
                ax.set_ylabel('counts', size=12)
                ax.grid()
            else:
                # Plot histogram that depends on pressure
                cts = np.zeros([21, len(eprs)-1])
                maxval = np.amax(np.abs(diff))
                bins = np.linspace(-maxval - (0.25*maxval), maxval + (0.25*maxval), 22)
                for k in range(len(eprs)-1):
                    idx = np.where(np.logical_and(prs <= eprs[k], prs > eprs[k+1]))[0]
                    if len(idx) == 0:
                        continue  
                    cts[:, k] = np.histogram(diff[idx], bins=bins)
                x2d, y2d = np.meshgrid((bins[1:] + bins[:-1]) / 2., (eprs[1:] + eprs[:-1]) / 2.)     
                cax = ax.contourf(x2d, y2d, cts)
                cbar = plt.colormap(cax, ax=ax)
                cbar.set_label('counts') 
                ax.set_title(o)
                ax.set_xlabel('%s (%s)' % (o, bufr_csv.meta[o]['units']), size=12)
                ax.set_ylabel('pressure (%s)' % bufr_csv['POB']['units'], size=12)
                ax.grid()

        plt.suptitle('Type = %d Difference Histograms' % t, size=16)
        plt.savefig('diff_hist_%d.png' % t)
        plt.close()


"""
End add_obs_errors.py 
"""
