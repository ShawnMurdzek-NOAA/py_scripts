"""
Plot Skew-T, logp Diagrams With Subsets of Lines

shawn.s.murdzek@noaa.gov
Date Created: 3 June 2023
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import metpy.calc as mc
from metpy.units import units
from metpy.plots import SkewT, Hodograph


#---------------------------------------------------------------------------------------------------
# Create Plots
#---------------------------------------------------------------------------------------------------

# Read in sample data
fname = 'bmx_201104280000'
snd_df = pd.read_csv(fname, skiprows=[0, 1, 3, 4, 5], delim_whitespace=True)
snd_df['HGHT_AGL'] = snd_df['HGHT'] - snd_df['HGHT'].iloc[0]

# Function to plot Skew-T, logp
def plot_skewt(df, dry=True, moist=True, mix=True, data=True):
    """
    Plot a skew-T, logp with various toggles to turn certain aspects of the plot on or off
    """

    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)
    skew.ax.set_xlim([-40, 40])
    skew.ax.set_ylim([1000, 100])
    skew.ax.set_xlabel('temperature ($^{\circ}$C)', size=14)
    skew.ax.set_ylabel('pressure (mb)', size=14)

    if data:
        skew.plot(df['PRES'], snd_df['TEMP'], 'r', lw=2.5)
        skew.plot(df['PRES'], snd_df['DWPT'], 'g', lw=2.5)

    if dry:
        skew.plot_dry_adiabats()

    if moist:
        skew.plot_moist_adiabats()

    if mix:
        skew.plot_mixing_lines()

    plt.show()

    return None

#plot_skewt(snd_df, dry=False, moist=False, mix=False, data=False)
#plot_skewt(snd_df, dry=False, moist=False, mix=False, data=True)
#plot_skewt(snd_df, dry=True, moist=False, mix=False, data=False)
#plot_skewt(snd_df, dry=True, moist=True, mix=False, data=False)
#plot_skewt(snd_df, dry=True, moist=True, mix=True, data=False)
#plot_skewt(snd_df, dry=True, moist=True, mix=True, data=True)

# Plot hodograph
# Note that the last wind ob is removed b/c it is unreasonably large
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
h = Hodograph(ax, component_range=100)
h.add_grid(increment=20)

ind_1 = np.where(snd_df['HGHT_AGL'].values <= 1000)[0][-1]
ind_3 = np.where(snd_df['HGHT_AGL'].values <= 3000)[0][-1]
ind_6 = np.where(snd_df['HGHT_AGL'].values <= 6000)[0][-1]
ind_9 = np.where(snd_df['HGHT_AGL'].values <= 9000)[0][-1]

wdir = snd_df['DRCT'].values[:-1] * units.deg
wspd = snd_df['SKNT'].values[:-1] * units.kt
u, v = mc.wind_components(wspd, wdir)

h.plot(u.magnitude[:(ind_1+1)], v.magnitude[:(ind_1+1)], c='r')
h.plot(u.magnitude[ind_1:(ind_3+1)], v.magnitude[ind_1:(ind_3+1)], c='b')
h.plot(u.magnitude[ind_3:(ind_6+1)], v.magnitude[ind_3:(ind_6+1)], c='y')
h.plot(u.magnitude[ind_6:(ind_9+1)], v.magnitude[ind_6:(ind_9+1)], c='purple')
h.plot(u.magnitude[ind_9:], v.magnitude[ind_9:], c='k')

ax.set_xlabel('u (kts)', size=14)
ax.set_ylabel('v (kts)', size=14)

plt.show()


"""
End skewt_hodo_plots.py
"""
