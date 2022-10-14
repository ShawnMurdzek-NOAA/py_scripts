"""
Class to Handle BUFR Output and Input CSVs

Output BUFR CSVs are created using the prepbufr_decode_csv.f90 utility in GSI-utils.

shawn.s.murdzek@noaa.gov
Date Created: 13 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------
 
import pandas as pd
import json
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


#---------------------------------------------------------------------------------------------------
# Define BUFR CSV Class
#---------------------------------------------------------------------------------------------------

class bufrCSV():
    """
    Class that handles CSV files that can easily be encoded to or decoded from BUFR files using 
    GSI-utils.

    Parameters
    ----------
    fname : string
        CSV file name

    """

    def __init__(self, fname):
   
        df = pd.read_csv(fname, usecols=list(range(35)), dtype={'SID':str})
  
        # Set missing values (1e11) to NaN
        self.df = df.where(df != 1e11)

        # Load metadata from JSON file
        self.meta = json.load(open('/mnt/lfs4/BMC/wrfruc/murdzek/py_scripts/modules/bufr_meta.json', 'r'))


#---------------------------------------------------------------------------------------------------
# Additional Functions
#---------------------------------------------------------------------------------------------------

def plot_obs_locs(x, y, fig=None, nrows=1, ncols=1, axnum=1, proj=ccrs.PlateCarree(), 
                  borders=True, **kwargs):
    """
    Plot observation locations using CartoPy.

    Parameters
    ----------
    x, y : float
        Observation locations (deg E, deg N)
    fig : matplotlib.pyplot.figure, optional
        Figure object to add axes to
    nrows, ncols, axnum : integer, optional
        Number of rows/columns of axes, and axes number for this plot
    proj : optional
        CartoPy map projection
    borders : boolean, optional
        Option to add country borders
    **kwargs : optional
        Other keyword arguments passed to matplotlib.pyplot.plot()

    Returns
    -------
    ax : matplotlib.axes
        Matplotlib.axes instance

    """   
  
    if fig == None:
        fig = plt.figure()

    ax = fig.add_subplot(nrows, ncols, axnum, projection=proj)

    # Add map features
    ax.coastlines('50m')
    if borders:
        ax.add_feature(cfeature.BORDERS)
    ax.gridlines()

    ax.plot(x, y, transform=proj, **kwargs)

    return ax


"""
End bufr.py
"""
