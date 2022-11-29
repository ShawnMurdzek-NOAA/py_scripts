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
import numpy as np


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

        # Remove space before nmsg
        self.df.rename(columns={' nmsg':'nmsg'}, inplace=True)

        # Load metadata from JSON file
        self.meta = json.load(open('/mnt/lfs4/BMC/wrfruc/murdzek/py_scripts/modules/bufr_meta.json', 
                                   'r'))

    def sample(self, fname):
        """
        Create a sample prepbufr CSV using only two lines from each unique prepbufr report type

        Parameters
        ----------
        fname : string
            Filename for sample bufr CSV file

        Returns
        -------
        None

        """

        # Extract the first two obs for each prepbufr report type
        typ = self.df['TYP'].unique()
        df = self.df.loc[self.df['TYP'] == typ[0]].iloc[:2]
        for t in typ[1:]:
            df = pd.concat([df, self.df.loc[self.df['TYP'] == t].iloc[:2]])

        # Replace NaNs with 1e11
        df = df.fillna(1e11)

        # Reorder message numbers
        nmsg = df['nmsg'].values
        for i, msg in enumerate(df['nmsg'].unique()):
            nmsg[np.where(nmsg == msg)] = i+1
        df['nmsg'] = nmsg        

        # Overwrite aircraft and ship SIDs
        sid = df['SID'].values
        sub = df['subset'].values
        for i, s in enumerate(sub):
            if s in ['AIRCAR', 'AIRCFT', 'SFCSHP']:
                sid[i] = 'SMPL%02d' % i
        df['SID'] = sid

        # Write to CSV
        df.to_csv(fname, index=False)
   
        # Add leading blank spaces and remove the .1 and .2 from the various NUL fields
        fptr = open(fname, 'r')
        contents = fptr.readlines()
        fptr.close()

        items = contents[0].split(',')
        for i, it in enumerate(items):
            if it[:3] == 'NUL':
                items[i] = 'NUL'

        contents[0] = ','.join(items)

        fptr = open(fname, 'w')
        for l in contents:
            fptr.write(' ' + l.strip() + ',\n')
        fptr.close() 


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


def df_to_csv(df, fname):
    """
    Write a DataFrame to a BUFR CSV file

    Parameters
    ----------
    df : pd.DataFrame
        Pandas DataFrame with the same format as a BUFR DataFrame
    fname : string
        Filename for bufr CSV file

    Returns
    -------
    None

    """

    # Replace NaNs with 1e11
    df = df.fillna(1e11)

    # Reorder message numbers
    nmsg = df['nmsg'].values
    for i, msg in enumerate(df['nmsg'].unique()):
        nmsg[np.where(nmsg == msg)] = i+1
    df['nmsg'] = nmsg        

    # Write to CSV
    df.to_csv(fname, index=False)
   
    # Add leading blank spaces and remove the .1 and .2 from the various NUL fields
    fptr = open(fname, 'r')
    contents = fptr.readlines()
    fptr.close()

    items = contents[0].split(',')
    for i, it in enumerate(items):
        if it[:3] == 'NUL':
            items[i] = 'NUL'

    contents[0] = ','.join(items)

    fptr = open(fname, 'w')
    for l in contents:
        fptr.write(' ' + l.strip() + ',\n')
    fptr.close() 


"""
End bufr.py
"""
