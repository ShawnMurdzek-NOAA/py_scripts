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
import scipy.interpolate as si

import meteo_util as mu


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

    def sample(self, fname, n=2):
        """
        Create a sample prepbufr CSV using only n lines from each unique prepbufr report type

        Parameters
        ----------
        fname : string
            Filename for sample bufr CSV file
        n : integer, optional
            Number of lines to use from each unique prepbufr report type

        Returns
        -------
        None

        """

        # Extract the first two obs for each prepbufr report type
        typ = self.df['TYP'].unique()
        df = self.df.loc[self.df['TYP'] == typ[0]].iloc[:n]
        for t in typ[1:]:
            df = pd.concat([df, self.df.loc[self.df['TYP'] == t].iloc[:n]])

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


def read_ob_errors(fname):
    """
    Parse out observation errors from an errtable file in GSI

    Parameters
    ----------
    fname : string
        Name of errtable text file

    Returns
    -------
    errors: dictionary
        A dictionary of pd.DataFrame objects containing the observation errors

    Notes
    -----
    More information about the errtable format in GSI can be found here: 
    https://dtcenter.ucar.edu/com-GSI/users/docs/users_guide/html_v3.7/gsi_ch4.html#conventional-observation-errors

    """

    # Extract contents of file
    fptr = open(fname, 'r')
    contents = fptr.readlines()
    fptr.close()

    # Loop over each line
    errors = {}
    headers = ['prs', 'Terr', 'RHerr', 'UVerr', 'PSerr', 'PWerr']
    for l in contents:
        if l[5:21] == 'OBSERVATION TYPE':
            key = int(l[1:4])
            errors[key] = {}
            for h in headers:
                errors[key][h] = []
        else:
            vals = l.strip().split(' ')
            for k, h in enumerate(headers):
                errors[key][h].append(float(vals[k]))

    # Convert to DataFrame
    for key in errors.keys():
        errors[key] = pd.DataFrame(errors[key])
        for h in headers[1:]:
            errors[key][h].where(errors[key][h] < 5e8, inplace=True)

    return errors


def add_obs_err_uncorr(df, errtable, typ='all'):
    """
    Add random, uncorrelated Gaussian errors to observations based on error standard deviations in 
    errtable

    Parameters
    ----------
    df : pd.DataFrame
        Pandas DataFrame with the same format as a BUFR DataFrame
    errtable : string
        Name of errtable text file
    typ: 'all' or list, optional
        List of observation types to apply errors to 

    Returns
    -------
    out_df : pd.DataFrame
        DataFrame of observations with random errors added

    Notes
    -----
    Observation errors are assumed to be uncorrelated in both space and time

    """
    
    # Make copy of DataFrame
    out_df = df.copy()

    # Determine list of all observation types
    if typ == 'all':
        typ = np.int32(out_df['TYP'].unique())

    # Read in error table file 
    etable = read_ob_errors(errtable)
    eprs = etable[100]['prs'].values

    # Convert specific humidities to relative humidities in BUFR CSV
    out_df['QOB'] = out_df['QOB'] * 1e-6
    mix = out_df['QOB']  / (1. - out_df['QOB'])
    out_df['RHOB'] = 10 * (mix / mu.equil_mix(out_df['TOB'] + 273.15, out_df['POB'] * 1e2))  

    # Convert surface pressures from Pa to hPa in BUFR CSV
    out_df['PRSS'] = out_df['PRSS'] * 1e-2

    # Loop over each observation type
    for t in typ:
        print()
        print('TYP = %d' % t)
        for ob, err in zip(['TOB', 'RHOB', 'UOB', 'VOB', 'PRSS', 'PWO'],
                           ['Terr', 'RHerr', 'UVerr', 'UVerr', 'PSerr', 'PWerr']):
            
            print(ob)
            # Check to see if errors are defined
            if np.all(np.isnan(etable[t][err].values)):
                print('no errors for typ = %d, ob = %s' % (t, ob))
                continue

            # Determine indices where errors are not NaN
            eind = np.where(np.logical_not(np.isnan(etable[t][err])))[0]
            
            # Determine indices in out_df for this ob type
            oind = np.where((out_df['TYP'] == t) & (out_df['POB'] <= eprs[eind[0]]) & 
                            (out_df['POB'] >= eprs[eind[-1]]))[0]

            # Determine if errors vary with pressure. If so, create a function for interpolation
            # Then add observation errors
            if len(np.unique(etable[t].loc[eind, err])) == 1:
                error = etable[t].loc[eind[0], err]
                out_df.loc[oind, ob] = (out_df.loc[oind, ob] +
                                        np.random.normal(scale=error, size=len(oind)))
            else:
                fct = si.interp1d(eprs[eind], etable[t].loc[eind, err])
                for k in oind:
                    error = fct(out_df.loc[k, 'POB'])
                    out_df.loc[k, ob] = out_df.loc[k, ob] + np.random.normal(scale=error)

    # Set relative humidities > 100% to 100%
    out_df.loc[out_df['RHOB'] > 10, 'RHOB'] = 10.

    # Compute specific humidity from relative humidity
    mix = 0.1 * out_df['RHOB'] * mu.equil_mix(out_df['TOB'] + 273.15, out_df['POB'] * 1e2)
    out_df['QOB'] = 1e6 * (mix / (1. + mix))
    out_df.drop(label='RHOB', axis=1, inplace=True)    
    
    # Convert surface pressure back to Pa
    out_df['PRSS'] = out_df['PRSS'] * 1e-2
    
    return out_df


"""
End bufr.py
"""
