"""
Extract an Area-Averaged Model Sounding

shawn.s.murdzek@noaa.gov
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import datetime as dt
import sys
import argparse
import copy
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.calc as mc
from metpy.units import units
from metpy.plots import SkewT, Hodograph



#---------------------------------------------------------------------------------------------------
# Main Program
#---------------------------------------------------------------------------------------------------

def parse_in_args(argv):
    """
    Parse input arguments

    Parameters
    ----------
    argv : list
        Command-line arguments from sys.argv[1:]

    Returns
    -------
    Parsed input arguments

    """

    parser = argparse.ArgumentParser(description='Extract an area-averaged sounding from a CM1 or \
                                                  WRF output file.')

    # Positional arguments
    parser.add_argument('in_file',
                        help='Model output netCDF file (e.g., cm1out_rst or wrfout file. \
                              Note that only CM1 restart files are supported)',
                        type=str)

    parser.add_argument('out_txt_file',
                        help='Output sounding text file',
                        type=str)

    parser.add_argument('out_image_file',
                        help='Output file with sounding image',
                        type=str)

    # Optional arguments
    parser.add_argument('-m',
                        dest='model',
                        default='CM1',
                        help='Model. Options are "CM1" or "WRF"',
                        type=str)

    parser.add_argument('-s',
                        dest='snd_fname',
                        default='./input_sounding',
                        help='Base-state CM1 sounding',
                        type=str)

    parser.add_argument('--istart',
                        dest='istart',
                        default=0,
                        help='Starting i index for averaging domain',
                        type=int)

    parser.add_argument('--iend',
                        dest='iend',
                        default=-1,
                        help='Ending i index for averaging domain',
                        type=int)

    parser.add_argument('--jstart',
                        dest='jstart',
                        default=0,
                        help='Starting j index for averaging domain',
                        type=int)

    parser.add_argument('--jend',
                        dest='jend',
                        default=-1,
                        help='Ending j index for averaging domain',
                        type=int)

    return parser.parse_args(argv)


def read_cm1_input_sounding(fname):
    """
    Read in CM1 base state from input_sounding
    """

    return pd.read_csv(fname,
                       sep='\s+',
                       names=['hgt', 'theta', 'qv', 'u', 'v'],
                       skiprows=1)


def read_cm1_rst(param):
    """
    Read CM1 restart file
    """

    out_dict = {}

    # Read in base-state sounding and restart file
    base_df = read_cm1_input_sounding(param.snd_fname)
    ds = xr.open_dataset(param.in_file, decode_timedelta=False)

    # Extract height (m AGL)
    out_dict['hgt'] = ds['zh'].values[:, np.newaxis, np.newaxis]

    # Extract theta perturbations and add base-state theta
    th_base = np.interp(np.squeeze(out_dict['hgt']), base_df['hgt'].values, base_df['theta'].values)
    out_dict['theta'] = (ds['tha'][0, :, param.jstart:param.jend, param.istart:param.iend].values +
                         th_base[:, np.newaxis, np.newaxis])

    # Extract water vapor mixing ratio and convert from kg/kg to g/kg
    out_dict['qv'] = 1e3 * ds['qv'][0, :, param.jstart:param.jend, param.istart:param.iend].values

    # Extract winds and interpolate to mass grid
    u = 0.5 * (ds['ua'][0, :, :, :-1].values + ds['ua'][0, :, :, 1:].values)
    v = 0.5 * (ds['va'][0, :, :-1, :].values + ds['va'][0, :, 1:, :].values)
    out_dict['u'] = u[:, param.jstart:param.jend, param.istart:param.iend]
    out_dict['v'] = v[:, param.jstart:param.jend, param.istart:param.iend]

    # Extract pressure and convert to hPa
    out_dict['p'] = 1e-2 * ds['prs'][0, :, param.jstart:param.jend, param.istart:param.iend].values

    # Extract surface conditions. Use 2-m theta and qv as surface values
    out_dict['psfc'] = 1e-2 * ds['psfc'][0, param.jstart:param.jend, param.istart:param.iend].values
    out_dict['thsfc'] = ds['th2'][0, param.jstart:param.jend, param.istart:param.iend].values
    out_dict['qvsfc'] = 1e3 * ds['q2'][0, param.jstart:param.jend, param.istart:param.iend].values

    return out_dict


def read_wrfout(param):
    """
    Read WRF output file
    """

    out_dict = {}
    ds = xr.open_dataset(param.in_file, decode_timedelta=False)

    # Extract height (m AGL)
    geopotential = (ds['PH'][0, :, param.jstart:param.jend, param.istart:param.iend].values +
                    ds['PHB'][0, :, param.jstart:param.jend, param.istart:param.iend].values)
    hgt = mc.geopotential_to_height(geopotential * units.m**2 / units.s**2).to(units.m).magnitude
    out_dict['hgt'] = 0.5 * (hgt[:-1] + hgt[1:])

    # Extract theta perturbations and add base-state value (300 K)
    out_dict['theta'] = ds['T'][0, :, param.jstart:param.jend, param.istart:param.iend].values + 300.

    # Extract water vapor mixing ratio and convert from kg/kg to g/kg
    out_dict['qv'] = 1e3 * ds['QVAPOR'][0, :, param.jstart:param.jend, param.istart:param.iend].values

    # Extract winds and interpolate to mass grid
    u = 0.5 * (ds['U'][0, :, :, :-1].values + ds['U'][0, :, :, 1:].values)
    v = 0.5 * (ds['V'][0, :, :-1, :].values + ds['V'][0, :, 1:, :].values)
    out_dict['u'] = u[:, param.jstart:param.jend, param.istart:param.iend]
    out_dict['v'] = v[:, param.jstart:param.jend, param.istart:param.iend]

    # Extract pressure
    out_dict['p'] = 1e-2 * (ds['P'][0, :, param.jstart:param.jend, param.istart:param.iend].values +
                            ds['PB'][0, :, param.jstart:param.jend, param.istart:param.iend].values)

    # Extract surface conditions. Use 2-m theta and qv as surface values
    out_dict['psfc'] = 1e-2 * ds['PSFC'][0, param.jstart:param.jend, param.istart:param.iend].values
    out_dict['thsfc'] = ds['TH2'][0, param.jstart:param.jend, param.istart:param.iend].values
    out_dict['qvsfc'] = 1e3 * ds['Q2'][0, param.jstart:param.jend, param.istart:param.iend].values

    return out_dict


def extract_avg_sounding(param, verbose=1):
    """
    Extract area-averaged sounding from model output
    """

    # Extract required 3D model fields
    if verbose > 0: print('extracting model data')
    if param.model == 'CM1':
        model_dict = read_cm1_rst(param)
    elif param.model == 'WRF':
        model_dict = read_wrfout(param)
    else:
        raise ValueError(f"model value {model} is not supported")

    # Compute area-averaged sounding
    if verbose > 0: print('computing area-averaged sounding')
    fields3d = ['hgt', 'theta', 'qv', 'u', 'v', 'p']
    snd_dict = {}
    for f in fields3d:
        snd_dict[f] = np.mean(model_dict[f], axis=(1, 2))

    # Compute area-averaged surface conditions
    fields2d = ['psfc', 'thsfc', 'qvsfc']
    sfc_dict = {}
    for f in fields2d:
        sfc_dict[f] = np.mean(model_dict[f], axis=(0, 1))

    return pd.DataFrame(snd_dict), sfc_dict


def write_sounding(snd_df, sfc_dict, fname):
    """
    Write sounding to an output text file
    """

    with open(fname, 'w') as fptr:

        # Write surface conditions
        fptr.write(f"{sfc_dict['psfc']:9.2f}{sfc_dict['thsfc']:11.3f}{sfc_dict['qvsfc']:11.3f}\n")

        # Write above-ground profile
        for i in range(len(snd_df)):
            fptr.write(f"{snd_df.loc[i,'hgt']:9.2f}{snd_df.loc[i,'theta']:11.3f}{snd_df.loc[i,'qv']:11.3f}{snd_df.loc[i,'u']:11.4f}{snd_df.loc[i,'v']:11.4f}\n")


def plot_sounding(snd_df):
    """
    Make a skew-T, logp plot of the sounding
    """

    # Compute T and Td
    p_units = snd_df['p'].values * units.hPa
    th_units = snd_df['theta'].values * units.K
    qv_units = snd_df['qv'].values * units.g / units.kg
    T = mc.temperature_from_potential_temperature(p_units, th_units).to('degC')
    RH = mc.relative_humidity_from_mixing_ratio(p_units, T, qv_units)
    Td = mc.dewpoint_from_relative_humidity(T, RH).to('degC')
    p = p_units.magnitude
    T = T.magnitude
    Td = Td.magnitude

    # Create figure
    fig = plt.figure(figsize=(8, 8))
    skew = SkewT(fig, rotation=45)
    skew.ax.set_xlim(-30, 60)
    skew.ax.set_ylim(1000, 100)

    # Plot thermodynamic profile
    skew.plot(p, T, 'r-', lw=3)
    skew.plot(p, Td, 'b-', lw=3)
    skew.plot_dry_adiabats(linewidth=0.75)
    skew.plot_moist_adiabats(linewidth=0.75)
    skew.plot_mixing_lines(linewidth=0.75)
    
    # Add hodograph
    hod_ax = inset_axes(skew.ax, '35%', '35%', loc=1)
    h = Hodograph(hod_ax, component_range=40)
    h.add_grid(increment=10)
    h.plot(snd_df['u'].values, snd_df['v'].values, c='b', ls='-')

    # Labels
    skew.ax.set_xlabel(r'temperature ($^{\circ}$C)', size=14)
    skew.ax.set_ylabel('pressure (hPa)', size=14)

    return fig


if __name__ == '__main__':

    start = dt.datetime.now()
    print('Starting extract_area_avg_model_sounding.py')
    print(f"Time = {start.strftime('%Y%m%d %H:%M:%S')}")

    # Parse arguments
    param = parse_in_args(sys.argv[1:])

    # Extract sounding
    snd_df, sfc_dict = extract_avg_sounding(param, verbose=1)

    # Write sounding to text file
    print('Writing sounding text file')
    write_sounding(snd_df, sfc_dict, param.out_txt_file)

    # Plot sounding
    print('Creating skew-T, logp diagram')
    fig = plot_sounding(snd_df)
    plt.savefig(param.out_image_file)

    print('\nProgram finished!')
    print(f"Elapsed time = {(dt.datetime.now() - start).total_seconds()} s")


"""
End extract_area_avg_model_sounding.py
"""
