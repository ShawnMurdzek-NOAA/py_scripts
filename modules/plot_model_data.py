"""
Real-Data Model Plotting Class

CartoPy Note: CartoPy downloads shapefiles to plot coastlines, state borders, etc. The package used
to do this does not work on Jet, so these files must be downloaded manually (the only way I've found
to do this is to download them onto my Mac and upload them to Jet via scp). These shapefiles can
be downloaded from the following website:

https://www.naturalearthdata.com/downloads/

and must be uploaded to the following directory, then unzipped::

~/.local/share/cartopy/shapefiles/natural_earth/<category>/

where <category> is specified in cfeature.NaturalEarthFeature (usually physical or cultural). 

Shawn Murdzek
shawn.s.murdzek@noaa.gov
Date Created: 4 October 2022
"""

#---------------------------------------------------------------------------------------------------
# Import Modules
#---------------------------------------------------------------------------------------------------

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import xarray as xr
from metpy.plots import SkewT, Hodograph

try:
    import wrf
except ImportError:
    print('cannot load WRF-python module')


#---------------------------------------------------------------------------------------------------
# Define Plotting Class
#---------------------------------------------------------------------------------------------------

class PlotOutput():
    """
    Class that handles plotting real-data model output from either WRF or FV3 as well as gridded
    observational data (e.g., Stage IV precip).

    Parameters
    ----------
    data : List of NetCDF4 filepointers (for WRF) or XArray DataSets (all other options)
        Actual data to plot
    dataset : string
        Dataset type ('wrf', 'fv3', 'upp', or 'stage4')
    fig : matplotlib.figure
        Figure that axes are added to
    nrows : integer
        Number of rows of subplots
    ncols : integer
        Number of columns of subplots
    axnum : integer
        Subplot number
    proj : CartoPy projection, optional
        Map projection

    Notes
    -----

    'wrf' datasets can be opened using the following command with netCDF4:
    `nc.Dataset(filename)`

    'upp' and 'stage4' datasets can be opened using the following xarray command:
    `xr.open_dataset(filename, engine='pynio')`

    """

    def __init__(self, data, dataset, fig, nrows, ncols, axnum, proj=ccrs.PlateCarree()):

        self.outtype = dataset
        self.fig = fig
        self.nrows = nrows
        self.ncols = ncols
        self.n = axnum
        self.proj = proj

        # Dictionary used to hold metadata
        self.metadata = {}

        # Extract first dataset
        if self.outtype == 'wrf':
            self.fptr = data[0]
        elif self.outtype == 'fv3':
            raise ValueError('Raw FV3 output is not supported yet')
        elif (self.outtype == 'upp' or self.outtype == 'stage4'):
            self.ds = data[0]

        # Extract second dataset (if applicable)
        if len(data) > 1:
            if self.outtype == 'wrf':
                self.fptr2 = data[1]
            elif self.outtype == 'fv3':
                raise ValueError('Raw FV3 output is not supported yet')
            elif (self.outtype == 'upp' or self.outtype == 'stage4'):
                self.ds2 = data[1]


        # Extract time
        if (self.outtype == 'stage4' or self.outtype == 'upp'):
            sample = list(self.ds.keys())[0]
            itime = dt.datetime.strptime(self.ds[sample].attrs['initial_time'], '%m/%d/%Y (%H:%M)')
            if self.ds[sample].attrs['forecast_time_units'] == 'hours':
                delta = dt.timedelta(hours=int(self.ds[sample].attrs['forecast_time'][0]))
            elif self.ds[sample].attrs['forecast_time_units'] == 'minutes':
                delta = dt.timedelta(minutes=int(self.ds[sample].attrs['forecast_time'][0]))
            elif self.ds[sample].attrs['forecast_time_units'] == 'days':
                delta = dt.timedelta(days=int(self.ds[sample].attrs['forecast_time'][0]))
            self.time = (itime + delta).strftime('%Y%m%d %H:%M:%S UTC')

            
    def _ingest_data(self, var, zind=np.nan, units=None, interp_field=None, interp_lvl=None, 
                     ptype='none0', diff=False):
        """
        Extract a single variable to plot and interpolate if needed.

        Parameters
        ----------
        var : string
            Variable from model output file to plot
        zind : integer, optional
            Index in z-direction
        units : string, optional
            Units for var (WRF only)
        interp_field : string, optional
            Interpolate var to a surface with a constant value of interp_field (WRF only)
        interp_lvl : string, optional
            Value of constant-interp_field used during interpolaton (WRF only)
        ptype : string, optional
            Plot type. Used as the key to store metadata
        diff : boolean, optional
            Is this a difference plot?

        Returns
        -------
        data : array
            An array containing var
        coords : list
            List of coordinates for data.
                For a horizontal cross section, this list is [lats, lons]
        ptype : string
            Plot type. Used as the key to store metadata

        """

        coords = []

        # Append an integer to ptype if already in use
        n = 1
        while ptype in self.metadata.keys():
            ptype = ptype[:-1] + str(n)
            n = n + 1
        self.metadata[ptype] = {}

        if self.outtype == 'wrf':

            # Extract raw variable from WRF output file
            if units != None:
                raw = wrf.getvar(self.fptr, var, units=units) 
            else: 
                raw = wrf.getvar(self.fptr, var)
            
            # Save metadata
            self.metadata[ptype]['var'] = var
            self.metadata[ptype]['name'] = raw.description
            self.metadata[ptype]['units'] = raw.units

            # Interpolate, if desired
            if interp_field != None:
                ifield = wrf.getvar(self.fptr, interp_field)
                data = wrf.interplevel(raw, ifield, interp_lvl)
                self.metadata[ptype]['interp'] = '%d-%s ' % (interp_lvl, ifield.units)
            else:
                data = raw
                self.metadata[ptype]['interp'] = ''
  
            # Get lat/lon coordinates
            lat, lon = wrf.latlon_coords(data)
            coords = [wrf.to_np(lat), wrf.to_np(lon)]

            # Extract time if not done so already
            # This should really be done in __init__, but it's a bit tricky finding the time in the 
            # NetCDF4 object
            if not hasattr(self, 'time'):
                self.time = np.datetime_as_string(data.Time.values)[:-10] + ' UTC'

            data = wrf.to_np(data)

        elif (self.outtype == 'upp' or self.outtype == 'stage4'):
            if np.isnan(zind):
                data = self.ds[var]
            else:
                data = self.ds[var][zind, :, :]

            # Save metadata
            self.metadata[ptype]['var'] = var
            self.metadata[ptype]['name'] = data.attrs['long_name']
            self.metadata[ptype]['units'] = data.attrs['units']
            self.metadata[ptype]['interp'] = ''

            # Get lat/lon coordinates
            coords = [self.ds['gridlat_0'].values, self.ds['gridlon_0'].values]
            
            # Create difference fields
            if diff:
                if np.isnan(zind):
                    data2 = self.ds2[var]
                else:
                    data2 = self.ds2[var][zind, :, :]
                data = data - data2
 
        return data, coords, ptype


    def _closest_gpt(self, lon, lat):
        """
        Determine the indices for the model gridpoint closest to the given (lat, lon) coordinate

        Parameters
        ----------
        lon, lat : float
            (lat, lon) coordinate with units of (deg N, deg E)

        Returns
        -------
        i, j : integer
            Indices of the gridpoint closest to (lat, lon)

        """

        if self.outtype == 'wrf':
            lat2d = wrf.getvar(self.fptr, 'lat')
            lon2d = wrf.getvar(self.fptr, 'lon')

        return np.unravel_index(np.argmin((wrflat.values - lat)**2 + (wrflon.values - lon)**2), wrflat.shape)


    def _create_hcrsxn_ax(self, data):
        """
        Add a matplotlib.axes instance to the figure for a horizontal cross sections

        Parameters
        ----------
        data : Array
            Output from _ingest_data() method

        """

        self.ax = self.fig.add_subplot(self.nrows, self.ncols, self.n, projection=self.proj) 


    def config_ax(self, coastlines=True, states=True, grid=True):
        """
        Add cartopy features to plotting axes. Axes must be defined first.

        Parameters
        ----------
        coastlines : boolean, optional
            Option to add coastlines
        states : boolean, optional
            Option to add state borders

        """
         
        if coastlines:
            self.ax.coastlines('50m')

        if states:
            borders = cfeature.NaturalEarthFeature(category='cultural',
                                                   scale='50m',
                                                   facecolor='none',
                                                   name='admin_1_states_provinces')
            self.ax.add_feature(borders, linewidth=0.5, edgecolor='k')

        if grid:
            self.ax.gridlines()


    def contourf(self, var, ingest_kw={}, cntf_kw={}, cbar_kw={}, label_kw={}):
        """
        Plot data using a filled contour plot

        Parameters
        ----------
        var : string
            Variable to plot
        ingest_kw : dict, optional
            Other keyword arguments passed to _ingest_data (key must be a string)
        cntf_kw : dict, optional
            Other keyword arguments passed to contourf (key must be a string)
        cbar_kw : dict, optional
            Other keyword arguments passed to colorbar (key must be a string)
        label_kw : dict, optional
            Other keyword arguments passed to colorbar.set_label (key must be a string)

        """

        data, coords, ptype = self._ingest_data(var, ptype='contourf0', **ingest_kw)

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        self.cax = self.ax.contourf(coords[1], coords[0], data, transform=self.proj, **cntf_kw)

        self.cbar = plt.colorbar(self.cax, ax=self.ax, **cbar_kw)
        self.cbar.set_label('%s%s (%s)' % (self.metadata[ptype]['interp'], 
                                           self.metadata[ptype]['name'], 
                                           self.metadata[ptype]['units']), **label_kw)
    

    def plot_diff(self, var, ingest_kw={}, cntf_kw={}, cbar_kw={}, label_kw={}, auto=True):
        """
        Plot data using a filled contour plot

        Parameters
        ----------
        var : string
            Variable to plot
        ingest_kw : dict, optional
            Other keyword arguments passed to _ingest_data (key must be a string)
        cntf_kw : dict, optional
            Other keyword arguments passed to contourf (key must be a string)
        cbar_kw : dict, optional
            Other keyword arguments passed to colorbar (key must be a string)
        label_kw : dict, optional
            Other keyword arguments passed to colorbar.set_label (key must be a string)
        auto : boolean, optional
            Automatically use the 'bwr' colormap and scale the contour levels so they are centered
            on zero and include the max differences

        """

        data, coords, ptype = self._ingest_data(var, ptype='contourf0', diff=True, **ingest_kw)

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        if auto:
            mx = np.amax(np.abs(data))
            lvls = np.linspace(-mx, mx, 20) 
            self.cax = self.ax.contourf(coords[1], coords[0], data, lvls, transform=self.proj, 
                                        cmap='bwr', **cntf_kw)
        else:
            self.cax = self.ax.contourf(coords[1], coords[0], data, transform=self.proj, **cntf_kw)

        # Compute RMSD
        rmsd = np.sqrt(np.mean(data*data))

        self.cbar = plt.colorbar(self.cax, ax=self.ax, **cbar_kw)
        self.cbar.set_label('diff %s%s (%s)\n[RMSD = %.2e]' % (self.metadata[ptype]['interp'], 
                                                               self.metadata[ptype]['name'], 
                                                               self.metadata[ptype]['units'],
                                                               rmsd), **label_kw)


    def contour(self, var, ingest_kw={}, cnt_kw={}):
        """
        Plot data using contours

        Parameters
        ----------
        var : string
            Variable from model output file to plot
        ingest_kw : dict, optional
            Other keyword arguments passed to _ingest_data (key must be a string)
        cnt_kw : dict, optional
            Other keyword arguments passed to contour (key must be a string)

        """

        data, coords, ptype = self._ingest_data(var, ptype='contour0', **ingest_kw)

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        self.cax = self.ax.contour(coords[1], coords[0], data, transform=self.proj, **cnt_kw)

    
    def barbs(self, xvar, yvar, thin=1, ingest_kw={}, barb_kw={}):
        """
        Plot data using wind barbs

        Parameters
        ----------
        xvar, yvar : string
            Variables from model output file to plot
        thin : integer, optional
            Option to plot every nth barb
        ingest_kw : dict, optional
            Other keyword arguments passed to _ingest_data (key must be a string)
        barb_kw : dict, optional
            Other keyword arguments passed to barb (key must be a string)

        """

        xdata, coords, ptype = self._ingest_data(xvar, ptype='barb0', **ingest_kw)
        ydata, coords, ptype = self._ingest_data(yvar, ptype='barb0', **ingest_kw)

        self.metadata.pop('barb1')

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        self.cax = self.ax.barbs(coords[1][::thin, ::thin], coords[0][::thin, ::thin], 
                                 xdata[::thin, ::thin], ydata[::thin, ::thin], 
                                 transform=self.proj, **barb_kw)

    
    def quiver(self, xvar, yvar, thin=1, ingest_kw={}, qv_kw={}):
        """
        Plot data using vectors

        Parameters
        ----------
        xvar, yvar : string
            Variables from model output file to plot
        thin : integer, optional
            Option to plot every nth barb
        ingest_kw : dict, optional
            Other keyword arguments passed to _ingest_data (key must be a string)
        qv_kw : dict, optional
            Other keyword arguments passed to quiver (key must be a string)

        """

        xdata, coords, ptype = self._ingest_data(xvar, ptype='vector0', **ingest_kw)
        ydata, coords, ptype = self._ingest_data(yvar, ptype='vector0', **ingest_kw)

        self.metadata.pop('vector1')

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        self.cax = self.ax.quiver(coords[1][::thin, ::thin], coords[0][::thin, ::thin], 
                                  xdata[::thin, ::thin], ydata[::thin, ::thin], 
                                  transform=self.proj, **qv_kw)


    def plot(self, lon, lat, plt_kw={}):
        """
        Plot (lat, lon) coordinates

        Parameters
        ----------
        lon, lat : float
            Longitude and latitude coordinates to plot
        plt_kw : dict, optional
            Other keyword arguments passed to plot (key must be a string)

        """

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        self.cax = self.ax.plot(lon, lat, transform=self.proj, **plt_kw)


    def skewt(self, lon, lat, hodo=True, barbs=True, thin=5):
        """
        Plot a Skew-T, log-p diagram for the gridpoint closest to (lat, lon)

        Parameters
        ----------
        lon, lat : float
            Longitude and latitude coordinates for Skew-T (Skew-T is plotted for model gridpoint
            closest to this coordinate)
        hodo : boolean, optional
            Option to plot hodograph inset
        barbs : boolean, optional
            Option to plot wind barbs
        thin : integer, optional
            Plot every x wind barb, where x = thin
    
        """

        # Determine indices closest to (lat, lon) coordinate
        i, j = self._closest_gpt(lon, lat)

        # Extract variables
        p = wrf.getvar(self.fptr, 'p', units='mb')[:, i, j] 
        T = wrf.getvar(self.fptr, 'temp', units='degC')[:, i, j] 
        Td = wrf.getvar(self.fptr, 'td', units='degC')[:, i, j] 
        if (barbs or hodo):
            u = wrf.getvar(self.fptr, 'ua', units='m s-1')[:, i, j]
            v = wrf.getvar(self.fptr, 'va', units='m s-1')[:, i, j]
        if hodo:
            z = wrf.getvar(self.fptr, 'height_agl', units='m')[:, i, j]

        # Create figure
        skew = SkewT(self.fig, rotation=45)

        skew.plot(p, T, 'r', linewidth=2.5)        
        skew.plot(p, Td, 'b', linewidth=2.5)        

        skew.plot_dry_adiabats(linewidth=0.75)
        skew.plot_moist_adiabats(linewidth=0.75)
        skew.plot_mixing_lines(linewidth=0.75)

        skew.ax.set_xlim(-40, 60)
        skew.ax.set_ylim(1000, 100)

        if hodo:

            # Create hodograph axes
            hod = inset_axes(skew.ax, '35%', '35%', loc=1) 
            h = Hodograph(hod, component_range=50.)
            h.add_grid(increment=10) 
 
            # Color-code hodograph based on height AGL
            zbds = [0, 1000, 3000, 6000, 9000]
            colors = ['k', 'r', 'b', 'g']
            for zi, zf, c in zip(zbds[:-1], zbds[1:], colors):
                ind = np.where(np.logical_and(z >= zi, z < zf))[0]
                ind = np.append(ind, ind[-1]+1)
                h.plot(u[ind], v[ind], c=c, linewidth=2)
            ind = np.where(z >= zbds[-1])[0]
            h.plot(u[ind], v[ind], c='goldenrod', linewidth=2) 

        if barbs:
            imax = np.where(p < 100)[0][0]
            skew.plot_barbs(p[:imax:thin], u[:imax:thin], v[:imax:thin])

        # Add title
        loc = '(%.3f $^{\circ}$N, %.3f $^{\circ}$E)' % (lat, lon)
        time = np.datetime_as_string(p.Time.values)[:-10] + ' UTC:\n'
        if (hodo or barbs):
            ttl = r'%s$T$ ($^{\circ}$C), $T_{d}$ ($^{\circ}$C), wind (m s$^{-1}$) at %s' % (time, loc)
        else:
            ttl = r'%s$T$ ($^{\circ}$C), $T_{d}$ ($^{\circ}$C) at %s' % (time, loc)
        skew.ax.set_title(ttl)
 

    def set_lim(self, minlat, maxlat, minlon, maxlon):
        """
        Set plotting limits

        Parameters
        ----------
        minlat, maxlat : float
            Latitude limits
        minlon, maxlon : float
            Longitude limits

        """

        self.ax.set_extent([minlon, maxlon, minlat, maxlat], crs=self.proj)


    def ax_title(self, txt='', **kwargs):
        """
        Create a title for the axes

        Parameters
        ----------
        txt : string, optional
            Text to add to the beginning of the title

        """

        s = '%s %s' % (txt, self.time)
        for k in self.metadata.keys():
            if k[:-1] != 'contourf':
                s = s + '\n%s: %s%s (%s)' % (k, self.metadata[k]['interp'], 
                                             self.metadata[k]['name'], self.metadata[k]['units']) 

        self.ax.set_title(s, **kwargs)
        

"""
End plot_model_data.py
"""
