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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import wrf


#---------------------------------------------------------------------------------------------------
# Define Plotting Class
#---------------------------------------------------------------------------------------------------

class PlotOutput():
    """
    Class that handles plotting real-data model output from either WRF or FV3.

    Parameters
    ----------
    fname : string
        Output file name
    model : string
        Numerical model ('wrf' or 'fv3')
    fig : matplotlib.figure
        Figure that axes are added to
    nrows : integer
        Number of rows of subplots
    ncols : integer
        Number of columns of subplots
    axnum : integer
        Subplot number

    """

    def __init__(self, fname, model, fig, nrows, ncols, axnum, **kwargs):

        self.model = model
        self.fig = fig
        self.nrows = nrows
        self.ncols = ncols
        self.n = axnum

        # Dictionary used to hold metadata
        self.metadata = {}

        if self.model == 'wrf':
            self.fptr = nc.Dataset(fname)
        elif self.model == 'fv3':
            raise ValueError('FV3 model is not supported yet')

            
    def _ingest_data(self, var, units=None, interp_field=None, interp_lvl=None, ptype='none0'):
        """
        Extract a single variable to plot and interpolate if needed.

        Parameters
        ----------
        var : string
            Variable from model output file to plot
        units : string, optional
            Units for var
        interp_field : string, optional
            Interpolate var to a surface with a constant value of interp_field
        interp_lvl : string, optional
            Value of constant-interp_field used during interpolaton
        ptype : string, optional
            Plot type. Used as the key to store metadata

        Returns
        -------
        data : array
            An array containing var
        coords : list
            List of coordinates for data.
                For a horizontal cross section, this list is [lats, lons]

        """

        coords = []

        # Append an integer to ptype if already in use
        n = 1
        while ptype in self.metadata.keys():
            ptype = ptype[:-1] + str(n)
            n = n + 1
        self.metadata[ptype] = {}

        if self.model == 'wrf':

            # Extract raw variable from WRF output file
            if units != None:
                raw = wrf.getvar(self.fptr, var, units=units) 
            else: 
                raw = wrf.getvar(self.fptr, var)

            # Interpolate, if desired
            if interp_field != None:
                ifield = wrf.getvar(self.fptr, var)
                data = wrf.interplevel(raw, ifield, interp_lvl)
                self.metadata[ptype]['interp'] = '%d-%s' % (interp_lvl, ifield.units)
            else:
                data = raw
 
            # Get lat/lon coordinates
            lat, lon = wrf.latlon_coords(data)
            coords = [lat, lon]

        # Extract time if not done so already
        # This should really be done in __init__, but it's a bit tricky finding the time in the 
        # NetCDF4 object
        if not hasattr(self, 'time'):
           self.time = np.datetime_as_string(data.Time.values)[:-10] + ' UTC'
 
        return data, coords


    def _create_hcrsxn_ax(self, data):
        """
        Add a matplotlib.axes instance to the figure for a horizontal cross sections

        Parameters
        ----------
        data : Array
            Output from _ingest_data() method

        """

        proj = wrf.get_cartopy(data)
        self.ax = self.fig.add_subplot(self.nrows, self.ncols, self.n, projection=proj) 


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


    def contourf(self, var, units=None, interp_field=None, interp_lvl=None, cmap='viridis', 
                 lvls=None, extend='neither'):
        """
        Plot data using a filled contour plot

        Parameters
        ----------
        var : string
            Variable from model output file to plot
        units : string, optional
            Units for var
        interp_field : string, optional
            Interpolate var to a surface with a constant value of interp_field
        interp_lvl : string, optional
            Value of constant-interp_field used during interpolaton
        cmap : optional
            Colormap for filled contour plot 
        lvls : array
            Contour levels (setting to None uses default levels)
        extend : string
            Option to extend colorbar. Options are 'max', 'min', 'both', or 'neither'

        """

        data, coords = self._ingest_data(var, units=units, interp_field=interp_field, 
                                         interp_lvl=interp_lvl, ptype='contourf0')

        if not hasattr(self, 'ax'):
            self._create_hcrsxn_ax(data)

        if self.model == 'wrf':
            self.cax = self.ax.contourf(wrf.to_np(coords[1]), wrf.to_np(coords[0]), 
                                        wrf.to_np(data), levels=lvls, cmap=cmap, extend=extend,
                                        transform=ccrs.PlateCarree())

        self.cbar = plt.colorbar(self.cax, ax=self.ax, orientation='horizontal', aspect=30)
        if interp_field != None:
            self.cbar.set_label('%s = %.1f %s (%s)' % (interp_field, interp_lvl, data.description, 
                                data.units))
        else:
            self.cbar.set_label('%s (%s)' % (data.description, data.units))

    
    def set_lim(self, minlat, maxlat, minlon, maxlon):
        """
        Set plotting limits

        Parameters
        ----------
        minlat : float
            Minimum latitude
        maxlat : float
            Maximum latitude
        minlon : float
            Minimum longitude
        maxlon : float
            Maximum longitude

        """

        self.ax.set_xlim([minlon, maxlon])
        self.ax.set_ylim([minlat, maxlat])


    def ax_title(self, size=14):
        """
        Create a title for the axes

        """

        s = self.time 
        for p in self.plots:
            s = '%s %s' % (s, p)

        self.ax.set_title(s, size=size)
        

"""
End plot_realdata.py
"""
