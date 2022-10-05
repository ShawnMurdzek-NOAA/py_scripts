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

    def __init__(self, fname, model, fig, nrows, ncols, axnum):

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

        if self.model == 'wrf':

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
            coords = [lat, lon]

        # Extract time if not done so already
        # This should really be done in __init__, but it's a bit tricky finding the time in the 
        # NetCDF4 object
        if not hasattr(self, 'time'):
           self.time = np.datetime_as_string(data.Time.values)[:-10] + ' UTC'
 
        return data, coords, ptype


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


    def contourf(self, var, ingest_kw={}, cntf_kw={}, cbar_kw={}, label_kw={}):
        """
        Plot data using a filled contour plot

        Parameters
        ----------
        var : string
            Variable from model output file to plot
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

        if self.model == 'wrf':
            self.cax = self.ax.contourf(wrf.to_np(coords[1]), wrf.to_np(coords[0]), 
                                        wrf.to_np(data), transform=ccrs.PlateCarree(), **cntf_kw)

        self.cbar = plt.colorbar(self.cax, ax=self.ax, **cbar_kw)
        self.cbar.set_label('%s%s (%s)' % (self.metadata[ptype]['interp'], 
                                           self.metadata[ptype]['name'], 
                                           self.metadata[ptype]['units']), **label_kw)


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

        if self.model == 'wrf':
            self.cax = self.ax.contour(wrf.to_np(coords[1]), wrf.to_np(coords[0]), 
                                       wrf.to_np(data), transform=ccrs.PlateCarree(), **cnt_kw)

    
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

        if self.model == 'wrf':
            self.cax = self.ax.barbs(wrf.to_np(coords[1])[::thin, ::thin], 
                                     wrf.to_np(coords[0])[::thin, ::thin], 
                                     wrf.to_np(xdata)[::thin, ::thin], 
                                     wrf.to_np(ydata)[::thin, ::thin], transform=ccrs.PlateCarree(), 
                                     **barb_kw)

    
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

        if self.model == 'wrf':
            self.cax = self.ax.quiver(wrf.to_np(coords[1])[::thin, ::thin], 
                                      wrf.to_np(coords[0])[::thin, ::thin], 
                                      wrf.to_np(xdata)[::thin, ::thin], 
                                      wrf.to_np(ydata)[::thin, ::thin], transform=ccrs.PlateCarree(), 
                                      **qv_kw)

    
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

        self.ax.set_xlim([minlon, maxlon])
        self.ax.set_ylim([minlat, maxlat])


    def ax_title(self, **kwargs):
        """
        Create a title for the axes

        """

        s = self.time
        for k in self.metadata.keys():
            if k[:-1] != 'contourf':
                s = s + '\n%s: %s%s (%s)' % (k, self.metadata[k]['interp'], 
                                             self.metadata[k]['name'], self.metadata[k]['units']) 

        self.ax.set_title(s, **kwargs)
        

"""
End plot_realdata.py
"""
