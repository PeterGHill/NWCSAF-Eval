''' Read GPM data in default HDF5 format

'''

import numpy as np
import h5py
import datetime
import glob
import sys

gpm_dir = '/badc/gpm/data/GPM-IMERG-v6/' # Change as appropriate


def read_gpm(timelist, lon_min=-23., lon_max=58., lat_min=-15., lat_max=26., varname='precipitationCal'):
    '''Reads GPM IMERG data for specified times and lat-lon limits

    Args:
        timelist (list): times to be read
        lon_min (float, optonal): Longitude of Western boundary of region of interest. Defaults to -23.
        lon_max (float, optional): Longitude of Eastern boundary of region of interest. Defaults to 58.
        lat_min (float, optional): Latitude of Southern boundary of region of interest. Defaults to -15.
        lat_max (float, optional): Latitude of Northern boundary of region of interest. Defaults to26.
        varname (string, optional): Name of IMERG variable to read. Defaults to precipitationCal

    Returns:
        lon (ndarray, shape(nlon)): Array of longitude values
        lat (ndarray, shape(nlat)): Array of latitude values
        rain (ndarray, shape(ntimes, nlon, nlat)): Array of values for varname

    '''
    if ((varname == 'precipitationNoIRCal') or (varname == 'precipitationNoIRUncal')):
        lon, lat, rain = read_gpm_no_ir(timelist, lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, varname=varname)
        return lon, lat, rain      
    rain = []
    for i, time in enumerate(timelist):
        f = get_gpm_filename(time)
        dataset = h5py.File(f, 'r')
        lon = dataset['Grid']['lon'][:]
        lat = dataset['Grid']['lat'][:]
        ind_lon = np.where((lon >= lon_min) & (lon <= lon_max))[0]
        ind_lat = np.where((lat >= lat_min) & (lat <= lat_max))[0]
        if dataset['Grid'][varname].ndim == 3:
            rain += [dataset['Grid'][varname][0,ind_lon[0]:ind_lon[-1]+1, ind_lat[0]:ind_lat[-1]+1]]
        else:
            print(("dataset['Grid'][varname].ndim=", dataset['Grid'][varname].ndim))
            sys.exit()
    rain = np.ma.masked_array(np.array(rain), mask=(np.array(rain) < 0.0))
    return lon[ind_lon], lat[ind_lat], rain


def read_gpm_no_ir(timelist, lon_min=-23., lon_max=58., lat_min=-15., lat_max=26., varname='precipitationNoIRUnCal'):
    '''Reads rainfall data from GPM, excluding points based on IR measurements.

    Reads rainfall data from GPM including only points that are based on either
    microwave measurements directly or morphing of microwave data. Points where
    the weighting of rainfall based on IR measrurements is nonzero are 
    excluded.

    Args:
        timelist (list): times to be read
        lon_min (float, optonal): Longitude of Western boundary of region of interest. Defaults to -23.
        lon_max (float, optional): Longitude of Eastern boundary of region of interest. Defaults to 58.
        lat_min (float, optional): Latitude of Southern boundary of region of interest. Defaults to -15.
        lat_max (float, optional): Latitude of Northern boundary of region of interest. Defaults to26.
        varname (string, optional): Rainfall variable to read. Either precipitationNoIRUnCal (before calibration with gauge data) or  precipitationNoIRCal (calibrated with gauge data). Defaults to precipitationNoIRUnCal, which is completely independent of IR.

    Returns:
        lon (ndarray, shape(nlon)): Array of longitude values
        lat (ndarray, shape(nlat)): Array of latitude values
        rain (ndarray, shape(ntimes, nlon, nlat)): Array of values for varname

   '''
    if (varname == 'precipitationNoIRCal'):
        varname = 'precipitationCal'
    elif (varname == 'precipitationNoIRUncal'):
        varname = 'precipitationUncal'
    rain = []
    for i, time in enumerate(timelist):
        f = get_gpm_filename(time)
        dataset = h5py.File(f)
        lon = dataset['Grid']['lon'][:]
        lat = dataset['Grid']['lat'][:]
        ind_lon = np.where((lon >= lon_min) & (lon <= lon_max))[0]
        ind_lat = np.where((lat >= lat_min) & (lat <= lat_max))[0]
        if dataset['Grid'][varname].ndim == 3:
            rain_temp = dataset['Grid'][varname][0,ind_lon[0]:ind_lon[-1]+1, ind_lat[0]:ind_lat[-1]+1]
            ir_weight = dataset['Grid']['IRkalmanFilterWeight'][0,ind_lon[0]:ind_lon[-1]+1, ind_lat[0]:ind_lat[-1]+1]
            rain_temp[ir_weight > 0.0] = -9999
            rain += [rain_temp]
        else:
            print(("dataset['Grid'][varname].ndim=", dataset['Grid'][varname].ndim)) # This has been an issue for IMERG netCDF files, may not be for HDF.
            sys.exit()
    rain = np.ma.masked_array(np.array(rain), mask=(np.array(rain) < 0.0))
    return lon[ind_lon], lat[ind_lat], rain


def get_gpm_filename(time):
    '''Identify GPM IMERG HDF file corresponding to given time

    '''
    f = glob.glob(gpm_dir + time.strftime('%Y/%j/*%Y%m%d-S%H%M*.HDF5'))
    if len(f) != 1:
        print(("gpm_dir + time.strftime('%Y/%j/*%Y%m%d-S%H%M*.HDF5')=", gpm_dir + time.strftime('%Y/%j/*%Y%m%d-S%H%M*.HDF5')))
        print(("f=", f))
        sys.exit()
    return f[0]
