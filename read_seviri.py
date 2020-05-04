'''Code to read SEVIRI data in HRIT format.

'''

import satpy # Requires virtual environmen for reading native (.nat)  and hrit files.
import numpy as np
import datetime
import glob
from pyresample import geometry, bilinear
import os

# The following are necessary to read compressed HRIT files
os.environ['XRIT_DECOMPRESS_PATH']='/home/users/phill/PublicDecompWT/2.06/xRITDecompress/xRITDecompress' # Necessary to read compressed files.
my_tmpdir = '/home/users/phill/tmp/' # should match the unix enviornment variable TMPDIR, setting this to homespace facillitates deleting decompressed files after they have been read, which can otherwise fill up tmp directory.

sev_data_dir1='/gws/nopw/j04/swift/earajr/HRIT_archive/'
sev_data_dir2='/gws/nopw/j04/swift/SEVIRI/' # Second directory to check. Necessary as first directory is incomplete

unavailable_times = (
                     [datetime.datetime(2014,3,2,12,00,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2014,3,3,12,00,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2014,3,4,12,00,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2014,3,26,2,45,0)]
                    +[datetime.datetime(2014,12,8,6,45,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2015,3,1,12,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2015,3,2,12,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2015,3,3,9,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(4)]
                    +[datetime.datetime(2015,3,3,12,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(3)]
                    +[datetime.datetime(2015,3,4,12,15,0)]
                    +[datetime.datetime(2015,3,16,11,30,0)]
                    +[datetime.datetime(2015,7,1,0,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(4)]
                    +[datetime.datetime(2015,10,10,11,45,0)]
                    +[datetime.datetime(2015,10,11,11,45,0)]
                    +[datetime.datetime(2015,10,12,11,45,0)]
                    +[datetime.datetime(2015,10,21,10,0,0)]
                    +[datetime.datetime(2015,11,15,3,30,0)+datetime.timedelta(seconds=60*15*n) for n in range(19)] 
                    +[datetime.datetime(2015,11,16,8,30,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2015,11,25,12,15,0)]
                    +[datetime.datetime(2016,2,29,12,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2016,6,8,14,0,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2016,7,29,2,45,0)]
                    +[datetime.datetime(2017,2,28,12,15,0)]
                    +[datetime.datetime(2017,3,17,21,45,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2017,4,22,22,15,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2017,6,22,11,45,0)]
                    +[datetime.datetime(2017,11,7,7,0,0) + datetime.timedelta(seconds=60*15*n) for n in range(5)]
                    +[datetime.datetime(2018,3,6,12,15,0,0)]
                    +[datetime.datetime(2018,3,21,12,15,0,0)]
                    +[datetime.datetime(2018,6,20,9,15,0) + datetime.timedelta(seconds=60*15*n) for n in range(3)]
                    +[datetime.datetime(2018,7,3,13,30,0)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2018,7,4,4,0,0)]
                    +[datetime.datetime(2018,7,4,4,45,0)]
                    +[datetime.datetime(2018,7,10,23,30,0)+datetime.timedelta(seconds=60*15*n) for n in range(16)]
                    +[datetime.datetime(2018, 9, 24, 12, 30)]
                    +[datetime.datetime(2018, 9, 27, 7, 30)+datetime.timedelta(seconds=60*15*n) for n in range(2)]
                    +[datetime.datetime(2018, 10, 10, 11, 45)]
                    +[datetime.datetime(2018, 10, 12, 11, 45)]
                    +[datetime.datetime(2018, 10, 13, 11, 45)]
                    +[datetime.datetime(2018, 10, 15, 11, 45)]
                    )# List of missing times - these are not available from the EumetSat website, or have multiple lines of missing data.

file_dict = {
             0.6 : 'VIS006',
             0.8 : 'VIS008',
             1.6 : 'IR_016',
             3.9 : 'IR_039',
             6.2 : 'WV_062',
             7.3 : 'WV_073',
             8.7 : 'IR_087',
             9.7 : 'IR_097',
             10.8 : 'IR_108',
             12.0 : 'IR_120',
             13.4 : 'IR_134'
            }


def reflectance_correction(rad_arr, lons, lats):
    '''Make satpy "reflectances" consistent with usual definition

    Reflectance is normally defined as outgoing radiation/incoming radiation.
    In satpy, the denominator is set to a fixed value - the incoming radiation 
    for a solar zenith angle of 0 and earth-sun distance of one. This function
    corrects this using the actual earth-sun distance and the solar zenith 
    angle appropriate for the observation time. See discussion at 
    https://github.com/pytroll/satpy/issues/536 for further details.
  
    Args:
        rad_arr (Xarray): Xarray for one of the two visible channels (0.6, 0.8)                           from satpy scene
        lons (ndarray, shape(nlat,nlon)): Array of longitude values 
        lats (ndarray, shape(nlat,nlon)): Array of longitude values

    Returns:
        rad_arr (Xarray): Input array, with reflectance corrected to account for
                          solar zenith angle and earth-sun distance. 

    '''
    from astropy.units import AU # Requires virtual environment
    from sunpy.sun import sunearth_distance # Requires virtual environment
    from pyorbital.astronomy import sun_zenith_angle
    nx = rad_arr.sizes['x']
    ny = rad_arr.sizes['y']
    mu0 = np.ma.zeros((ny,nx))
    dist_sun_earth = np.ma.zeros(ny)
    Tacq = rad_arr.attrs["end_time"] - rad_arr.attrs["start_time"]
    for j in range(ny) :
        tacq = rad_arr.attrs["start_time"] + datetime.timedelta( seconds=(j/float(ny))*Tacq.seconds )
        mu0[j,:] = np.ma.masked_outside( np.pi*sun_zenith_angle( tacq, lons[j,:], lats[j,:])/180.,0.035, 1, copy=False) # in degrees
        dist_sun_earth[j] = float(sunearth_distance(tacq) / AU)
    rad_arr.values *= ((dist_sun_earth[:,None]**2) / np.cos(mu0)) # sun earth distance in AU.
    return rad_arr


def read_seviri_channel(channel_list, time, subdomain=(), regrid=False, my_area=geometry.AreaDefinition('Tropical Africa', 'Tropical Africa', 'Hatano Equal Area',{'proj' : 'hatano', 'lon_0' : 15.0}, 1732, 1510, (-4330000.,-3775000.,4330000., 3775000.)), interp_coeffs=(), reflectance_correct=False):
    '''Read SEVIRI data for given channels and time

    Includes functionality to subsample or regrid. Requires satpy.
    Assumes SEVIRI files are located in sev_data_dir1 set above, with 
    directory structure sev_data_dir1/Year/YearMonthDay/Hour/

    Args:
        channel_list (list): list of channels to read, see file_dict for 
                             possible values
        time (datetime): SEVIRI file date and time, every 00, 15, 30 or 
                         45 minutes exactly,
        subdomain (tuple, optional): If not empty and regrid is False, then 
                                     tuple values are (West boundary, 
                                     South boundary, East boundary, 
                                     North boundary) Defaults to empty tuple.
        regrid (bool, optional): If True, then data is regriddedonto grid 
                                 defined by my_area. Defaults to False.
        my_area (AreaDefinition, optional): pyresmaple.geometry.AreaDefinition
                                            Only used if regrid=True
                                            Defaults to a Hatano equal area 
                                            projection ~4.5 km resolution
                                            extending from ~33W to ~63E and
                                            ~29S to ~29N.
        interp_coeffs (tuple, optional): Interpolation coefficients that may be
                                         used for bilinear interpolation onto 
                                         new grid. Faccilitates use of same 
                                         coeffcients when regridding operation 
                                         is repeated in multiple calls to 
                                         read_seviri_channel. 
                                         Defaults to empty tuple.
        reflectance_correct(bool, optional): Correct visible reflectances for 
                                             variation in solar zenith angle and
                                             earth-sun distances.
                                             Defaults to False.

    Returns:
        data (dict): Dictionary containing following entries:
                     lons (ndarray, shape(nlat,nlon)): Array of longitude values
                     lats (ndarray, shape(nlat,nlon)): Array of latitude values
                     interp_coeffs (tuple): If regrid is True, then the 
                                            interpolation coefficients are 
                                            returned in this variable to 
                                            speed up future regridding
                     channel (ndarray, shape(nlat,nlon)): Dictionary contains 
                                                          separate entry for 
                                                          each channel in 
                                                          channel_list

    '''
    filenames = []
    sat_names = ['MSG4', 'MSG3', 'MSG2', 'MSG1']
    sat_ind = -1
    if time in unavailable_times:
        raise UnavailableFileError("SEVIRI observations for "+time.strftime("%Y/%m/%d_%H%M")+" are not available")
    while ((len(filenames) == 0) & (sat_ind < len(sat_names)-1)): # Sometimes have data from multiple instruments (e.g. 20160504_1045 has MSG3 and MSG1), this ensures most recent is prioritised.
        sat_ind += 1
        filenames=glob.glob(sev_data_dir1+time.strftime("%Y/%Y%m%d/%H/*")+sat_names[sat_ind]+time.strftime("*EPI*%Y%m%d%H%M-*"))+ glob.glob(sev_data_dir1+time.strftime("%Y/%Y%m%d/%H/*")+sat_names[sat_ind]+time.strftime("*PRO*%Y%m%d%H%M-*"))# PRO and EPI files necessary in all scenarios
        sev_dir = sev_data_dir1+time.strftime("%Y/%Y%m%d/%H/*")+sat_names[sat_ind]
    if len(filenames) == 0: # Try alternative directory for SEVIRI data.
        filenames=glob.glob(sev_data_dir2+time.strftime("%Y/%Y%m%d/%H/*")+time.strftime("*EPI*%Y%m%d%H%M-*"))+ glob.glob(sev_data_dir2+time.strftime("%Y/%Y%m%d/%H/*")+time.strftime("*PRO*%Y%m%d%H%M-*"))# PRO and EPI files necessary in all scenarios
        sev_dir = sev_data_dir2+time.strftime("%Y/%Y%m%d/%H/*")
    if len(filenames) == 0:
        print 'sev_data_dir2+time.strftime("%Y/%Y%m%d/%H/*")+sat_names[sat_ind]+time.strftime("*EPI*%Y%m%d%H%M-*")=', sev_data_dir2+time.strftime("%Y/%Y%m%d/%H/*")+sat_names[sat_ind]+time.strftime("*EPI*%Y%m%d%H%M-*")
        raise MissingFileError("SEVIRI observations for "+time.strftime("%Y/%m/%d_%H%M")+" are missing. Please check if they can be downloaded and if not, add to the list of unavailable times.")
    else:
        for channel in channel_list:
            filenames=filenames + glob.glob(sev_dir+'*'+file_dict[channel]+time.strftime("*%Y%m%d%H%M-*")) # add channels required
        scene = satpy.Scene(reader="seviri_l1b_hrit", filenames=filenames)
        data = {}
        scene.load(channel_list)
        if reflectance_correct:
            lons, lats = scene[channel_list[0]].area.get_lonlats()
            if 0.8 in channel_list:
                scene[0.8] = reflectance_correction(scene[0.8], lons, lats) 
            if 0.6 in channel_list:
                scene[0.6] = reflectance_correction(scene[0.6], lons, lats) 
        if regrid != False:
            lons, lats = my_area.get_lonlats()
            if len(interp_coeffs) == 0:
                interp_coeffs = bilinear.get_bil_info(scene[channel_list[0]].area, my_area, radius=50e3, nprocs=1)
                data.update({'interp_coeffs': interp_coeffs})
            for channel in channel_list:
                data.update({str(channel): bilinear.get_sample_from_bil_info(scene[channel].values.ravel(), interp_coeffs[0], interp_coeffs[1], interp_coeffs[2], interp_coeffs[3], output_shape=my_area.shape)})
        else:
            if len(subdomain) > 0:
                scene = scene.crop(ll_bbox=subdomain)
            lons, lats = scene[channel_list[0]].area.get_lonlats()
            lons = lons[:,::-1] # Need to invert y-axis to get longitudes increasing.
            lats = lats[:,::-1]
            for channel in channel_list:
                data.update({str(channel) : scene[channel].values[:,::-1]})
        data.update({'lons' : lons, 'lats' : lats, 'interp_coeffs' : interp_coeffs})
        # Compressed files are decompressed to TMPDIR. Now tidy up
        delete_list = glob.glob(my_tmpdir+time.strftime("*%Y%m%d%H%M-*"))
        for d in delete_list: os.remove(d)
        return data


class FileError(Exception):
    """Base class for other exceptions

    """
    pass


class UnavailableFileError(FileError):
    """Raised when the file is not available from EumetSat

    """
    pass


class MissingFileError(FileError):
    """Raised when the file is missing, but we still need to
    check whether it is definitely not available from EumetSat

    """
    pass
