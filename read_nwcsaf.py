'''Code to read variables from the netCDF files from a selection of the NWCSAF
products

'''

from netCDF4 import Dataset
import datetime
import glob
import numpy as np

rdt_unreadable_file_list = ['RDT/20190319/S_NWC_RDT-CW_MSG4_GuineaCoast-VISIR_20190319T223000Z.nc', 'RDT/20190319/S_NWC_RDT-CW_MSG4_Senegal-VISIR_20190319T214500Z.nc', 'RDT/20190319/S_NWC_RDT-CW_MSG4_GuineaCoast-VISIR_20190319T223000Z_060.nc'] # List of RDT files that are corrputed.
rdt_unreadable_file_list = rdt_unreadable_file_list + [f.replace('.nc', s+'.nc') for f in rdt_unreadable_file_list for s in ['_015', '_030', '_045', '_060']]
rdt_diff_domain_file_list = ['RDT/20190301/S_NWC_RDT-CW_MSG4_GuineaCoast-VISIR_'+datetime.datetime.strftime(datetime.datetime(2019,3,1,0,0,0,0) +datetime.timedelta(seconds=t*60*15), '%Y%m%dT%H%M')+'00Z.nc' for t in range(44)]+['RDT/20190301/S_NWC_RDT-CW_MSG4_Kenya-VISIR_'+datetime.datetime.strftime(datetime.datetime(2019,3,1,0,0,0,0) +datetime.timedelta(seconds=t*60*15), '%Y%m%dT%H%M')+'00Z.nc' for t in range(44)] # Domains had slightly different extents initially, which can confuse things.
rdt_diff_domain_file_list = rdt_diff_domain_file_list + [f.replace('.nc', s+'.nc') for f in rdt_diff_domain_file_list for s in ['_015', '_030', '_045', '_060']]
crr_unreadable_file_list = [
'/gws/nopw/j04/swift/earajr/NWCSAF_GEO/EXIM/20190501/S_NWC_CRR_MSG4_GuineaCoast-VISIR_20190501T223000Z_030.nc', # empty file!
'/gws/nopw/j04/swift/earajr/NWCSAF_GEO/EXIM/20190501/S_NWC_CRR_MSG4_Senegal-VISIR_20190501T193000Z_030.nc', # empty file
'/gws/nopw/j04/swift/earajr/NWCSAF_GEO/EXIM/20190506/S_NWC_CRR_MSG4_Senegal-VISIR_20190506T143000Z_030.nc', # Empty file
'/gws/nopw/j04/swift/earajr/NWCSAF_GEO/EXIM/20190517/S_NWC_CRR_MSG4_Senegal-VISIR_20190517T104500Z_030.nc', #  Empty file
                          ] # List of CRR files that are corrupted.


def set_nwcsaf_directory():
    '''Uses hostname to set directory where NWCSAF data are located

    Returns:
        nwcsaf_dir (str): Location of all NWCSAF data

    '''
    import socket
    myhost = socket.gethostname() # Use to set directory where GPM data is located
    if 'jc.rl.ac.uk' in myhost: # On JASMIN
        nwcsaf_dir = '/gws/nopw/j04/swift/earajr/NWCSAF_GEO/'
    elif 'jasmin' in myhost: # On JASMIN
        nwcsaf_dir = '/gws/nopw/j04/swift/earajr/NWCSAF_GEO/'
    elif 'racc-login' in myhost: # On Reading academic computing cluster
        nwcsaf_dir = '/storage/silver/csspclouds/SWIFT/NWCSAF/'
    else:
        nwcsaf_dir = 'Error! No setting for gpm_dir for hostname = ' + myhost
    return nwcsaf_dir


class nwcsaf_var():
    '''Determine file and variable names

    Convert variable names I use to those used in CRR netCDF files, and 
    identifies correct file to read from, which depends on time, variable name
    and domain

    Args:
        varname (str): My variable name including forecast lead time
        time (datetime object): Valid time
        domain (str): Region of interest, included in filename attribute

    Attributes:
        nc_varname (str): Variable name as in NCWSAF netCDF files
        filename (str): Full name of file to read including all directories

    '''
    def __init__(self,varname, time, domain):
        nwcsaf_dir = set_nwcsaf_directory()
        if varname == 'water path based rain rate':
            self.nc_varname  = 'crrph_intensity'
            self.filename = nwcsaf_dir + 'CRR/' + time.strftime('%Y%m%d/') + 'S_NWC_CRR-Ph_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
        elif varname == 'rain rate':
            self.nc_varname  = 'crr_intensity'
            self.filename = nwcsaf_dir + 'CRR/' + time.strftime('%Y%m%d/') + 'S_NWC_CRR_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
        elif 'extrapolated rain rate' in varname:
            self.nc_varname  = 'crr_intensity'
            self.filename = nwcsaf_dir + 'EXIM/' + time.strftime('%Y%m%d/') + 'S_NWC_CRR_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z_015.nc') #  assumes 15 minte and corrects to other times later
        elif 'extrapolated water path based rain rate' in varname:
            self.nc_varname  = 'crrph_intensity'
            self.filename = nwcsaf_dir + 'EXIM/' + time.strftime('%Y%m%d/') + 'S_NWC_CRR-Ph_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z_015.nc')
        elif varname == 'probability of convective initiation in next 30 minutes':
            self.nc_varname  = 'PC'
            self.filename = nwcsaf_dir + 'PC/' + time.strftime('%Y%m%d/') + 'S_NWC_PC_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
        elif varname == 'K index':
            self.nc_varname  = 'ishai_ki'
            self.filename = nwcsaf_dir + 'iSHAI/' + time.strftime('%Y%m%d/') + 'S_NWC_iSHAI_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
        elif varname == 'Lifted index':
            self.nc_varname  = 'ishai_li'
            self.filename = nwcsaf_dir + 'iSHAI/' + time.strftime('%Y%m%d/') + 'S_NWC_iSHAI_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
        elif varname == 'day flag':
            self.nc_varname  = 'crr_conditions'
            self.filename = nwcsaf_dir + 'CRR/' + time.strftime('%Y%m%d/') + 'S_NWC_CRR_MSG4_'+ domain + time.strftime('*%Y%m%dT%H%M00Z.nc')
#       Now for extrapolated products, return the data valid at "time" based on
#       extrapolating by as many minutes as requested. This requires a different
#       timestamp in the filename
        if '15 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 15)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_'))
        elif '30 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 30)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_')).replace('015.nc','030.nc')
        elif '45 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 45)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_')).replace('015.nc','045.nc')
        elif '60 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 60)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_')).replace('015.nc','060.nc')
        elif '75 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 75)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_')).replace('015.nc','075.nc')
        elif '90 minute' in varname:
            forecast_time = time - datetime.timedelta(seconds = 60 * 90)
            self.filename = self.filename.replace(time.strftime('%Y%m%d/'), forecast_time.strftime('%Y%m%d/')).replace(time.strftime('*%Y%m%dT%H%M00Z_'), forecast_time.strftime('*%Y%m%dT%H%M00Z_')).replace('015.nc','090.nc')


def read_rdt_stormlist(timelist=[datetime.datetime(2019,3,1,12,0,0,0)], domain='GuineaCoast', lead_time=0):
    '''Reads output from RDT into a list of storms at each point in time. 

    Args:
        timelist (list, optional): List of times to read in datetime format
            Defaults to 1 March 2019 12:00
        domain (str, optional)   : String identifying the domain of interest
            Defaults to GuineaCoast
        lead_time (int,optional): Lead time for RDT data, must be one of
            0,15,30,45, or 60. Defaults to 0.

    Returns:
        lons (ndarray, shape(nlat, nlon)): Array of longitude values
        lats (ndarray, shape(nlat, nlon)): Array of latitude values
        stormlist (list): List of storms. Each element of list is a dictionary
            of storm properties corresponding to storms at a single point in
            time. 

    '''
    nwcsaf_dir = set_nwcsaf_directory()
    rdt_dir = nwcsaf_dir + 'RDT/'
    stormlist = []
    lon_max, lon_min, lat_max, lat_min = [], [], [], []
    rdt_skip_file_list = rdt_unreadable_file_list + rdt_diff_domain_file_list
    rdt_skip_file_list = [nwcsaf_dir+f for f in rdt_skip_file_list]
    for t in timelist:
        storm_dict = {}
        if lead_time == 0:
            filename = glob.glob(rdt_dir + t.strftime('%Y%m%d/*')+domain+t.strftime('*%Y%m%dT%H%M00Z.nc'))
        else:
            filename = glob.glob(rdt_dir + t.strftime('%Y%m%d/*')+domain+t.strftime('*%Y%m%dT%H%M00Z_'+'{:03}'.format(lead_time)+'.nc')) #  Note that unlike the other variables, for rdt, the file read matches the forecast time, not the valid time. This facillitates matching of cells between files.
        if (len(filename) == 1):
            if (not(filename[0] in rdt_skip_file_list)):
                ncfile = Dataset(filename[0])
                storm_id_birth = ncfile.variables['NumIdBirth'][:]
                ind = np.where(storm_id_birth == 0)
                storm_id_birth[ind] = ncfile.variables['NumIdCell'][ind]
                initiation_time = np.array([(t - datetime.timedelta(seconds=float(s))) for s in ncfile.variables['Duration'][:]])
                storm_dict.update({'InitiationTime' : initiation_time})
                storm_dict.update({'MyStormID' : np.array([10000*float(s.strftime('%Y%m%d%H%M'))+i for s,i in zip(initiation_time, storm_id_birth)])})
                storm_dict.update({'StormPhase' : ncfile.variables['PhaseLife'][:]})
                storm_dict.update({'LonContour' : ncfile.variables['LonContour'][:]})
                storm_dict.update({'LatContour' : ncfile.variables['LatContour'][:]})
                storm_dict.update({'LatG' : ncfile.variables['LatG'][:]})
                storm_dict.update({'LonG' : ncfile.variables['LonG'][:]})
                storm_dict.update({'Area' : ncfile.variables['Surface'][:]})
                storm_dict.update({'Strength' : ncfile.variables['SeverityIntensity'][:]})
                if isinstance(ncfile.variables['SeverityType'][:], np.ma.MaskedArray): # Sometimes this is a masked array, sometimes not. If masked, read data.
                    storm_dict.update({'SeverityType' : ncfile.variables['SeverityType'][:].data})
                else:
                    storm_dict.update({'SeverityType' : ncfile.variables['SeverityType'][:]})
                if isinstance(ncfile.variables['MvtSpeed'][:], np.ma.MaskedArray):# Sometimes this is a masked array, sometimes not. If masked, read data.
                    storm_dict.update({'MvtSpeed' : ncfile.variables['MvtSpeed'][:].data})
                    storm_dict['MvtSpeed'][ncfile.variables['MvtSpeed'][:].mask] = 0.0
                else:
                    storm_dict.update({'MvtSpeed' : ncfile.variables['MvtSpeed'][:]})
                if isinstance(ncfile.variables['MvtSpeed'][:], np.ma.MaskedArray):# Sometimes this is a masked array, sometimes not. If masked, read data.
                    storm_dict.update({'MvtDirection' : ncfile.variables['MvtDirection'][:].data})
                    storm_dict['MvtDirection'][ncfile.variables['MvtDirection'][:].mask] = 0.0
                else:
                    storm_dict.update({'MvtDirection' : ncfile.variables['MvtDirection'][:]})
                storm_dict.update({'MvtQuality' : ncfile.variables['MvtQuality'][:]})
                stormlist += [storm_dict]
                if lead_time != 0: ncfile = Dataset(filename[0].replace('_{:03}'.format(lead_time), ''))
                lons = ncfile.variables['lon'][:]
                lats = ncfile.variables['lat'][:]
                # Check for domain change (hopefully this only occurs for GuineaCoast region. Pain as it makes code much slower...
                lon_max += [lons.max()]
                lon_min += [lons.min()]
                lat_max += [lats.max()]
                lat_min += [lats.min()]
                if ((np.unique(lon_max).size > 1) | (np.unique(lon_min).size > 1) | (np.unique(lat_max).size > 1) | (np.unique(lat_min).size > 1)):
                    print("Domain size changed!")
                    print("filename=", filename)
                    print("lon_max, lon_min, lat_max, lat_min=", lon_max, lon_min, lat_max, lat_min)
                    quit()
            else:
                stormlist += [{'Missing File' : -9999.}]
        else:
            print("rdt_dir + t.strftime('%Y%m%d/*')+domain+t.strftime('*%Y%m%dT%H%M00Z.nc')=", rdt_dir + t.strftime('%Y%m%d/*')+domain+t.strftime('*%Y%m%dT%H%M00Z.nc'))
            stormlist += [{'Missing File' : -9999.}]
            lons = Dataset('/gws/nopw/j04/swift/earajr/NWCSAF_GEO/RDT/20190331/S_NWC_RDT-CW_MSG4_'+domain+'-VISIR_20190331T233000Z.nc').variables['lon'][:]
            lats = Dataset('/gws/nopw/j04/swift/earajr/NWCSAF_GEO/RDT/20190331/S_NWC_RDT-CW_MSG4_'+domain+'-VISIR_20190331T233000Z.nc').variables['lon'][:]
    return lons, lats, stormlist


def read_spatial_data(timelist=[datetime.datetime(2019,3,1,12,0,0,0)], domain='GuineaCoast', varname='rain rate', good_data_time=datetime.datetime(2019,3,1,15,0,0,0)):
    '''Reads all NWCSAF data except for RDT

    Includes rainfall rates, stability indeces, prop of precip
    and prob of convective initiation 

    Args:
        timelist (list, optional): List of times to read in datetime format
        domain (str, optional)   : String identifying the domain of interest
        varname (str, optional)  :
        good_data_time (datetime object, optional):

    Returns:
        lons (ndarray, shape(nlat, nlon)): Array of longitude values
        lats (ndarray, shape(nlat, nlon)): Array of latitude values
        data (ndarray, shape(ntimes, nlat, nlon)): Array of data values

    '''
    data = []
    for t in timelist:
        var = nwcsaf_var(varname, t, domain)
        filename = glob.glob(var.filename)
        if len(filename) == 1:
            if (not(filename[0] in crr_unreadable_file_list)):
                ncfile = Dataset(filename[0])
                if 'lon' in ncfile.variables.keys(): # In some files, lon is missing, so set lon if present, otherwise, read from file where definitely present
                    lon = ncfile.variables['lon'][:]
                    lat = ncfile.variables['lat'][:]
                else:
                   lon, lat, placeholder = read_spatial_data([datetime.datetime(2019,3,1,15,0,0,0)], domain=domain, varname=varname)
                if var.nc_varname in ncfile.variables.keys():
                    res = ncfile.variables[var.nc_varname][:].astype('float32')
                    res = np.ma.filled(res, -9999.)
                    data += [res]
                else:
                    data += ['Empty']
            else:
                data += ['Empty']
        else:
            data += ['Empty']
    if 'Empty' in data: # Corresponds to missing data. Replace this with array of same size as good data, filled with -9999s
        missing = next((v for v in data if v != 'Empty'), -1)
        if type(missing) == int:
            print('All NWCSAF files missing, unable to generate array of missing data')
            lon, lat, missing = read_spatial_data([good_data_time], domain=domain, varname=varname) # If haven't found any non-missing data to copy array size, then read file that deffinitely doesn't contain missing data.
            missing = missing[0,:,:]
        missing2 = np.copy(missing)
        missing2[:] = -9999.    
        data = [missing2 if d == 'Empty' else d for d in data]
    return lon, lat, np.array(data)
