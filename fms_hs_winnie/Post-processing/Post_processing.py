import datetime as dt
from operator import le

from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pandas as pd
import numpy as np
from cartopy import crs
import cartopy
import matplotlib.ticker as mticker
import os
import shutil

import pandas as pd
import datetime as dt
from netCDF4 import Dataset
import re



# =====================================================================================================================================================
# ============================================================== Obtaining output files Dry core GCM ==================================================
name = 'exp4_NCEPsymm_noSeason_noTop'
path_output = f'/scratch/brown/castanev/DryCore_Wu/output/{name}/'
path_data = f'{path_output}post_processed/'
path_output_final = f'{path_output}post_processed/output/'
path_figures = f'/home/castanev/DryCore_Wu/Post_processing/'
list_data = np.sort(os.listdir(f'{path_data}'))
list_data = np.delete(list_data, -1)
list_begindates = np.unique(np.sort(np.array([re.findall(r'\d+', list_data[i][:13])[0] for i in range(len(list_data))])))

datat = Dataset(f'{path_data}{list_data[0]}')
lons = np.array(datat['lon'])
levels = np.array(datat['pfull'])



for datee in list_begindates:
    cont_d = cont_m = 0
    for i in list_data: 
        if f'{datee}.atmos_daily.nc.' in i:        
            cont_d = cont_d + 1
            data_i_nc = Dataset(f'{path_data}{i}')
            datat_i = np.array(data_i_nc['temp'])  #(time, pfull, lat, lon) (365, 40, 4, 128)
            datau_i = np.array(data_i_nc['ucomp'])
            datav_i = np.array(data_i_nc['vcomp'])
            lats = np.array(data_i_nc['lat'])
            

            if cont_d == 1: 
                datat_d = datat_i; datau_d = datau_i; datav_d = datav_i; lat_d = lats
            else: 
                datat_d = np.concatenate((datat_d, datat_i), axis = 2)
                datau_d = np.concatenate((datau_d, datau_i), axis = 2)
                datav_d = np.concatenate((datav_d, datav_i), axis = 2)
                lat_d = np.concatenate((lat_d, lats))
        
        
        if f'{datee}.atmos_monthly.nc.' in i:
            cont_m = cont_m + 1
            data_i_nc = Dataset(f'{path_data}{i}')
            datat_i = np.array(data_i_nc['temp'])  #(time, pfull, lat, lon) (12, 40, 4, 128)
            datau_i = np.array(data_i_nc['ucomp'])
            datav_i = np.array(data_i_nc['vcomp'])
            lats = np.array(data_i_nc['lat'])

            if cont_m == 1: 
                datat_m = datat_i; datau_m = datau_i; datav_m = datav_i; lat_m = lats
            else: 
                datat_m = np.concatenate((datat_m, datat_i), axis = 2)
                datau_m = np.concatenate((datau_m, datau_i), axis = 2)
                datav_m = np.concatenate((datav_m, datav_i), axis = 2)
                lat_m = np.concatenate((lat_m, lats))
    

    # Saving daily files
    time = np.arange(1, 366, 1)
    nc_name = f'{datee}.atmos_daily.nc'
    ncfile = Dataset(f'{path_output_final}{nc_name}', 'w')

    ncfile.createDimension('lat', len(lat_d))
    ncfile.createDimension('lon', len(lons))
    ncfile.createDimension('time', len(time))
    ncfile.createDimension('lev', len(levels))

    var_lats = ncfile.createVariable('lat', 'f', ('lat'))
    var_lons = ncfile.createVariable('lon', 'f', ('lon'))
    var_time = ncfile.createVariable('time', 'f', ('time'))
    var_lev = ncfile.createVariable('lev', 'f', ('lev'))

    var_lats[:] = lat_d
    var_lons[:] = lons
    var_time[:] = time
    var_lev[:] = levels

    varr = ncfile.createVariable('temp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr[:,:,:,:] = datat_d[:,:,:,:]

    varr2 = ncfile.createVariable('ucomp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr2[:,:,:,:] = datau_d[:,:,:,:]

    varr3 = ncfile.createVariable('vcomp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr3[:,:,:,:] = datav_d[:,:,:,:]
    ncfile.close()

    # Saving monthly files
    time = np.arange(1, 13, 1)
    nc_name = f'{datee}.atmos_monthy.nc'
    ncfile = Dataset(f'{path_output_final}{nc_name}', 'w')

    ncfile.createDimension('lat', len(lat_m))
    ncfile.createDimension('lon', len(lons))
    ncfile.createDimension('time', len(time))
    ncfile.createDimension('lev', len(levels))

    var_lats = ncfile.createVariable('lat', 'f', ('lat'))
    var_lons = ncfile.createVariable('lon', 'f', ('lon'))
    var_time = ncfile.createVariable('time', 'f', ('time'))
    var_lev = ncfile.createVariable('lev', 'f', ('lev'))

    var_lats[:] = lat_m
    var_lons[:] = lons
    var_time[:] = time
    var_lev[:] = levels

    varr = ncfile.createVariable('temp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr[:,:,:,:] = datat_m[:,:,:,:]

    varr2 = ncfile.createVariable('ucomp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr2[:,:,:,:] = datau_m[:,:,:,:]

    varr3 = ncfile.createVariable('vcomp', 'f', ('time', 'lev', 'lat', 'lon'))
    varr3[:,:,:,:] = datav_m[:,:,:,:]
    ncfile.close()
      

#shutil.rmtree(f'{path_output}history/')


list_data = np.sort(os.listdir(f'{path_data}/output/'))

nc_name = f'00000000.atmos_daily.nc'
ncfile = Dataset(f'{path_output_final}{nc_name}')
lat_d = np.array(ncfile['lat'])
lev   = np.array(ncfile['lev'])
pos_300 = np.where(abs(lev-300) == np.min(abs(lev-300)))[0][0]


for num, i in enumerate(list_begindates): 
    if f'{i}.atmos_daily.nc' in list_data:        
        data_i_nc = Dataset(f'{path_data}/output/{i}.atmos_daily.nc')
        datat_i = np.array(data_i_nc['temp'])[:,-1,:,:]  #(time, lev, lat, lon) 

        if num == 0: datat_d = datat_i
        #elif num == 16: break
        else: datat_d = np.concatenate((datat_d, datat_i), axis = 0)


# Saving daily files
time = np.arange(1, datat_d.shape[0]+1, 1)
nc_name = f't.atmos_daily.nc'
ncfile = Dataset(f'{path_output_final}{nc_name}', 'w')

ncfile.createDimension('lat', len(lat_d))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('time', len(time))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_time = ncfile.createVariable('time', 'f', ('time'))

var_lats[:] = lat_d
var_lons[:] = lons
var_time[:] = time

varr = ncfile.createVariable('temp', 'f', ('time', 'lat', 'lon'))
varr[:,:,:] = datat_d[:,:,:]

ncfile.close()




for num, i in enumerate(list_begindates): 
    if f'{i}.atmos_daily.nc' in list_data:        
        data_i_nc = Dataset(f'{path_data}/output/{i}.atmos_daily.nc')
        datat_i = np.array(data_i_nc['vcomp'])[:,pos_300,:,:]  #(time, pfull, lat, lon) 

        if num == 0: datat_d = datat_i
        #elif num == 16: break
        else: datat_d = np.concatenate((datat_d, datat_i), axis = 0)
     

# Saving daily files
time = np.arange(1, datat_d.shape[0]+1, 1)
nc_name = f'v.atmos_daily.nc'
ncfile = Dataset(f'{path_output_final}{nc_name}', 'w')

ncfile.createDimension('lat', len(lat_d))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('time', len(time))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_time = ncfile.createVariable('time', 'f', ('time'))

var_lats[:] = lat_d
var_lons[:] = lons
var_time[:] = time

varr = ncfile.createVariable('v', 'f', ('time', 'lat', 'lon'))
varr[:,:,:] = datat_d[:,:,:]

ncfile.close()

