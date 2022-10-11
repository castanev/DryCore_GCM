from shutil import copyfile
from scipy import interpolate

import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pandas as pd
import numpy as np
from cartopy import crs
import cartopy
import matplotlib.ticker as mticker
import scipy.interpolate as scp
import os

import pandas as pd
import datetime as dt
from netCDF4 import Dataset
import matplotlib.colors as mcolors

def center_white_anom(cmap, num, bounds, limite):
    import matplotlib as mpl
    barcmap = mpl.cm.get_cmap(cmap, num)
    barcmap.set_bad(color='white', alpha=0.5)
    bar_vals = barcmap(np.arange(num))  # extract those values as an array
    pos = np.arange(num)
    centro = pos[(bounds >= -limite) & (bounds <= limite)]
    for i in centro:
        bar_vals[i] = [1, 1, 1, 1]  # change the middle value
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("new" + cmap, bar_vals)
    return newcmap

def maps1(x, y, minn, maxx, matriz,  cmap, path, norm, units):
    fig = plt.figure(figsize=[8, 6])
    ax = fig.add_subplot(1, 1, 1)

    im = ax.contourf(x, y[:], matriz.T[:,:], levels=20, cmap = cmap, vmin=minn, vmax=maxx, norm=norm)
    r = ax.contour(x, y, matriz.T, levels=20, colors='k', linewidths=0.5)
    ax.invert_yaxis()

    #cbaxes = fig.add_axes([0.2, -0.1, 0.6, 0.030])
    #cb = plt.colorbar(im, orientation="horizontal", pad=0.1, cax=cbaxes, format='%.1f')
    cb = plt.colorbar(im, orientation="horizontal", pad=0.1, format='%.1f', shrink=0.8)
    cb.set_label(units, fontsize=9, color='dimgrey')
    cb.outline.set_edgecolor(None)
    cb.ax.tick_params(labelcolor='dimgrey', color='dimgrey', labelsize=9)
    plt.savefig(path, dpi=200)
    #plt.show()
    plt.close()


source = 'NCEP'
path_home = '/home/castanev/'
path_figures = f'{path_home}DryCore_Wu/data/Figures/'



# ======================================================== Temperature =======================================================================
'''
# Downloading the data of surface temperature
years = pd.date_range(start=dt.date(1948, 1,1), end=dt.date(2020, 1, 1), freq='A')
for year in years.strftime('%Y'):
    print(f'ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/air.{year}.nc')
    filename = wget.download(f'ftp://ftp2.psl.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/air.{year}.nc', out = f'/scratch/brown/castanev/NCEP_data_T')


# Reading the data of t at 1000 hPa, globally, lats and lons in ascending order.
path_data = f'/scratch/brown/castanev/NCEP_data_T/'
list_data = np.sort(os.listdir(f'{path_data}'))

nc_name = f'air.1948.nc'
ncfile = Dataset(f'{path_data}{nc_name}')
lats = np.array(ncfile['lat'])
lons = np.array(ncfile['lon'])
lev = np.array(ncfile['level'])

contt = 0
t = np.zeros([len(list_data),len(lev),len(lats),len(lons)]) * np.NaN
for i in list_data:
    if 'air.' in i:
        contt = contt + 1
        datat_i_nc = Dataset(f'{path_data}{i}')  #(t, lev, lat, lon)
        datat_i = np.array(datat_i_nc['air'])[:, :, :, :]
        timei_t = np.array(datat_i_nc['time']) #Daily
        dates_i = np.array([dt.datetime(1800, 1, 1) + dt.timedelta(hours=int(timei_t[i])) for i in range(len(timei_t))])
        Month_i = np.array([ii.month for ii in dates_i])
        pos_summer = np.where([ii in [6, 7, 8] for ii in Month_i])[0]
        dates_i_summer = dates_i[pos_summer]
        datat_i_summer = datat_i[pos_summer]
        datat_i_summer = np.nanmean(datat_i_summer, axis = 0)
        t[contt - 1, :, :, :] = datat_i_summer
        print(contt)


t_zonal_mean_summer_cl = np.nanmean(t, axis = 0)

#t_zonal_mean_summer_cl_k = pd.DataFrame(data=np.reshape(t_zonal_mean_summer_cl, [t_zonal_mean_summer_cl.shape[0] * t_zonal_mean_summer_cl.shape[1]]))
#t_zonal_mean_summer_cl_k = t_zonal_mean_summer_cl_k - 273.15 # °C
#t_zonal_mean_summer_cl = t_zonal_mean_summer_cl_k.values.reshape(t_zonal_mean_summer_cl.shape[0], t_zonal_mean_summer_cl.shape[1])


li, ls, intervalos, limite, color = t_zonal_mean_summer_cl[0].min(), t_zonal_mean_summer_cl[0].max(), 8, 0, 'coolwarm'
bounds = np.round(np.linspace(li, ls, intervalos), 3)
colormap = center_white_anom(color, intervalos, bounds, limite)

norm = mcolors.DivergingNorm(vmin=t_zonal_mean_summer_cl[0].min(), vmax = t_zonal_mean_summer_cl[0].max(), vcenter=273)
path = path_figures + f'Mean_state_T1000_{source}_3D.png'
maps1(lons, lats, t_zonal_mean_summer_cl[0].min(), t_zonal_mean_summer_cl[0].max() + np.abs(bounds[0] - bounds[1]), t_zonal_mean_summer_cl[0].T, norm, 'coolwarm', path)


# Saving the mean state (T)
# Lats in ascending order
nc_name = f'Mean_state_T_{source}_3D.nc'
ncfile = Dataset(f'{path_home}Data Analysis/Inputs/{nc_name}', 'w')

ncfile.createDimension('lat', len(lats))
ncfile.createDimension('lon', len(lons))
ncfile.createDimension('lev', len(lev))

var_lats = ncfile.createVariable('lat', 'f', ('lat'))
var_lons = ncfile.createVariable('lon', 'f', ('lon'))
var_lev = ncfile.createVariable('lev', 'f', ('lev'))

var_lats[:] = lats[::-1]
var_lons[:] = lons
var_lev[:] = lev

vwnd = ncfile.createVariable('T', 'f', ('lev', 'lat', 'lon'))
vwnd[:, :, :] = t_zonal_mean_summer_cl[:,::-1, :]
ncfile.close()
'''


nc_name = f'Mean_state_T_NCEP_3D.nc'
ncfile = Dataset(f'{path_home}Data Analysis/Inputs/{nc_name}') 
lats_NCEP = np.array(ncfile['lat'])  #ascending (73)
lons_NCEP = np.array(ncfile['lon'])  #(144)
lev_NCEP  = np.array(ncfile['lev'])  #(17) [hPa] surface to top
t_NCEP  = np.array(ncfile['T'])


# ======================================================= teq_ite31_0.9.nc =====================================================
path_data = '/home/castanev/DryCore_Wu/data'
nc_name = f'teq_ite31_0.9.nc'   
ncfile = Dataset(f'{path_data}/{nc_name}')
lev  = np.array(ncfile['pfull'])  #34 [hPa] top to surface
lats = np.array(ncfile['lat'])   #64
lons = np.array(ncfile['lon'])   #128
t = np.array(ncfile['temp']) 


# SYMMETRICAL TEMPERATURE INPUT FILE
# Zonal mean
t_NCEP_zonal_m = np.nanmean(t_NCEP, axis = 2)  # (lev, lats)

# Reversing the levels (following teq_ite31_0.9.nc structure)
t_NCEP_zonal_m_top_sur = t_NCEP_zonal_m[::-1,:]
lev_NCEP_new = lev_NCEP[::-1]  #top to surface

f_interpolate = interpolate.interp2d(lats_NCEP,lev_NCEP_new,t_NCEP_zonal_m_top_sur,kind='cubic')
t_NCEP_zonal_m_top_sur_interp = f_interpolate(lats, lev)

t_NCEP_zonal_m_top_sur_interp_symm = np.repeat(t_NCEP_zonal_m_top_sur_interp[:, :, np.newaxis], len(lons), axis=2)
t_NCEP_zonal_m_top_sur_interp_symm = np.repeat(t_NCEP_zonal_m_top_sur_interp_symm[np.newaxis, :, :, :], 12, axis=0)

ploot =t_NCEP_zonal_m_top_sur_interp_symm[0,:,:,0]
norm = mcolors.DivergingNorm(vmin=ploot.min(), vcenter=273.15, vmax = ploot.max())
maps1(lats, lev,  ploot.min(), ploot.max(), ploot.T, path=f'{path_figures}NCEP_interpolated_zonalm_T.png', cmap='coolwarm', norm=norm, units='T [K]')

nc_name = f'Mean_state_T_{source}_symmetrical_4D.nc'
newfilename = f'{path_home}DryCore_Wu/data/{nc_name}'
copyfile('/home/castanev/DryCore_Wu/data/teq_ite31_0.9.nc',newfilename)

newfile     = Dataset(newfilename,'r+')
newfile.variables['temp'][:,:,:,:]   = t_NCEP_zonal_m_top_sur_interp_symm
newfile.close()



# ASYMMETRICAL TEMPERATURE INPUT FILE
lev_NCEP_new = lev_NCEP[::-1]  #top to surface
t_NCEP_new = t_NCEP[::-1,:,:]

t_NCEP_2d_top_sur_interp = np.zeros(([len(lev_NCEP_new), len(lats), len(lons)]))
for i in range(len(lev_NCEP_new)):
    t_NCEP_2d_top_sur = t_NCEP_new[i,:,:]
    f_interpolate = interpolate.interp2d(lons_NCEP, lats_NCEP,t_NCEP_2d_top_sur,kind='cubic')
    t_NCEP_2d_top_sur_interp[i,:,:] = f_interpolate(lons, lats)

t_NCEP_zonal_m_top_sur_interp_asymm = np.zeros(([len(lev), len(lats), len(lons)]))
for i in range(len(lons)):
    t_NCEP_2d_top_sur = t_NCEP_2d_top_sur_interp[:,:,i]
    f_interpolate = interpolate.interp2d(lats, lev_NCEP_new, t_NCEP_2d_top_sur,kind='cubic')
    t_NCEP_zonal_m_top_sur_interp_asymm[:,:,i] = f_interpolate(lats, lev)

t_NCEP_zonal_m_top_sur_interp_asymm = np.repeat(t_NCEP_zonal_m_top_sur_interp_asymm[np.newaxis, :, :, :], 12, axis=0)

t_NCEP_zonal_m_top_sur_interp_surf = np.mean(t_NCEP_zonal_m_top_sur_interp_asymm, axis=0)[-1, :, :]
norm = mcolors.DivergingNorm(vmin=t_NCEP_zonal_m_top_sur_interp_surf.min(), vcenter=273.15, vmax = t_NCEP_zonal_m_top_sur_interp_surf.max())
maps1(lons, lats,  t_NCEP_zonal_m_top_sur_interp_surf.min(), t_NCEP_zonal_m_top_sur_interp_surf.max(), t_NCEP_zonal_m_top_sur_interp_surf.T, path=f'{path_figures}NCEP_interpolated_asymm_surface_T.png', cmap='coolwarm', norm=norm, units='T [K]')

nc_name = f'Mean_state_T_{source}_asymmetrical_4D.nc'
newfilename = f'{path_home}DryCore_Wu/data/{nc_name}'
copyfile('/home/castanev/DryCore_Wu/data/teq_ite31_0.9.nc',newfilename)

newfile     = Dataset(newfilename,'r+')
newfile.variables['temp'][:,:,:,:]   = t_NCEP_zonal_m_top_sur_interp_asymm
newfile.close()