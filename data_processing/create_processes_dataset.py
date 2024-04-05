'''
Script to create the dataframe/dataset containing the (mixed-phase, mpc) cloud cases.
Change file names etc to where you want to save the file.
Select whether to do the selection for polar day and/or polar night.

The processes will be grouped into frozen and liquid processes.

Issue: Sometimes the netcdf reading in will cause problems for unknown reasons. 
Rerunning the scripts solves the problem so it must be some bug in the xarray mf_dataset.
'''
import numpy as np
import xarray as xr
import glob
import dask
import os

# for saturation computation
import sys
sys.path.append('../required_modules/')
import icon_mcrph_func as imf

dask.config.set({"array.slicing.split_large_chunks": True}) # this is to avoid memory overload

# set path to save file and whether to do the selection for polar day (pd) and/or polar night (pn)
do_pn = True
do_pd = True

# set whether MPC or all clouds are used
do_mpc = False
output_dir = '/data/processed/nc_files_revisions/'
input_dir = '/low_clouds_datasets/'

if do_pn and do_pd:
    if do_mpc:
        nc_name = 'df_all_tends_mpc_pn_and_pd_v2.nc'
    else:
        nc_name = 'df_all_tends_all_pn_and_pd_v2.nc'
    datapath_pd_wrapper = input_dir+'polar_day/LLC_wrapper_v6_PD/'
    datapath_pd_meteo = input_dir+'polar_day/LLC_meteogram_v6_PD/'
    datapath_pn = input_dir+'polar_night/'
    ds_in = xr.open_mfdataset(glob.glob(datapath_pd_wrapper+"LLC_wrapper_tend_*_ICONv1_v6.nc")+
                       glob.glob(datapath_pd_meteo+"LLC_meteo*_ICONv1_v6.nc")+
                       glob.glob(datapath_pn+"LLC_wrapper_tend_*_ICONv1_v6.nc")+
                       glob.glob(datapath_pn+"LLC_meteo*_ICONv1_v6.nc"),
                   combine='by_coords', engine='netcdf4', parallel=True, chunks={'time': 'auto'})
elif do_pd and not do_pn:
    if do_mpc:
        nc_name = 'df_all_tends_mpc_pd_v1.nc'
    else:
        nc_name = 'df_all_tends_all_pd_v1.nc'
    datapath_pd_wrapper = input_dir+'polar_day/LLC_wrapper_v6_PD/'
    datapath_pd_meteo = input_dir+'polar_day/LLC_meteogram_v6_PD/'
    ds_in = xr.open_mfdataset(glob.glob(datapath_pd_wrapper+"LLC_wrapper_tend_*_ICONv1_v6.nc")+
                       glob.glob(datapath_pd_meteo+"LLC_meteo*_ICONv1_v6.nc"),
                       combine='by_coords', engine='netcdf4', parallel=True, chunks={'time': 'auto'})
elif do_pn and not do_pd:
    if do_mpc:
        nc_name = 'df_all_tends_mpc_pn_v1.nc'
    else:
        nc_name = 'df_all_tends_all_pn_v1.nc'
    datapath_pn = input_dir+'polar_night/'
    ds_in = xr.open_mfdataset(glob.glob(datapath_pn+"LLC_wrapper_tend_*_ICONv1_v6.nc")+
                       glob.glob(datapath_pn+"LLC_meteo*_ICONv1_v6.nc"),
                   combine='by_coords', engine='netcdf4', parallel=True, chunks={'time': 'auto'})
else:
    raise ValueError('You need to select at least one of the options (do_pn or do_pd)')

# set output location
output_file = output_dir+nc_name

# if file exists, then through an error
if os.path.isfile(output_file):
    raise ValueError('File already exists!', output_file)
else:
    print('File does not exist, will proceed to create it')
    pass

# select only levels below 2.5km  
level_idx = 111
ds = ds_in.isel(height_2=slice(level_idx,None), height=slice(level_idx,-1))
# interpolate w from full level to half level
ds['W'] = ds_in.W.interp(height=ds_in.height_2)

# Convert xarray DataArray to numpy arrays because the xarray is massively slow
QI = ds.QI.values
QS = ds.QS.values
QG = ds.QG.values
QH = ds.QH.values
QC = ds.QC.values
QR = ds.QR.values

if do_mpc:
    # Perform the operations using numpy
    mask = ((QI + QS + QG + QH) > 1e-8) & ((QC + QR) > 1e-8)
    # Convert the numpy array back to xarray DataArray
    mask_da = xr.DataArray(mask, dims=ds.QI.dims, coords=ds.QI.coords)

    # Use the mask in the where function
    ds = ds.where(mask_da, drop=True)
else:
    # set all non cloud values to nan
    mask = ((QI + QS + QG + QH + QC + QR) > 1e-8)
    mask_da = xr.DataArray(mask, dims=ds.QI.dims, coords=ds.QI.coords)
    ds = ds.where(mask_da, drop=True)

print("Splitting and summation of processes")
# split sublimation and deposition
# store positive values of depsub_tend in a new variable
ds['deposition_ice_tend']  = ds['depsub_ice_tend'].where(ds['depsub_ice_tend'] > 0)
ds['sublimation_ice_tend'] = ds['depsub_ice_tend'].where(ds['depsub_ice_tend'] < 0)
# do same for snow, graupel and hail
ds['deposition_snow_tend']     = ds['depsub_snow_tend'].where(ds['depsub_snow_tend'] > 0)
ds['sublimation_snow_tend']    = ds['depsub_snow_tend'].where(ds['depsub_snow_tend'] < 0)
ds['deposition_graupel_tend']  = ds['depsub_graupel_tend'].where(ds['depsub_graupel_tend'] > 0)
ds['sublimation_graupel_tend'] = ds['depsub_graupel_tend'].where(ds['depsub_graupel_tend'] < 0)
ds['deposition_hail_tend']     = ds['depsub_hail_tend'].where(ds['depsub_hail_tend'] > 0)
ds['sublimation_hail_tend']    = ds['depsub_hail_tend'].where(ds['depsub_hail_tend'] < 0)
ds['deposition_atmo_tend']     = ds['depsub_atmo_tend'].where(ds['depsub_atmo_tend'] > 0)
ds['sublimation_atmo_tend']    = ds['depsub_atmo_tend'].where(ds['depsub_atmo_tend'] < 0)

# split condensation and evaporation
# store positive values of satad_cloud_tend in a new variable
ds['condensation_cloud_tend'] = ds['satad_cloud_tend'].where(ds['satad_cloud_tend'] > 0)
ds['evaporation_cloud_tend']  = ds['satad_cloud_tend'].where(ds['satad_cloud_tend'] < 0)
ds['condensation_atmo_tend']  = ds['satad_atmo_tend'].where(ds['satad_atmo_tend'] > 0)
ds['evaporation_atmo_tend']   = ds['satad_atmo_tend'].where(ds['satad_atmo_tend'] < 0)

# get all tendencies which mention ice, snow, graupel or hail and not total
frozen_tend = [tend for tend in ds.data_vars if ('ice' in tend or 'snow' in tend or 'graupel' in tend \
                                                        or 'hail' in tend) and 'total' not in tend \
                                                        and 'depsub' not in tend ]
print("Frozen tendencies:", frozen_tend)
# compute sum over all tendencies which mention ice, snow, graupel or hail
#ds['total_mass_tend']= ds[frozen_tend].to_array().sum(dim='variable')
# same for liquid water so with rain and cloud
liquid_tend = [tend for tend in ds.data_vars if ('rain' in tend or 'cloud' in tend) \
                                                       and 'total' not in tend and 'satad' not in tend]
print("Liquid tendencies:",liquid_tend)

atmo_tend = [tend for tend in ds.data_vars if ('atmo' in tend) and 'total' not in tend \
                                 and 'satad' not in tend and 'depsub' not in tend]
print("Vapour tendencies:",atmo_tend)

# frozen processes
deposition_fr_vars = [var for var in frozen_tend if 'deposition' in var ]
sublimation_fr_vars = [var for var in frozen_tend if 'sublimation' in var ]
rime_fr_vars = [var for var in frozen_tend if 'rime' in var ]
homhet_fr_vars = [var for var in frozen_tend if 'homhet' in var ]
c_homfr_fr_vars = [var for var in frozen_tend if 'c_homfr' in var ]
fr_col_fr_vars = [var for var in frozen_tend if 'col_' in var]
fr_eva_fr_vars = [var for var in frozen_tend if 'eva' in var ]
g_to_h_fr_vars = [var for var in frozen_tend if 'g_to_h' in var  ]
r_freeze_fr_vars = [var for var in frozen_tend if 'r_freeze' in var ]
melt_fr_vars = [var for var in frozen_tend if 'melt' in var ]

# liquid processes
au_li_vars= [var for var in liquid_tend if 'auSB' in var ]
ac_li_vars = [var for var in liquid_tend if 'acSB' in var ]
eva_li_vars = [var for var in liquid_tend if 'eva' in var ] # rain evaporation and cloud evaporation both (so includes satad)
c_homfr_li_vars = [var for var in liquid_tend if 'homfr' in var ]
homhet_li_vars = [var for var in liquid_tend if 'homhet' in var ]
rime_li_vars = [var for var in liquid_tend if 'rime' in var ]
r_freeze_li_vars = [var for var in liquid_tend if 'r_freeze' in var ]
melt_li_vars = [var for var in liquid_tend if 'melt' in var ]
cond_li_vars = [var for var in liquid_tend if 'condensation' in var ]
ccn_act_li_vars = [var for var in liquid_tend if 'ccn_act' in var ]

# create a new xarray dataset which contains all of the tendencies for the mixed-phase cloud cases but for frozen and liquid summed up.
ds_tends = xr.Dataset()
ds_tends['deposition_fr'] = ds[deposition_fr_vars].to_array().sum(dim='variable')
ds_tends['sublimation_fr'] = np.abs(ds[sublimation_fr_vars].to_array().sum(dim='variable')) # take absolute value
ds_tends['rime_fr'] = ds[rime_fr_vars].to_array().sum(dim='variable')
ds_tends['homhet_fr'] = ds[homhet_fr_vars].to_array().sum(dim='variable')
ds_tends['c_homfr_fr'] = ds[c_homfr_fr_vars].to_array().sum(dim='variable')
ds_tends['fr_col_fr'] = ds[fr_col_fr_vars].to_array().sum(dim='variable')
ds_tends['fr_eva_fr'] = np.abs(ds[fr_eva_fr_vars].to_array().sum(dim='variable'))
ds_tends['g_to_h_fr'] = ds[g_to_h_fr_vars].to_array().sum(dim='variable')
ds_tends['r_freeze_fr'] = ds[r_freeze_fr_vars].to_array().sum(dim='variable')
ds_tends['melt_fr'] = np.abs(ds[melt_fr_vars].to_array().sum(dim='variable'))

# for liquid the same
ds_tends['au_li'] = ds[au_li_vars].to_array().sum(dim='variable')
ds_tends['ac_li'] = ds[ac_li_vars].to_array().sum(dim='variable')
ds_tends['evaporation_li'] = np.abs(ds[eva_li_vars].to_array().sum(dim='variable'))
ds_tends['c_homfr_li'] = np.abs(ds[c_homfr_li_vars].to_array().sum(dim='variable'))
ds_tends['homhet_li'] = np.abs(ds[homhet_li_vars].to_array().sum(dim='variable'))
ds_tends['rime_li'] = np.abs(ds[rime_li_vars].to_array().sum(dim='variable'))
ds_tends['r_freeze_li'] = np.abs(ds[r_freeze_li_vars].to_array().sum(dim='variable'))
ds_tends['melt_li'] = np.abs(ds[melt_li_vars].to_array().sum(dim='variable'))
ds_tends['condensation_li'] = ds[cond_li_vars].to_array().sum(dim='variable')
ds_tends['ccn_act_li'] = ds[ccn_act_li_vars].to_array().sum(dim='variable')

# now add the masses to the dataset
ds_tends['QI'] = ds.QI
ds_tends['QS'] = ds.QS
ds_tends['QG'] = ds.QG
ds_tends['QH'] = ds.QH
ds_tends['QC'] = ds.QC
ds_tends['QR'] = ds.QR
ds_tends['QV'] = ds.QV
ds_tends['frozen_mass'] = ds.QI + ds.QS + ds.QG + ds.QH
ds_tends['liquid_mass'] = ds.QC + ds.QR

print("Adding Wegener-Bergeron-Findeisen process for MPC cases...")
deposition_fr = ds_tends.deposition_fr.values
evaporation_li = ds_tends.evaporation_li.values
frozen_mass = ds_tends.frozen_mass.values
liquid_mass = ds_tends.liquid_mass.values
ds_tends['wbf'] = xr.DataArray(np.min([deposition_fr, evaporation_li], axis=0),
                             dims=ds_tends.deposition_fr.dims, coords=ds_tends.deposition_fr.coords)
mask = (deposition_fr > 1e-18) & (evaporation_li > 1e-18) & \
                ((frozen_mass > 1e-8) & (liquid_mass > 1e-8))
mask_da = xr.DataArray(mask, dims=ds_tends.frozen_mass.dims, coords=ds_tends.frozen_mass.coords)
ds_tends['wbf'] = ds_tends.wbf.where(mask_da, drop=True)

print("Adding supercooled liquid fraction (SLF)...")
ds_tends['SLF'] = ds_tends.liquid_mass / (ds_tends.liquid_mass + ds_tends.frozen_mass)

print("Adding Saturation...")
### add saturation with respect to ice to the ds
ds_tends['sat_i'] = xr.DataArray(np.zeros(ds['QV'].shape),
                                dims=('time', 'height_2'))
ds_tends['sat_i'].values[:-1,:] = imf.saturation_ice(ds['T'].values[1:,:],
                                                  ds['QV'].values[:-1,:], ds['RHO'].values[:-1,:])
### add saturation with respect to water to the ds
ds_tends['sat_w'] = xr.DataArray(np.zeros(ds['QV'].shape),
                                dims=('time', 'height_2'))
ds_tends['sat_w'].values[:-1,:] = imf.saturation_water(ds['T'].values[1:,:],
                                                    ds['QV'].values[:-1,:], ds['RHO'].values[:-1,:])

# convert T to Celsius
ds_tends['T'] = ds['T'] - 273.15
ds_tends['W'] = ds['W']


# now save the dataset
print("Saving dataset to netCDF file...")
print("File path:", output_file)
ds_tends.to_netcdf(output_file, mode='w', format='NETCDF4', engine='netcdf4')
#df.to_csv(csv_name)


