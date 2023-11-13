"""
Create a datasets containing cases with low-level clouds.
Cases with clouds from model output and save them to a netcdf file.
"""

# Import libraries
import xarray as xr
import os
import glob
import numpy as np

# directory where the data is stored
#datadir='../../../data/'
datadir='scratch/low_clouds_datasets/polar_night/'
wrapper_dir = '/scratch/wrapper_output/'
date = '202202'
icon_version='v1'
cloudsel_version='v6'
#output_file = datadir+'/processed/meteo_low_clouds_202211.nc'
output_file = datadir+'LLC_meteo'+date+'_ICON'+icon_version+'_'+cloudsel_version+'.nc'
wr_mass_outfile = datadir+'LLC_wrapper_mass_'+date+'_ICON'+icon_version+'_'+cloudsel_version+'.nc'
wr_tend_outfile = datadir+'LLC_wrapper_tend_'+date+'_ICON'+icon_version+'_'+cloudsel_version+'.nc'
wr_perc_outfile = datadir+'LLC_wrapper_perc_'+date+'_ICON'+icon_version+'_'+cloudsel_version+'.nc'

## 1. Model cloud set
# -------------------

# level index which is used to determine "low clouds"
lev_idx = 111 # ~ 2537m
# minimum number of timesteps in the dataset (10 min)
ts_min = int(10*60/9)

file_string = '/icon_burga_'+icon_version+'/'+date+'*/cn_classification*.nc'
for file in glob.glob(file_string):
# loop over all files in the folder
#for file in os.listdir('../../data/cloudnet_model/'):
    print(file)

    try:
        ds = xr.open_dataset(file)
        # open the corresponding meteo file which is named "METEOGRAM_patch01_YYYYMMDD_awipev.nc"
        date = file.split("_")[2][3:11]
        # replace the cn_classification with METEOGRAM_patch001_${DATE}_awipev.nc
        meteo_file = file.replace("cn_classification"+date+".nc", "METEOGRAM_patch001_"+date+"_awipev.nc")
        ds_meteo = xr.open_dataset(meteo_file)
    except:
        print('Cannot read file: '+file)
        continue

    try:
        # check if wrapper file exists
        wrapper_mass_file = wrapper_dir + '/wrapper_mass/wrapper_mass_'+date+'.nc'
        wrapper_tend_file = wrapper_dir + '/wrapper_tend/wrapper_tend_'+date+'.nc'
        #wrapper_perc_file = wrapper_dir + '/wrapper_tend_perc/percentage_'+date+'.nc'
        ds_wr_mass = xr.open_dataset(wrapper_mass_file)
        ds_wr_tend = xr.open_dataset(wrapper_tend_file)
        #ds_wr_perc = xr.open_dataset(wrapper_perc_file)
    except:
        print('Cannot read wrapper output: '+wrapper_mass_file)
        continue

    print('File opened')

     # check if existing data already contains the new data (same dates) and skip if yes. Only check if year, month and day are the same
    if os.path.isfile(output_file):
        # Open the existing netCDF file with writing permission
        existing_data = xr.open_dataset(output_file, mode="a")
        ds_date = ds.time.dt.strftime('%Y%m%d').values[0]  # Convert date to YYYYMMDD format

        # Check if the new_date exists in the existing_data
        if (ds_date in existing_data.time.dt.strftime('%Y%m%d').values):
            print('Date already in concatenated file')
            continue

    # create a mask where every value above 2 (cloudy) is set to 1 and every value which is 2 (clear) is set to 0
    mask = np.where(ds.nclass.values > 2, 1, 0)
    # higher level - lower level (cloud top is negative -1 and cloud bottom is positive 1)
    mask_change = np.zeros((mask.shape[0], mask.shape[1]))
    mask_change[:,1:] = mask[:,:-1] - mask[:,1:]

    ds["mask_change"] = xr.DataArray(mask_change, dims=['time', 'height_2'],
                            coords={'time': ds.time.values, 'height_2': ds.height_2.values})

    # Criteria: Cloud top above 2.5km and next cloud layer minimum 500m higher.
    to_delete_time_list = []
    for i, time_val in enumerate(ds.time.values):
        # no cloud base above the low cloud found. 
        low_vs_middle_dist = 10000
    
        # if there is no low cloud (no cloud top below 2.5km) delete the timestep
        if np.all(ds.mask_change[i, lev_idx:] != -1):
            to_delete_time_list.append(time_val)
            continue
    
        # Get the next cloud top by going downwards from 2.5km on through the vertical levels
        for lev in range(lev_idx, len(ds.height_2) - 1):
            if ds.mask_change.values[i, lev] == -1:
                # then we found the cloud top height
                # now find the next cloud base height
                lev2 = lev - 1
                while lev2 >= 0:
                    if ds.mask_change.values[i, lev2] == 1:
                        # get distance between the cloud layers
                        low_vs_middle_dist = ds.height_2.values[lev2] - ds.height_2.values[lev]
                        break
                    lev2 -= 1
                break
    
        if low_vs_middle_dist < 500:
            # delete it from the ds_sel_low_clouds_times array
            to_delete_time_list.append(time_val)
        # else set all values above the low cloud to 0 (clear sky) so that they don't affect the analysis.
        #else:
            # replace all values from all variables above the low cloud with 0
         #   ds_meteo = ds_meteo.where(ds_meteo.height_2 > ds.height_2.values[lev])
            #ds_meteo[i, lev+1:] = 0
          #  ds_wr_mass[i, lev+1:] = 0
           # ds_wr_tend[i, lev+1:] = 0
            #ds_wr_perc[i, lev+1:] = 0

    to_delete_time_arr = np.array(to_delete_time_list)

    # Drop all times which aren't low clouds times
    ds_sel_low_clouds_times = ds.drop_sel(time=to_delete_time_arr)

    # Criteria: Filter out small cloud gaps/cloud holes (this didn't work well and wasn't used in the end)
    #if ds_mask.time.size > ts_10min:
    #    ds_mask = ds_mask.where(ds_mask.rolling(time=ts_10min, center=True).sum() > 0, drop=True)
    
    # only select times with hour > 2
    ds_sel_low_clouds_times = ds_sel_low_clouds_times.where(ds_sel_low_clouds_times.time.dt.hour > 2, drop=True)
    
    # From the meteogram file only select these times
    ds_meteo = ds_meteo.sel(time=ds_sel_low_clouds_times.time.values)
    ds_wr_mass = ds_wr_mass.sel(time=ds_sel_low_clouds_times.time.values)
    ds_wr_tend = ds_wr_tend.sel(time=ds_sel_low_clouds_times.time.values)
    #ds_wr_perc = ds_wr_perc.sel(time=ds_sel_low_clouds_times.time.values)

    # if dataset is empty or too small, or doesn't contain the right amount of vertical levels, continue with the next file
    if (ds_sel_low_clouds_times.time.size < ts_min) or \
            (ds_sel_low_clouds_times.height_2.size < 150):
        print('Dataset empty, too small or not enough vertical levels')
        continue

    # store time variable to csv file and append if it exists already
    if os.path.isfile(datadir+'/processed/low_cloud_times_'+cloudsel_version+'.csv'):
        ds_sel_low_clouds_times.time.to_dataframe().to_csv(datadir+'low_cloud_times_'+cloudsel_version+'.csv',
                                                    mode='a', header=False)
    else:
        ds_sel_low_clouds_times.time.to_dataframe().to_csv(datadir+'low_cloud_times_'+cloudsel_version+'.csv')
    
    # store the wrapper tendency files to a netCDF file
    if os.path.isfile(wr_mass_outfile):
        # Open the existing netCDF file with writing permission
        existing_data = xr.open_dataset(wr_mass_outfile, mode="a")
        # Concatenate the new data with the existing data along the time dimension
        concatenated_data = xr.concat([existing_data, ds_wr_mass], dim="time")
        existing_data.close()
        # sort concatenated data by time
        concatenated_data = concatenated_data.sortby('time')
        # Save the concatenated data back to the netCDF file (overwrite)
        concatenated_data.to_netcdf(wr_mass_outfile, unlimited_dims='time', mode='w')
        # Close the datasets
        concatenated_data.close()
    else:
        # Save the new data to a new netCDF file which contains all variables which were in the original meteogram file
        ds_wr_mass.to_netcdf(wr_mass_outfile, unlimited_dims='time', mode='w')

    if os.path.isfile(wr_tend_outfile):
        # Open the existing netCDF file with writing permission
        existing_data = xr.open_dataset(wr_tend_outfile, mode="a")
        # Concatenate the new data with the existing data along the time dimension
        concatenated_data = xr.concat([existing_data, ds_wr_tend], dim="time")
        existing_data.close()
        # sort concatenated data by time
        concatenated_data = concatenated_data.sortby('time')
        # Save the concatenated data back to the netCDF file (overwrite)
        concatenated_data.to_netcdf(wr_tend_outfile, unlimited_dims='time', mode='w')
        # Close the datasets
        concatenated_data.close()
    else:
        # Save the new data to a new netCDF file which contains all variables which were in the original meteogram file
        ds_wr_tend.to_netcdf(wr_tend_outfile, unlimited_dims='time', mode='w')

   
    if os.path.isfile(output_file):
        # Open the existing netCDF file with writing permission
        existing_data = xr.open_dataset(output_file, mode="a")
        # Concatenate the new data with the existing data along the time dimension
        concatenated_data = xr.concat([existing_data, ds_meteo], dim="time")
        existing_data.close()
        # sort concatenated data by time
        concatenated_data = concatenated_data.sortby('time')
        # Save the concatenated data back to the netCDF file (overwrite)
        concatenated_data.to_netcdf(output_file, unlimited_dims='time', mode='w')
        # Close the datasets
        concatenated_data.close()
    else:
        # Save the new data to a new netCDF file which contains all variables which were in the original meteogram file
        ds_meteo.to_netcdf(output_file, unlimited_dims='time', mode='w')

"""
    if os.path.isfile(wr_perc_outfile):
        # Open the existing netCDF file with writing permission
        existing_data = xr.open_dataset(wr_perc_outfile, mode="a")
        # Concatenate the new data with the existing data along the time dimension
        concatenated_data = xr.concat([existing_data, ds_wr_perc], dim="time")
        existing_data.close()
        # sort concatenated data by time
        concatenated_data = concatenated_data.sortby('time')
        # Save the concatenated data back to the netCDF file (overwrite)
        concatenated_data.to_netcdf(wr_perc_outfile, unlimited_dims='time', mode='w')
        # Close the datasets
        concatenated_data.close()
    else:
        # Save the new data to a new netCDF file which contains all variables which were in the original meteogram file
        ds_wr_perc.to_netcdf(wr_perc_outfile, unlimited_dims='time', mode='w')
"""
