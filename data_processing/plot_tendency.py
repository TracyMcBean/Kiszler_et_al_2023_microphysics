import xarray as xr
import os.path
import sys
sys.path.append('../..')
import scripts.helper_functions as hf

# This script is used to plot the tendency of hydrometeors
# and the percentage of each tendency.
#
# Example of usage:
# python plot_tendency.py 20190101 201901 wrapper_tend_20190101.nc wrapper_out_20190101.nc

# read in date from command line
YYYYMMDD = sys.argv[1]
YYYYMM = sys.argv[2]
# read in tendfile from command line
tendfile = sys.argv[3]
# read in file from wrapper output with hydrometeor masses
outfile = sys.argv[4]

# Set height limit for plotting 
hlim = 100
print("Haaaaalllloooooooo")
#read in the data from the wrapper tendency output
if not os.path.isfile(tendfile):
    print("File does not exist: " + tendfile)
    exit()

ds_tend = xr.open_dataset(tendfile)
print('THis is ds_tend')
print(ds_tend)

if not os.path.isfile(outfile):
    print("File does not exist: " + outfile)
    exit()
ds_out = xr.open_dataset(outfile)

# only keep periods where there is a cloud
cloud_lim = 1e-8
ds_tend = ds_tend.where(ds_out.QC + ds_out.QS + ds_out.QI + \
                        ds_out.QR + ds_out.QH + ds_out.QG >= cloud_lim)

# Select all variables containing "rain" in their name from ds_tend
qr_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "rain" in v)
qi_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "ice" in v)
# Do same for snow
qs_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "snow" in v)
# Do same for graupel
qg_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "graupel" in v)
# Do same for cloud water
qc_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "cloud" in v)
# Do same for hail
qh_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "hail" in v)
# Do same for specific humidity
qv_tend = ds_tend.filter_by_attrs(name=lambda v: v is not None and "atmo" in v)

# Only keep times where there is a cloud (in model cloud_lim=1e-8)
# cloud_lim = 1e-8
# qi_tend = qi_tend.where(ds_out.QI >= cloud_lim)
# qc_tend = qc_tend.where(ds_out.QC >= cloud_lim)
# qr_tend = qr_tend.where(ds_out.QR >= cloud_lim)
# qs_tend = qs_tend.where(ds_out.QS >= cloud_lim)
# qg_tend = qg_tend.where(ds_out.QG >= cloud_lim)
# qh_tend = qh_tend.where(ds_out.QH >= cloud_lim)
print(qi_tend)

# Get percentage of each tendency 
qi_tend_perc, qi_perc_name = hf.compute_percentage(qi_tend)
qc_tend_perc, qc_perc_name = hf.compute_percentage(qc_tend)
qr_tend_perc, qr_perc_name = hf.compute_percentage(qr_tend)
qs_tend_perc, qs_perc_name = hf.compute_percentage(qs_tend)
qg_tend_perc, qg_perc_name = hf.compute_percentage(qg_tend)
qh_tend_perc, qh_perc_name = hf.compute_percentage(qh_tend)

print(qi_perc_name)
### Plots

# 1. Plot total masses
save_path_mass = "../../../plots/mass_plots/"+YYYYMM+"/"
if not os.path.exists(save_path_mass):
    os.makedirs(save_path_mass)

hf.plot_masses(ds_out, hlim=hlim, is_save=True,
               save_path=save_path_mass+YYYYMMDD+"_masses_.png")


# 2. Plot the tendency in percentage for all hydrometeors
save_path_perc = "../../../plots/percent_plots/"+YYYYMM+"/"
# check if save_path_perc exists, if not create it
if not os.path.exists(save_path_perc):
    os.makedirs(save_path_perc)

hf.plot_tend_perc(qi_tend, qi_perc_name, "ice", hlim=hlim, 
                  is_save=True,
                  save_path=save_path_perc+YYYYMMDD+"_percent_ice.png")
# same for cloud droplets
hf.plot_tend_perc(qc_tend, qc_perc_name, "cloud_droplets", hlim=hlim, 
                  is_save=True,
                  save_path=save_path_perc+YYYYMMDD+"_percent_cloud_droplets.png")
# same for rain
hf.plot_tend_perc(qr_tend, qr_perc_name, "rain", hlim=hlim, 
                  is_save=True,
                  save_path=save_path_perc+YYYYMMDD+"_percent_rain.png")
# same for snow
hf.plot_tend_perc(qs_tend, qs_perc_name, "snow", hlim=hlim,
                    is_save=True,
                    save_path=save_path_perc+YYYYMMDD+"_percent_snow.png")
# same for graupel
hf.plot_tend_perc(qg_tend, qg_perc_name, "graupel", hlim=hlim,
                    is_save=True,
                    save_path=save_path_perc+YYYYMMDD+"_percent_graupel.png")
# same for hail
hf.plot_tend_perc(qh_tend, qh_perc_name, "hail", hlim=hlim,
                    is_save=True,
                    save_path=save_path_perc+YYYYMMDD+"_percent_hail.png")


# 3. Plot the changes of mass
save_path_tend="../../../plots/tendency_plots/"+YYYYMM+"/"
if not os.path.exists(save_path_tend):
    os.makedirs(save_path_tend)

hf.plot_tend(qi_tend, hyd_name="ice", hlim=hlim,
             vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_ice.png")

hf.plot_tend(qc_tend, hyd_name="cloud_droplets", hlim=hlim,
             vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_cloud_droplets.png")

hf.plot_tend(qr_tend, hyd_name="rain", hlim=hlim,
             vmin=-1e-7, vmax=1e-7, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_rain.png")

hf.plot_tend(qs_tend, hyd_name="snow", hlim=hlim,
             vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_snow.png")

hf.plot_tend(qg_tend, hyd_name="graupel", hlim=hlim,
             vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_graupel.png")

hf.plot_tend(qh_tend, hyd_name="hail", hlim=hlim,
             vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
             is_save=True, save_path=save_path_tend+YYYYMMDD+"_tendency_hail.png")