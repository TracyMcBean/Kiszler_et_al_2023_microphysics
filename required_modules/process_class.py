import xarray as xr
import numpy as np

"""
Creation of a class of which each process can be an instance to 
automatically compute statistics and have methods that are used for each process. 
"""

class Process:
    def __init__(self, data_array: xr.DataArray, threshold: float = 1e-18) -> None:
        self.stats_mean = {"total": 0, "PN": 0, "PD": 0}
        self.freq = {"total": 0, "PN": 0, "PD": 0}
        self.count = {"total": 0, "PN": 0, "PD": 0}
        self.da = data_array
        self.threshold = threshold

        self.__get_stats()

        return

    def __get_stats(self) -> None:
        " Get mean values of the process."

        da_sub = self.da.where(self.da >= self.threshold, drop=True)
        self.stats_mean["total"] = da_sub.mean()
        self.stats_mean["PN"] = da_sub.sel(time=da_sub['time.month'].isin([11, 12, 1, 2])).mean()
        self.stats_mean["PD"] = da_sub.sel(time=da_sub['time.month'].isin([5, 6, 7, 8])).mean()
        return

    def get_perc(self, frozen_mass : xr.DataArray, liquid_mass: xr.DataArray,
                 cloud_type: str = 'mpc') -> None:
        " Get percentage of total occurence of the process."

        # filter only cloudy cases because the xarray structure 
        # fills missing values so that the dimensions are the same
        if cloud_type == 'mpc':
            # assuming only mpc is present
            da_sub = self.da.where(self.da >= self.threshold, drop=True)
            da_cloudy = frozen_mass + liquid_mass
            da_cloudy = da_cloudy.where(da_cloudy > 1e-8, drop=True)
        elif cloud_type == 'liq':
            da_sub = self.da.where((self.da >= self.threshold) & (liquid_mass > 1e-8) & 
                                    (frozen_mass <= 1e-8), drop=True)
            da_cloudy = liquid_mass.where((liquid_mass > 1e-8) & (frozen_mass <= 1e-8), drop=True)
        elif cloud_type == 'frozen':
            da_sub = self.da.where((self.da >= self.threshold) & 
                                    (frozen_mass > 1e-8) & (liquid_mass <= 1e-8), drop=True)
            da_cloudy = frozen_mass.where((frozen_mass > 1e-8) & (liquid_mass <= 1e-8), drop=True)

        self.count["total"] = da_sub.count(dim=['time', 'height_2']).values
        self.count["PN"] = da_sub.sel(time=da_sub['time.month'].isin([11, 12, 1, 2])).count(dim=['time', 'height_2']).values
        self.count["PD"] = da_sub.sel(time=da_sub['time.month'].isin([5, 6, 7, 8])).count(dim=['time', 'height_2']).values
        
        self.freq["total"] = np.round(self.count["total"] \
                                      /da_cloudy.count(dim=['time', 'height_2']).values*100, 2)
        self.freq["PN"] = np.round(self.count["PN"] \
                                   /da_cloudy.sel(time=da_cloudy['time.month'].isin([11, 12, 1, 2])).count(dim=['time', 'height_2']).values*100, 2)
        self.freq["PD"] = np.round(self.count["PD"] \
                                   /da_cloudy.sel(time=da_cloudy['time.month'].isin([5, 6, 7, 8])).count(dim=['time', 'height_2']).values*100, 2)

    def get_normheight_freq(self) -> dict:
        " Get normalized frequency of occurrence with height."

        da_sub = self.da.where(self.da >= self.threshold, drop=True)
        pn_pd_arr = da_sub.count(dim='time').values/da_sub.count(dim='time').values.max()
        pn_arr = da_sub.sel(time=da_sub['time.month'].isin([11, 12, 1, 2])).count(dim='time').values/ \
            da_sub.sel(time=da_sub['time.month'].isin([11, 12, 1, 2])).count(dim='time').values.max()
        pd_arr = da_sub.sel(time=da_sub['time.month'].isin([5, 6, 7, 8])).count(dim='time').values/ \
            da_sub.sel(time=da_sub['time.month'].isin([5, 6, 7, 8])).count(dim='time').values.max()
        
        return {"total": pn_pd_arr, "PN": pn_arr, "PD": pd_arr}

    def get_normtemp_freq(self, temp: xr.DataArray, \
                          T_bins: np.ndarray = np.arange(-40, 5, 1) ) -> dict:
        " Get normalized frequency of occurrence with temperature."
        proc_binned = np.zeros(len(T_bins)-1)
        proc_binned_pn = np.zeros(len(T_bins)-1)
        proc_binned_pd = np.zeros(len(T_bins)-1)

        da_sub = self.da.where(self.da >= self.threshold, drop=True)
        temp = temp.where(self.da >= self.threshold, drop=True)

        proc_values = da_sub.values.flatten()
        T_values = temp.values.flatten()

        proc_values_pn = da_sub.sel(time=da_sub['time.month'].isin([11, 12, 1, 2])).values.flatten()
        T_values_pn = temp.sel(time=temp['time.month'].isin([11, 12, 1, 2])).values.flatten()

        proc_values_pd = da_sub.sel(time=da_sub['time.month'].isin([5, 6, 7, 8])).values.flatten()
        T_values_pd = temp.sel(time=temp['time.month'].isin([5, 6, 7, 8])).values.flatten()

        for i in range(len(T_bins)-1):
            mask = (T_values >= T_bins[i]) & (T_values < T_bins[i+1])
            proc_binned[i] = proc_values[mask].shape[0]

            mask_pn = (T_values_pn >= T_bins[i]) & (T_values_pn < T_bins[i+1])
            proc_binned_pn[i] = proc_values_pn[mask_pn].shape[0]
            
            mask_pd = (T_values_pd >= T_bins[i]) & (T_values_pd < T_bins[i+1])
            proc_binned_pd[i] = proc_values_pd[mask_pd].shape[0]
            
        proc_norm_binned = proc_binned / proc_binned.sum() * 100
        proc_norm_binned_pn = proc_binned_pn / proc_binned_pn.sum() * 100
        proc_norm_binned_pd = proc_binned_pd / proc_binned_pd.sum() * 100

        return {"total": proc_norm_binned, "PN": proc_norm_binned_pn, "PD": proc_norm_binned_pd}

    
