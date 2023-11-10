import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import matplotlib.dates as mdates
import seaborn as sns
import math


# Adapted from: https://stackoverflow.com/questions/47222585/matplotlib-generic-colormap-from-tab10
def categorical_cmap(nc, nsc, cmap="tab10", continuous=False):
    """Function to create a categorical colormap

    Parameters
    ----------
    nc : integer
        number of categories
    nsc : integer
        number of subcategories
    cmap : string
        name of colormap
    continuous : boolean
        whether to create a continuous colormap
    
    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        colormap
    """
    if nc > plt.get_cmap(cmap).N:
        raise ValueError("Too many categories for colormap.")
    if continuous:
        ccolors = plt.get_cmap(cmap)(np.linspace(0,1,nc))
    else:
        ccolors = plt.get_cmap(cmap)(np.arange(nc, dtype=int))
    cols = np.zeros((nc*nsc, 3))
    for i, c in enumerate(ccolors):
        chsv = colors.rgb_to_hsv(c[:3])
        arhsv = np.tile(chsv,nsc).reshape(nsc,3)
        arhsv[:,1] = np.linspace(chsv[1],0.25,nsc)
        arhsv[:,2] = np.linspace(chsv[2],1,nsc)
        rgb = colors.hsv_to_rgb(arhsv)
        cols[i*nsc:(i+1)*nsc,:] = rgb
    cmap = colors.ListedColormap(cols)
    return cmap

def compute_percentage(tend):
    """Compute percentage of each tendency.
     For this compute the absolute sum
     over all processes but don't include the total tendency.
    
     Parameters
     ----------
     tend : xarray.Dataset
         Dataset containing the tendencies of single hydrometeor
    
     Returns
     -------
     tend : xarray.Dataset
         Dataset containing the tendencies and the percentage of each tendency
    """
    print(tend)
    tend_sub = tend.filter_by_attrs(name=lambda v: v is not None and not "total" in v)
    # Store all names of processes in a list except the time and height variable
    proc_name = list(tend_sub.variables)
    print(proc_name)
    if "time" in proc_name:
        proc_name.remove("time")
    if "height_2" in proc_name:
        proc_name.remove("height_2")
    # names for percentage variables    
    proc_name_perc = [i + "_perc" for i in proc_name]

    tend_sum = 0
    for i in range(len(proc_name)):
        tend_sum = tend_sum + abs(tend_sub[proc_name[i]])
  
    for i in range(len(proc_name)):
        tend[proc_name_perc[i]] = abs(tend_sub[proc_name[i]].where(tend_sum != 0)) / \
                             tend_sum.where(tend_sum != 0)
    
    # Set long name of percentage variables
    for i in range(len(proc_name)):
        tend[proc_name_perc[i]].attrs["name"] = proc_name_perc[i]
        tend[proc_name_perc[i]].attrs["long_name"] = "Percentage of " + proc_name[i] + " process"
    print(proc_name_perc)
    print(proc_name)
        
    return tend, proc_name_perc

def plot_masses(ds, hlim=100, is_save=False, save_path=None):
    """Function to plot the mass of each hydrometeor
    
    Parameters
    ----------
    ds: xarray.Dataset 
        containing hydrometeors
    hlim: integer
        height limit for plotting
    is_save: boolean
        Whether to save the figure
    save_path: string
        path to save the figure
    
    Returns
    -------
    None
    """
    sns.set_style("whitegrid")

    cloud_lim = 1e-8
    max_val = 1e-3
    min_val = 1e-8

    fig, ax = plt.subplots(6,1, figsize=(8,14))
    my_cmap=categorical_cmap(5,3, cmap="viridis_r", continuous=True)

    # plot wrapper output
    cbar = ax[0].pcolormesh(ds.time, ds.height_2[hlim:], ds.QI.where(ds.QI >= cloud_lim)[:,hlim:].transpose(),
                          norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                          shading='auto',
                          cmap=my_cmap)
    ax[0].set_title("Cloud ice mass")
    ax[0].set_ylabel("Height m")
    ax[0].get_xaxis().set_ticklabels([])
    fig.colorbar(cbar, ax=ax[0], label=" kg kg-1")

    # same for qr, qs, qg, qh, qc
    cbar = ax[1].pcolormesh(ds.time, ds.height_2[hlim:], ds.QC.where(ds.QC >= cloud_lim)[:,hlim:].transpose(),
                      norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                      shading='auto',
                      cmap=my_cmap)
    ax[1].set_title("Cloud droplet mass")
    ax[1].set_ylabel("Height m")
    ax[1].get_xaxis().set_ticklabels([])
    fig.colorbar(cbar, ax=ax[1], label=" kg kg-1")

    cbar = ax[2].pcolormesh(ds.time, ds.height_2[hlim:], ds.QR.where(ds.QR >= cloud_lim)[:,hlim:].transpose(),
                        norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                        shading='auto',
                        cmap=my_cmap)
    ax[2].set_title("Rain mass")
    ax[2].set_ylabel("Height m")
    ax[2].get_xaxis().set_ticklabels([])
    fig.colorbar(cbar, ax=ax[2], label=" kg kg-1")

    cbar = ax[3].pcolormesh(ds.time, ds.height_2[hlim:], ds.QS.where(ds.QS >= cloud_lim)[:,hlim:].transpose(),
                        norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                        shading='auto',
                        cmap=my_cmap)
    ax[3].set_title("Snow mass")
    ax[3].set_ylabel("Height m")
    ax[3].get_xaxis().set_ticklabels([])
    fig.colorbar(cbar, ax=ax[3], label=" kg kg-1")

    cbar = ax[4].pcolormesh(ds.time, ds.height_2[hlim:], ds.QG.where(ds.QG >= cloud_lim)[:,hlim:].transpose(),
                        norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                        shading='auto',
                        cmap=my_cmap)
    ax[4].set_title("Graupel mass")
    ax[4].set_ylabel("Height m")
    ax[4].get_xaxis().set_ticklabels([])
    fig.colorbar(cbar, ax=ax[4], label=" kg kg-1")
    
    cbar = ax[5].pcolormesh(ds.time, ds.height_2[hlim:], ds.QH.where(ds.QH >= cloud_lim)[:,hlim:].transpose(),
                        norm=colors.LogNorm(vmin=min_val,vmax=max_val),
                        shading='auto',
                        cmap=my_cmap)
    ax[5].set_title("Hail mass")
    ax[5].set_ylabel("Height m")
    ax[5].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.colorbar(cbar, ax=ax[5], label=" kg kg-1")

    if is_save:
        plt.savefig(save_path, format="png", dpi=150, bbox_inches="tight")
    else:
        plt.show()

    return
    

def plot_tend_perc(tend, proc_name_perc, hyd_name, hlim=100,
                    is_save=False, save_path=None):
    """Function to plot the percentage of each process contributing to the total tendency
    
    Parameters
    ----------
    tend: xarray.Dataset 
        containing the time and height dimension
    proc_name_perc: list of strings 
        containing the name of the processes
    hyd_name: string
        name of the hydrometeor
    hlim: integer
        height limit for plotting
    is_save: boolean 
        to save the figure
    save_path: string
        path to save the figure
        
    Returns
    -------
    None
       
    Example
    -------
    plot_tend_perc(tend_perc, tend_name, hyd_name,
                   is_save=True, save_path="tend_perc.png")
    """


    fig, ax = plt.subplots(math.ceil(len(proc_name_perc)/2),2, figsize=(15,10))
    for i in range(len(proc_name_perc)):
        cbar = ax[i//2,i%2].pcolormesh(tend.time, tend.height_2[hlim:],
                                tend[proc_name_perc[i]].values[:,hlim:].transpose(), 
                                shading='auto', cmap="tab20b", vmin=0, vmax=1)
        ax[i//2,i%2].set_title(proc_name_perc[i])
        if i < len(proc_name_perc)-2:
          ax[i//2,i%2].get_xaxis().set_ticklabels([])
        else:
            ax[i//2,i%2].set_xlabel("Time UTC")
            ax[i//2,i%2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        if i%2 == 0:
            ax[i//2,i%2].set_ylabel("Height m")
        
    
    fig.tight_layout()
    plt.colorbar(cbar, ax=ax.ravel().tolist(), orientation="horizontal",
                    label="Percentage of total tendency of " + hyd_name )
    
    if is_save:
        fig.savefig(save_path, bbox_inches="tight", pad_inches=0.1,
                     format="png", dpi=150)
    else:
        plt.show()

    return

def plot_tend(tend, hyd_name, hlim=100,
              vmin=-1e-8, vmax=1e-8, my_cmap="PiYG",
               is_save=False, save_path=None):
    """ Function to plot the tendency of each process contributing to the total tendency
    Parameters
    ----------
    tend: xarray.Dataset
        Containing the time and height dimension
    hyd_name: string
        names of the hydrometeor
    hlim: integer
        height limit for plotting
    vmin: float
        minimum value for colorbar
    vmax: float
        maximum value for colorbar
    my_cmap: string
        colormap for plotting
    is_save: boolean 
        to save the figure
    save_path: string
        path to save the figure
    
    Returns
    -------
    None
    """

    # Get process names for plotting
    tend_sub = tend.filter_by_attrs(name=lambda v: v is not None and not "total" in v and not "perc" in v)
    # Store all names of processes in a list except the time and h
    proc_name = list(tend_sub.variables)
    if "time" in proc_name:
        proc_name.remove("time")
    if "height_2" in proc_name:
        proc_name.remove("height_2")


    fig, ax = plt.subplots(math.ceil(len(proc_name)/2),2, figsize=(15,10))
    for i in range(len(proc_name)):
        cbar = ax[i//2,i%2].pcolormesh(tend.time, tend.height_2[hlim:],
                                tend[proc_name[i]].values[:,hlim:].transpose(), 
                                shading='auto', cmap=my_cmap,
                                norm=colors.SymLogNorm(linthresh=1e-15, vmin=vmin, vmax=vmax))
        ax[i//2,i%2].set_title(proc_name[i])
        if i < len(proc_name)-2:
          ax[i//2,i%2].get_xaxis().set_ticklabels([])
        else:
            ax[i//2,i%2].set_xlabel(mdates.DateFormatter('%Y:%m:%D'))
            ax[i//2,i%2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
        if i%2 == 0:
            ax[i//2,i%2].set_ylabel("Height m")
    
    fig.tight_layout()
    plt.colorbar(cbar, ax=ax.ravel().tolist(), orientation="horizontal",
                    label="Tendency of " + hyd_name )
    
    if is_save:
        fig.savefig(save_path, bbox_inches="tight", pad_inches=0.1,
                     format="png", dpi=150)
    else:
        plt.show()
    
    return

