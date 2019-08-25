#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.cm as mpl_cm
import glob
import xesmf as xe
from netCDF4 import Dataset
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import cartopy.feature as cfeature

import cdo as cd

from scripts.functions import *

import cartopy.feature     as cf
import cartopy             as cartopy

import os
import mmap
import re

import warnings
warnings.filterwarnings("ignore")
get_ipython().run_line_magic('matplotlib', 'inline')


# # plot C20 Reanalysis

# In[94]:


# UPDATE FIG 7

# import file
C20 = xr.open_dataset('test2/C20-Reanalysis.cvdp_data.1871-2012_corr2.nc', decode_times=False)

# select variable
C20_nao = C20.variables["nao_pattern_djf"][:]
C20_psl = C20.variables["psl_spatialmean_djf"][:]
C20_tas = C20.variables["tas_spatialmean_djf"][:]
C20_pr = C20.variables["pr_spatialmean_djf"][:]

# select grid
grid_out = xe.util.grid_global(1,1)

# create grid
regridder_M1 = xe.Regridder(C20, grid_out, 'nearest_s2d', reuse_weights=False)

# use new grid
C20_nao_1 = regridder_M1(C20_nao.values)
C20_psl_2 = regridder_M1(C20_psl.values)
C20_tas_3 = regridder_M1(C20_tas.values)
C20_pr_4 = regridder_M1(C20_pr.values)

list_11 = [C20_nao_1, C20_psl_2, C20_tas_3, C20_pr_4]
names_11 = ["C20_nao_1", "C20_psl_2", "C20_tas_3", "C20_pr_4"]


#draw_models([C20_nao_1], [" "], "C20_reanalysis_nao", lons, lats, 1, 1, -5,5)
#draw_models([C20_tas_3], [" "], "C20_reanalysis_tas", lons, lats, 1, 1, -35,35) # coolwarm
draw_models([C20_pr_4], [" "], "C20_reanalysis_pr", lons, lats, 1, 1, 0,8)


# # name lists + lats and lons

# In[8]:


"""lists"""
experiments = ["lig127k", "lgm", "midHolocene", "piControl", "1pctCO2"] #


"""lats and lons needed for plotting"""

lats= np.arange(-90,90,1)
lons= np.arange(-180,180,1)


# # correct inversed models

# In[3]:


"""inverse latitudes in models that are wrong"""

"""def reverse_data(data, file):
        if False: #"AWI-ESM_" in file: #or "MPI-ESM-P" in file:     # add: or "file_name" in file for other models with same problem
        data = data[::-1,:]
        print("PLAPLAPLAPLAPLAPLAPL")
    return data"""


# In[14]:


"""inverse all models for std calculation"""

def reverse_data(data, file):
    if "nc" in file:
        data = data[::-1,:]
    return data


# # FUNCTION: draw models

# In[182]:


""" draw models """
# ACCORDING TO NEED:
# CMAP -->
# nao = bwr
# psl =
# tas = PuOr instead of Spectral
# pr = BrBG
# psl = RdBu
# figsize: 15,9, lig alone 15,4

#letters = ["a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)",
           #"l)","m)","n)","o)","p)","q)","r)","s)","t)","u)"] # FOR MULTIPLOTS

#letters = ["i)","j)","k)","l)"] # FOR SINGLE PLOTS

def draw_models(models, names, name_out, lons, lats, rows, cols, min_val=0, max_val=0):
    
    if(min_val == 0 and max_val == 0):   # if min and max values are not defined when calling on function, find model's min + max
        print("su")
        old_max=-999999
        old_min=999999
        for model in models:
            for a in model:
                new_min=min(a)
                new_max=max(a)
                if old_min>new_min:
                    old_min=new_min
                if old_max<new_max:
                    old_max=new_max
    else:
        old_max=max_val
        old_min=min_val
        
    
    fig = plt.figure(figsize=(15, 10), dpi=100) # one_plot=15,4 or 5, multiplot=15,10 
    levelArray = []
    levelArray = np.linspace(min_val,max_val,29) # set this manually - note necessarily linspace!! do not include 0 in colourbar
    proj  =  ccrs.PlateCarree() # or ccrs.Orthographic(central_longitude=0, central_latitude=90, globe=None) # ccrs.PlateCarree() 
    plt.gcf().subplots_adjust(hspace=0.2, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925) # (hspace=0.2, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925), one plot: bottom=0.15
    print("hello")
    
    i=1
    for model, name, num in zip(models, names, range(len(models))): # letter, letters
        ax = fig.add_subplot(rows, cols, i, projection=proj) # ccrs.Orthographic(0,90) OR ccrs.PlateCarree(central_longitude=0)
        ax.set_extent([-100,30,0,80], proj) #[-100,30,0,80] [-25, 50, 25, 75]
        im = plt.contourf(lons, lats, model[:], transform=proj, vmin=old_min, vmax=old_max, extend="both",
                          cmap=mpl_cm.bwr, levels=levelArray) #pcolormesh
        
        #im = plt.contour(im:colors="k")
        
        # USE ONLY FOR SINGLE PLOTS OTHERWISE COMMENT
        """gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.xlabels_top = False
        gl.xlabels_bottom = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.ylocator = mticker.FixedLocator([20, 60, 80])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 12, 'color': 'black'}"""
        
        #plt.title(letter, loc='left', fontweight='bold', fontsize=15)
        plt.title(name, loc='center', fontsize='large') # size 16 for discussion plots
        plt.gca().coastlines()
        colorbar = fig.colorbar(im, ax=ax, cax=fig.add_axes([0.27,0.13,0.65,0.02]),
                     orientation='horizontal', format='%.0f') #CHANGE 0f om inga decimaler, #cax = [left, bottom, width, height][[0.15,0.0,0.70,0.03]], one plot= [0.27,0.0,0.46,0.03]
        colorbar.ax.xaxis.set_label_position('top')
        colorbar.ax.set_xlabel("hPa", fontsize=15) # CHANGE ACCORDINGLY °C/σ
        i=i+1
        
    
    plt.savefig(name_out+".png", bbox_inches="tight")


# In[54]:


""" draw models """
# ACCORDING TO NEED:
# CMAP -->
# nao = bwr
# tas = PuOr instead of Spectral
# pr = BrBG
# psl = RdBu
# figsize: 15,9, lig alone 15,4


letters = [" "]

def draw_models_one(models, names, name_out, lons, lats, rows, cols, min_val=0, max_val=0):
    
    if(min_val == 0 and max_val == 0):   # if min and max values are not defined when calling on function, find model's min + max
        print("su")
        old_max=-999999
        old_min=999999
        for model in models:
            for a in model:
                new_min=min(a)
                new_max=max(a)
                if old_min>new_min:
                    old_min=new_min
                if old_max<new_max:
                    old_max=new_max
    else:
        old_max=max_val
        old_min=min_val
        
    
    fig = plt.figure(figsize=(15, 5), dpi=100) # one_plot=15,4 or 5
    levelArray = []
    levelArray = np.linspace(min_val,max_val,29) # set this manually - note necessarily linspace!! do not include 0 in colourbar
    proj  =  ccrs.PlateCarree() # or ccrs.Orthographic(central_longitude=0, central_latitude=90, globe=None) # ccrs.PlateCarree() 
    plt.gcf().subplots_adjust(hspace=0.2, wspace=0.05, top=0.95, bottom=0.07, left=0.075, right=0.925) # (hspace=0.2, wspace=0.05, top=0.95, bottom=0.05, left=0.075, right=0.925), one plot: bottom=0.15
    print("hello")
    
    i=1
    for model, name, letter, num in zip(models, names, letters, range(len(models))):
        ax = fig.add_subplot(rows, cols, i, projection=proj) # ccrs.Orthographic(0,90) OR ccrs.PlateCarree(central_longitude=0)
        ax.set_extent([-15,35,30,50], proj) #[-100,30,0,80] [-25,50,25,75] NEU: [-15,35,50,70], MED: [-15,35,30,50]
        im = plt.contourf(lons, lats, model[:], transform=proj, vmin=old_min, vmax=old_max, extend="both",
                          cmap=mpl_cm.bwr, levels=levelArray) #pcolormesh
        
        #im = plt.contour(im:colors="k")
        
        # USE ONLY FOR SINGLE PLOTS OTHERWISE COMMENT
        """gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
        gl.xlabels_top = False
        gl.xlabels_bottom = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.ylocator = mticker.FixedLocator([20, 60, 80])
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.ylabel_style = {'size': 12, 'color': 'black'}"""
        
        plt.title(letter, loc='left', fontweight='bold', fontsize=15)
        plt.title(name, loc='center', fontsize='large')
        plt.gca().coastlines()
        colorbar = fig.colorbar(im, ax=ax, cax=fig.add_axes([0.15,0.0,0.70,0.02]),
                     orientation='horizontal', format='%.1f') #CHANGE 0f om inga decimaler, #cax = [left, bottom, width, height][[0.15,0.0,0.70,0.03]], one plot= [0.27,0.0,0.46,0.03]
        colorbar.ax.xaxis.set_label_position('top')
        colorbar.ax.set_xlabel("hPa", fontsize=15) # CHANGE ACCORDINGLY °C/σ
        i=i+1
        
    
    plt.savefig(name_out+".png", bbox_inches="tight")


# # FUNCTION: model evaluation: identify difference field between C20 and each model in piControl

# In[7]:


""" identify difference field for model evaluation - between C20 and each model in PiControl"""

def identify_difference_field(exp_name, not_regression=True):

    C20 = xr.open_dataset('test2/C20-Reanalysis.cvdp_data.1871-2012_corr2.nc', decode_times=False)
    
    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    var_dir = {}
    # r=root, d=directories, f = files
    nao= {}
    tas= {}
    pr= {}
    nao_tas = {}
    nao_pr = {}
    
    grid_out = xe.util.grid_global(1,1)
    
    regridder = xe.Regridder(C20, grid_out, 'nearest_s2d', reuse_weights=False) # remember to change to False
    
    C20_nao = regridder(C20.variables.get("nao_pattern_djf").values)
    C20_tas = regridder(C20.variables.get("tas_spatialmean_djf").values)
    C20_pr = regridder(C20.variables.get("pr_spatialmean_djf").values)
    C20_nao_tas = regridder(C20.variables.get("nao_tas_regression_djf").values)
    C20_nao_pr = regridder(C20.variables.get("nao_pr_regression_djf").values)
    
    for r, d, f in os.walk(path):
        for file in f:
            if exp_name in file:
                pathtofile=os.path.join(r, file)
               
                open_file = xr.open_dataset(pathtofile, decode_times=False)
                
                regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False) # remember to change to False
                
                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)
                
                pathtofile = pathtofile[1:-5]
                
                pathtofile = pathtofile.replace("_", " ")
                
                #print(pathtofile)
                
                if not_regression:
                    
                    temp = open_file.variables.get("nao_pattern_djf") # THIS ONE NOT NECESSARILY NEEDED (DEPENDS)
                    if temp is not None:
                        data = regridder(open_file.variables.get("nao_pattern_djf").values)-C20_nao
                        data = reverse_data(data, file)
                        nao[pathtofile] = data

                    temp = open_file.variables.get("tas_spatialmean_djf")
                    if temp is not None:
                        data = regridder(open_file.variables.get("tas_spatialmean_djf").values)-C20_tas
                        data = reverse_data(data, file)
                        tas[pathtofile] = data

                    temp = open_file.variables.get("pr_spatialmean_djf")
                    if temp is not None:
                        data = regridder(open_file.variables.get("pr_spatialmean_djf").values)-C20_pr
                        data = reverse_data(data, file)
                        pr[pathtofile] = data
                    
                else:
                    
                    temp = open_file.variables.get("nao_tas_regression_djf")
                    if temp is not None:
                        data = regridder(open_file.variables.get("nao_tas_regression_djf").values)-C20_nao_tas
                        data = reverse_data(data, file)
                        nao_tas[pathtofile] = data
                        
                    temp = open_file.variables.get("nao_pr_regression_djf")
                    if temp is not None:
                        data = regridder(open_file.variables.get("nao_pr_regression_djf").values)-C20_nao_pr
                        data = reverse_data(data, file)
                        nao_pr[pathtofile] = data
     
    if not_regression:
        var_dir["nao_pattern_djf"]=nao
        var_dir["tas_spatialmean_djf"]=tas
        var_dir["pr_spatialmean_djf"]=pr
    
    else:
        var_dir["nao_tas_regression_djf"]=nao_tas
        var_dir["nao_pr_regression_djf"]=nao_pr
    
    print("workwork...")
    
    return var_dir

# add if clause for models with "error"


# In[41]:


### model evaluation plots ###

variables_dir = identify_difference_field("piControl")
draw_models(variables_dir["nao_pattern_djf"].values(), variables_dir["nao_pattern_djf"].keys(), "nao_mod_evaluation", lons, lats, 5, 5, -1, 2)
#draw_models(variables_dir["tas_spatialmean_djf"].values(), variables_dir["tas_spatialmean_djf"].keys(), "tas_mod_evaluation", lons, lats, 5, 5, -6, 7)
#draw_models(variables_dir["pr_spatialmean_djf"].values(), variables_dir["pr_spatialmean_djf"].keys(), "pr_mod_evaluation", lons, lats, 5, 5, -7, 8)


# In[34]:


draw_models(variables_dir["nao_pattern_djf"].values(), variables_dir["nao_pattern_djf"].keys(), "nao_mod_evaluation", lons, lats, 5, 5, -1, 2)


# In[119]:


### regression evaluation plots ###
variables_dir = identify_difference_field("piControl", False)

draw_models(variables_dir["nao_tas_regression_djf"].values(), variables_dir["nao_tas_regression_djf"].keys(), "nao_tas_mod_evaluation", lons, lats, 5, 5,-1,2)
draw_models(variables_dir["nao_pr_regression_djf"].values(), variables_dir["nao_pr_regression_djf"].keys(), "nao_pr_mod_evaluation", lons, lats, 5, 5,-1,1)


# In[97]:


### NORMAL model evaluation with ensemble means ###

means_list2 = []
keys_list1 = []
difference_dir = identify_difference_field("piControl")

#print(difference_dir)

for var in difference_dir:
    #print(difference_dir[var].values())
    time_period_np = list(difference_dir[var].values())
    time_period_mean = np.nanmean(time_period_np, axis=0)
    #print(time_period_mean.shape)
    time_period_np.append(time_period_mean)
    #np.append(time_period_np, [time_period_mean], axis=0)
    keys = list(difference_dir[var].keys())
    #print(keys)
    keys.append("Ensemble Mean") #var+" mean"

    means_list2.append(time_period_np)
    keys_list1.append(keys)


# In[114]:


# draw above cell

#draw_models(means_list2[0], keys_list1[0], "nao_mod_evaluation_ens", lons, lats, 5, 5, min_val=-3.5, max_val=3.5)
#draw_models(means_list2[1], keys_list1[1], "tas_mod_evaluation_ens", lons, lats, 5, 5, min_val=-8.2, max_val=8.2)
draw_models(means_list2[2], keys_list1[2], "pr_mod_evaluation_ens", lons, lats, 5, 5, min_val=-4.7, max_val=4.7)


# In[11]:


### REGRESSION model evaluation with ensemble means ###

means_list2 = []
keys_list1 = []
difference_dir = identify_difference_field("piControl", False)

for var in difference_dir:  
    time_period_np = list(difference_dir[var].values())
    time_period_mean = np.nanmean(time_period_np, axis=0)
    #print(time_period_mean.shape)
    time_period_np.append(time_period_mean)
    #np.append(time_period_np, [time_period_mean], axis=0)
    keys = list(difference_dir[var].keys())
    #print(keys)
    keys.append(var+" mean")

    means_list2.append(time_period_np)
    keys_list1.append(keys)
        

draw_models(means_list2[0], keys_list1[0], "nao_tas_mod_evaluation_ens", lons, lats, 5, 5, min_val=-1, max_val=2)
draw_models(means_list2[1], keys_list1[1], "nao_pr_mod_evaluation_ens", lons, lats, 5, 5, min_val=-0.5, max_val=0.6)


# # FUNCTION: identify all models for all experiments and variables

# In[11]:


"""identify all models for each time period and each variable""" 

def identify_ensemble_members(exp_name, variable="all", not_regression=True, not_timeseries=True):

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    var_dir = {}
    # r = root, d = directories, f = files
    nao= {}
    psl= {}
    pr= {}
    tas= {}
    nao_pr = {}
    nao_tas = {}
    nao_time = {}
    
    grid_out = xe.util.grid_global(1,1)
    
    for r, d, f in os.walk(path):
        for file in f:
            if exp_name in file:
                
                pathtofile = os.path.join(r, file)
               
                open_file = xr.open_dataset(pathtofile, decode_times=False)
                    
                regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)
                
                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)
                
                pathtofile = pathtofile[1:-5]
                
                #pathtofile = pathtofile.replace("_", " ")
                
                print(pathtofile)
                
                if not_regression and not_timeseries:
                    if variable == "nao_pattern_djf" or variable == "all":
                        temp = open_file.variables.get("nao_pattern_djf")
                        if temp is not None:
                            percent = open_file.variables.get("nao_pattern_djf").attrs["pcvar"][:-1]
                            data = regridder(open_file.variables.get("nao_pattern_djf").values)
                            data = reverse_data(data, file)
                            nao[pathtofile+"    "+percent] = data
                            
                    if variable == "psl_spatialmean_djf" or variable == "all":
                        temp = open_file.variables.get("psl_spatialmean_djf")
                        if temp is not None:
                            data = regridder(open_file.variables.get("psl_spatialmean_djf").values)
                            data = reverse_data(data, file)
                            psl[pathtofile]=data 
                            
                    if variable == "tas_spatialmean_djf" or variable == "all":
                        temp = open_file.variables.get("tas_spatialmean_djf")
                        if temp is not None:
                            data = regridder(open_file.variables.get("tas_spatialmean_djf").values)
                            data = reverse_data(data, file)
                            tas[pathtofile]=data
                            
                    if variable == "pr_spatialmean_djf" or variable == "all":
                        temp = open_file.variables.get("pr_spatialmean_djf")
                        if temp is not None:
                            data = regridder(open_file.variables.get("pr_spatialmean_djf").values)
                            data = reverse_data(data, file)
                            pr[pathtofile]= data
                        
                elif not_timeseries:
                    if variable == "nao_pr_regression_djf" or variable == "all":
                        temp = open_file.variables.get("nao_pr_regression_djf")
                        if temp is not None:
                            data = regridder(open_file.variables.get("nao_pr_regression_djf").values)
                            data = reverse_data(data, file)
                            nao_pr[pathtofile]= data
                    if variable == "nao_tas_regression_djf" or variable == "all":   
                        temp = open_file.variables.get("nao_tas_regression_djf")
                        if temp is not None:
                            data = regridder(open_file.variables.get("nao_tas_regression_djf").values)
                            data = reverse_data(data, file)
                            nao_tas[pathtofile]= data
                else:
                    
                    temp = open_file.variables.get("nao_timeseries_djf")
                    if temp is not None:
                        data = open_file.variables.get("nao_timeseries_djf").values
                        nao_time[pathtofile]= data
                
            
    if not_regression and not_timeseries:
        
        var_dir["nao_pattern_djf"]=nao
        var_dir["psl_spatialmean_djf"]=psl
        var_dir["tas_spatialmean_djf"]=tas
        var_dir["pr_spatialmean_djf"]=pr
    
    elif not_timeseries:
        
        var_dir["nao_pr_regression_djf"]=nao_pr
        var_dir["nao_tas_regression_djf"]=nao_tas
        
    else:
        
        var_dir["nao_timeseries_djf"]=nao_time

    
    return var_dir


# In[8]:


### ------ draw all models in each time period and all variables ------- ###

for a in experiments:
    print("\t",a)
    variables_dir = identify_ensemble_members(a)
    draw_models(variables_dir["nao_pattern_djf"].values(), variables_dir["nao_pattern_djf"].keys(), "moimoi", lons, lats, 4, 5, min_val=-5, max_val=5)
    draw_models(variables_dir["psl_spatialmean_djf"].values(), variables_dir["psl_spatialmean_djf"].keys(), "moimoi2", lons, lats, 4, 5, min_val=980, max_val=1040)
    draw_models(variables_dir["tas_spatialmean_djf"].values(), variables_dir["tas_spatialmean_djf"].keys(), "moimoi3", lons, lats, 4, 5, min_val=-40, max_val=30)
    draw_models(variables_dir["pr_spatialmean_djf"].values(), variables_dir["pr_spatialmean_djf"].keys(), "moimoi4", lons, lats, 4, 5, min_val=0, max_val=15)


# In[43]:


### RUN HERE FOR CHANGES IN CELL ABOVE TO FIT NEEDS ###

variables_dir = identify_ensemble_members("piControl")

draw_models(variables_dir["nao_pattern_djf"].values(), variables_dir["nao_pattern_djf"].keys(), "moimoi", lons, lats, 4, 5, min_val=-5, max_val=5)
#draw_models(variables_dir["psl_spatialmean_djf"].values(), variables_dir["psl_spatialmean_djf"].keys(), "moimoi2", lons, lats, 4, 5, min_val=980, max_val=1040)
#draw_models(variables_dir["tas_spatialmean_djf"].values(), variables_dir["tas_spatialmean_djf"].keys(), "moimoi3", lons, lats, 4, 5, min_val=-40, max_val=30)
#draw_models(variables_dir["pr_spatialmean_djf"].values(), variables_dir["pr_spatialmean_djf"].keys(), "moimoi4", lons, lats, 4, 5, min_val=0, max_val=15)


# In[324]:


draw_models(variables_dir["nao_pattern_djf"].values(), variables_dir["nao_pattern_djf"].keys(), "moimoi", lons, lats, 4, 5, min_val=-5, max_val=5)


# In[183]:


### ----- draw all models along with their ensemble means ----- ###

### UPDATE HÄR - skriv denna som funktion eller lägg till "mean" i funktionen identify_ensemble_members(exp)

for exp in experiments:
    print(exp)
    means_list1 = []
    keys_list = []
    variables_dir = identify_ensemble_members("piControl") # set to exp for all experiments

    for var in variables_dir:    
        time_period_np = list(variables_dir["nao_pattern_djf"].values())
        time_period_mean = np.nanmean(time_period_np, axis=0)
        #print(time_period_mean.shape)
        time_period_np.append(time_period_mean)
        #np.append(time_period_np, [time_period_mean], axis=0)
        keys = list(variables_dir["nao_pattern_djf"].keys())
        print(keys)
        keys.append("Ensemble Mean")

        means_list1.append(time_period_np)
        keys_list.append(keys)

        #break
    #break

    draw_models(means_list1[0], keys_list[0], "nao" + exp, lons, lats, 5, 5, min_val=-5, max_val=5) # nao_pattern
    #draw_models(means_list1[1], keys_list[1], "psl" + exp, lons, lats, 5, 5, min_val=980, max_val=1040)
    #draw_models(means_list1[2], keys_list[2], "tas" + exp, lons, lats, 5, 5, min_val=-40, max_val=30)
    #draw_models(means_list1[3], keys_list[3], "pr" + exp, lons, lats, 5, 5, min_val=0, max_val=15)
    break


# # ONLY ENSEMBLE MEANS OF VARIABLES

# In[9]:


def calculate_means(list_of_numbers):
    n=0
    average=0
    for element in list_of_numbers:
        average=(n*average+element)/(n+1)
        n=n+1
    return average


# In[10]:


def identify_ensemble_means(exp, variable, specific_exp="all"):  ## istället för specific_exp --> namn på experiment

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    time_dir = {}
    mean_dir = {}
    # r = root, d = directories, f = files
    for name in exp:
        time_dir[name]=[]
    
    grid_out = xe.util.grid_global(1,1)
    
    for r, d, f in os.walk(path):
        for file in f:
                
            pathtofile = os.path.join(r, file)
            if not file.startswith("."):
                
                open_file = xr.open_dataset(pathtofile, decode_times=False)

                regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]

            
                if specific_exp != "all":
                    if specific_exp in pathtofile:
                        temp = open_file.variables.get(variable)
                        if temp is not None:
                            data = regridder(temp.values)
                            time_dir[specific_exp].append(data)
                else:
                    for name in exp:
                        if name in pathtofile:
                            temp = open_file.variables.get(variable)
                            if temp is not None:
                                data = regridder(temp.values)
                                time_dir[name].append(data)            
                            break
                open_file.close()
            
    for experiment in time_dir.keys():
        if len(time_dir[experiment]) != 0:
            mean_dir[experiment]=calculate_means(time_dir[experiment])
            
    return mean_dir


# In[46]:


### PLOT ONLY ENSEMBLE MEANS ### 

# UPDATE

nao_tas_test = identify_ensemble_means(experiments, "nao_tas_regression_djf", "piControl") # add specific experiment if needed, e.g. "piControl"


# In[74]:


draw_models(nao_tas_test.values(), " ", "nao_tas_regression_ens_mean_pi", lons, lats, 1, 1, min_val=-3.5, max_val=3.5)


# In[57]:


# UPDATE

nao_pr_test = identify_ensemble_means(experiments, "nao_pr_regression_djf", "piControl") # add specific experiment if needed, e.g. "piControl"


# In[70]:


draw_models(nao_pr_test.values(), " ", "nao_pr_regression_ens_mean_pi", lons, lats, 1, 1, min_val=-1.2, max_val=1.2)


# In[137]:


# FIGURE 1 IN DISCUSSION - UPDATE

# NAO ensemble means in all experiments

all_means_nao = identify_ensemble_means(experiments, "nao_pattern_djf")


# In[144]:



draw_models(all_means_nao.values(), all_means_nao.keys(), "nao_all_means", lons, lats, 1, 5, min_val=-5, max_val=5)


# # ENSEMBLE MEANS FOR CMIP5 AND CMIP6 SEPARATELY

# In[11]:


def identify_CMIP_ensemble_means(models, variable, name, specific_exp="all"):  ## istället för specific_exp --> namn på experiment

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    model_list = []
    mean_dir = {}
    # r = root, d = directories, f = files
    
    grid_out = xe.util.grid_global(1,1)
    
    for r, d, f in os.walk(path):
        for file in f:
                
            pathtofile = os.path.join(r, file)
            if not file.startswith("."):
                
                open_file = xr.open_dataset(pathtofile, decode_times=False)

                regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]
                
                index = pathtofile.rfind("_")
                    
                model_name = pathtofile[:index]

                if specific_exp != "all":
                    if specific_exp in pathtofile and (len(list(filter(lambda x: model_name in x,models)))>0):
                        temp = open_file.variables.get(variable)
                        if temp is not None:
                            data = regridder(temp.values)
                            model_list.append(data)
            
                open_file.close()
    

    mean_dir[name]=calculate_means(model_list) # can make into a variable 
            
    return mean_dir


# In[68]:


CMIP6_models = ["AWI-CM", "HADGEM3-GC31", "IPSL-CM6A-LR",
                "BCC-CSM2-MR", "CNRM-CM6", "CNRM-ESM2", "GISS-E2-1-G",
                "MIROC6", "MRI-ESM2"]
CMIP5_models = ["BCC-CSM1","CCSM4", "CNRM-CM5","COSMOS-ASO","CSIRO-Mk3", "FGOALS-g2", "GISS-E2-R", "IPSL-CM5A", "MIROC-ESM", "MPI-ESM-P", "MRI-CGCM3"]

CMIP6_average = identify_CMIP_ensemble_means(CMIP6_models, "nao_pattern_djf", "CMIP6 MEAN", "piControl")
CMIP6_tas= identify_CMIP_ensemble_means(CMIP6_models, "nao_tas_regression_djf", " ", "piControl")
CMIP6_pr = identify_CMIP_ensemble_means(CMIP6_models, "nao_pr_regression_djf", " ", "piControl")

CMIP5_average = identify_CMIP_ensemble_means(CMIP5_models, "nao_pattern_djf", "CMIP5 MEAN", "piControl")
CMIP5_tas = identify_CMIP_ensemble_means(CMIP5_models, "nao_tas_regression_djf", " ", "piControl")
CMIP5_pr = identify_CMIP_ensemble_means(CMIP5_models, "nao_pr_regression_djf", " ", "piControl")


# In[92]:


#print(CMIP6_tas["CMIP6_MEAN"].shape)
draw_models_one(CMIP6_pr.values(), CMIP6_pr.keys(), "CMIP6_pr_average", lons, lats, 1, 1, -1.8, 1.8)


# In[94]:


#print(CMIP5_average["CMIP5_MEAN"].shape)
draw_models_one(CMIP5_pr.values(), CMIP5_pr.keys(), "CMIP5_pr_average", lons, lats, 1, 1, -1.8, 1.8)


# In[131]:


### DIFFERENCES BETWEEN CMIP5 AND CMIP6 ###

#NAO
NAO_CMIP5 = CMIP5_average["CMIP5 MEAN"]
NAO_CMIP6 = CMIP6_average["CMIP6 MEAN"]
CMIP_NAO_DIFF = NAO_CMIP6-NAO_CMIP5

#TAS
TAS_CMIP5 = CMIP5_tas[" "]
TAS_CMIP6 = CMIP6_tas[" "]
CMIP_TAS_DIFF = TAS_CMIP6-TAS_CMIP5

#PR
PR_CMIP5 = CMIP5_pr[" "]
PR_CMIP6 = CMIP6_pr[" "]
CMIP_PR_DIFF = PR_CMIP6-PR_CMIP5


draw_models_one([CMIP_NAO_DIFF], [" "], "CMIP6-CMIP5_NAO", lons, lats, 1,1,-1.2,1.2)
#draw_models_one([CMIP_TAS_DIFF], [" "], "CMIP6-CMIP5_TAS", lons, lats, 1,1,-0.5,0.5)
#draw_models_one([CMIP_PR_DIFF], [" "], "CMIP6-CMIP5_PR", lons, lats, 1,1,-0.5,0.5)


# # all models regression

# In[86]:


### UPDATE IF NEEDED (APPENDIX?) ###
lats= np.arange(-90,90,1)
lons= np.arange(-180,180,1)

#var_dir["nao_pr_regression_djf"]

for a in experiments:
    print(a)
    variables_dir = identify_ensemble_members(a, False)
    draw_models(variables_dir["nao_pr_regression_djf"].values(), variables_dir["nao_pr_regression_djf"].keys(), "moimoi3", lons, lats, 5, 5)
    draw_models(variables_dir["nao_tas_regression_djf"].values(), variables_dir["nao_tas_regression_djf"].keys(), "moimoi4", lons, lats, 5, 5)


# # identify difference field between ensemble members

# In[12]:


def identify_ensemble_difference_field(experiments, not_regression=True):

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    variables = ["nao_pattern_djf", "psl_spatialmean_djf", "tas_spatialmean_djf", "pr_spatialmean_djf", "nao_tas_regression_djf", "nao_pr_regression_djf"]
    #models = ["AWI-CM", "HADGEM3-G", "IPSL-CM6A-LR", "BCC-CSM2-MR", "CNRM-CM6", "CNRM-ESM2", "GISS-E2-1-G", "MIROC6", "MRI-ESM2","CCSM4","COSMOS-ASO","CSIRO-Mk3L","EC-EARTH-2-2", "FGOALS-g2","MRI-CGCM3","MIROC-ESM","GISS-E2-R","CNRM-CM5","MPI-ESM-P"]
    models2 = []
    
    var_dir = {}
    # r = root, d = directories, f = files
    nao= {}
    psl= {}
    pr= {}
    tas= {}
    nao_pr = {}
    nao_tas = {}
    nao_time = {}
    exp_dir = {}
    result_dir = {}
    
    grid_out = xe.util.grid_global(1,1)
    
    for exp in experiments:
        
        for r, d, f in os.walk(path):
            for file in f:
                if exp in file:

                    pathtofile = os.path.join(r, file)

                    open_file = xr.open_dataset(pathtofile, decode_times=False)

                    #regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)

                    pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                    pathtofile = pathtofile[1:-5]

                    print(pathtofile)
                    
                    index = pathtofile.rfind("_")
                    
                    model_name = pathtofile[:index]
                    
                    if model_name not in models2:
                        models2.append(model_name)

                    if not_regression:

                        temp = open_file.variables.get("nao_pattern_djf")
                        if temp is not None:
                            data = open_file#.variables.get("nao_pattern_djf")
                            data = reverse_data(data, file)
                            nao[pathtofile] = data

                        temp = open_file.variables.get("psl_spatialmean_djf")
                        if temp is not None:
                            data = open_file#.variables.get("psl_spatialmean_djf")
                            data = reverse_data(data, file)
                            psl[pathtofile]=data     

                        temp = open_file.variables.get("tas_spatialmean_djf")
                        if temp is not None:
                            data = open_file#.variables.get("tas_spatialmean_djf")
                            data = reverse_data(data, file)
                            tas[pathtofile]=data

                        temp = open_file.variables.get("pr_spatialmean_djf")
                        if temp is not None:
                            data = open_file#.variables.get("pr_spatialmean_djf")
                            data = reverse_data(data, file)
                            pr[pathtofile]= data

                    else:

                        temp = open_file.variables.get("nao_pr_regression_djf")
                        if temp is not None:
                            data = open_file#.variables.get("nao_pr_regression_djf")
                            data = reverse_data(data, file)
                            nao_pr[pathtofile]= data

                        temp = open_file.variables.get("nao_tas_regression_djf")
                        if temp is not None:
                            data = open_file#.variables.get("nao_tas_regression_djf")
                            data = reverse_data(data, file)
                            nao_tas[pathtofile]= data
                            
                            
        if not_regression:
            var_dir["nao_pattern_djf"]=nao
            var_dir["psl_spatialmean_djf"]=psl
            var_dir["tas_spatialmean_djf"]=tas
            var_dir["pr_spatialmean_djf"]=pr

        else:
            var_dir["nao_pr_regression_djf"]=nao_pr
            var_dir["nao_tas_regression_djf"]=nao_tas
    
    
        exp_dir[exp]=var_dir
        
    
    print(models2)
    
    for timeperiod in exp_dir.keys():
        diff_dir = {}
        if timeperiod != "piControl":
            for var in variables:
                var_dir = {}
                for model in models2:
                    left=exp_dir.get(timeperiod).get(var)
                    right=exp_dir.get("piControl").get(var)
                    if left is not None and right is not None:
                        left = left.get(model +"_"+ timeperiod)
                        right = right.get(model +"_"+ "piControl")
                        if left is not None and right is not None:
                            diff =left[var]-right[var]
                            regridder = xe.Regridder(left, grid_out, 'nearest_s2d', reuse_weights=True)
                            var_dir[var+model]=regridder(diff)
                            
                diff_dir[var]=var_dir
                            
        result_dir[timeperiod]=diff_dir
    
    
    return result_dir


# In[82]:


regrid_diff = identify_ensemble_difference_field(experiments)


# In[135]:


### difference fields with ensemble means ###

means_list = []
keys_list = []

for timeperiod in regrid_diff.keys():
    if timeperiod != "piControl": 
        for variable in regrid_diff[timeperiod].keys():
            n=0
            average=0
            if len(list(regrid_diff[timeperiod][variable].values())) != 0:
                for model in regrid_diff[timeperiod][variable].keys():
                    average=(n*average+regrid_diff[timeperiod][variable][model])/(n+1)
                    n=n+1
                time_period_np = list(regrid_diff[timeperiod][variable].values())
                #time_period_mean = np.nanmean(time_period_np, axis=0)
                time_period_np.append(average)
                keys = list(regrid_diff[timeperiod][variable].keys())
                keys.append(variable+"_mean")
                draw_models(time_period_np, keys, "regrid_diff_and_ensemble", lons, lats, 4, 5, -10, 10)


# # specific difference

# In[25]:


def identify_ensemble_difference_field_specific(experiment, variable):

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    variables = ["nao_pattern_djf", "psl_spatialmean_djf", "tas_spatialmean_djf", "pr_spatialmean_djf", "nao_tas_regression_djf", "nao_pr_regression_djf"]
    #models = ["AWI-ESM", "HADGEM3-G", "IPSL-CM6A-LR", "BCC-CSM2-MR", "CNRM-CM6", "CNRM-ESM2", "GISS-E2-1-G", "MIROC6", "MRI-ESM2","CCSM4","COSMOS-ASO","CSIRO-Mk3L","EC-EARTH-2-2", "FGOALS-g2","MRI-CGCM3","MIROC-ESM","GISS-E2-R","CNRM-CM5","MPI-ESM-P"]
    models2 = []
    
    var_dir = {}
    exp_dir = {}
    # r = root, d = directories, f = files
    nao= {}
    pi = {}
    
    grid_out = xe.util.grid_global(1,1)
    
    for r, d, f in os.walk(path):
        for file in f:
            if experiment in file:

                pathtofile = os.path.join(r, file)

                open_file = xr.open_dataset(pathtofile, decode_times=False)

                #regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]
                
                #pathtofile = pathtofile.replace("_", " ")

                print(pathtofile)

                index = pathtofile.rfind("_")

                model_name = pathtofile[:index]


                if model_name not in models2:
                    models2.append(model_name)


                temp = open_file.variables.get(variable)
                if temp is not None:
                    data = open_file#.variables.get("nao_pattern_djf")
                    data = reverse_data(data, file)
                    nao[model_name] = data

            if "piControl" in file:

                pathtofile = os.path.join(r, file)

                open_file = xr.open_dataset(pathtofile, decode_times=False)

                #regridder = xe.Regridder(open_file, grid_out, 'nearest_s2d', reuse_weights=False)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]
                
                #pathtofile = pathtofile.replace("_", " ")

                print(pathtofile)

                index = pathtofile.rfind("_")

                model_name = pathtofile[:index]


                if model_name not in models2:
                    models2.append(model_name)


                temp = open_file.variables.get(variable)
                if temp is not None:
                    data = open_file#.variables.get("nao_pattern_djf")
                    data = reverse_data(data, file)
                    pi[model_name] = data
                    
               
    exp_dir[experiment]=nao        
    exp_dir["piControl"]=pi

    diff_dir = {}
    
    for model in nao.keys():
        left=exp_dir.get(experiment).get(model)
        right=exp_dir.get("piControl").get(model)
        
        if left is not None and right is not None:
            print(model)
            print(variable)
            print(left[variable].shape)
            print(right[variable].shape)
            diff =left[variable]-right[variable]
            regridder = xe.Regridder(left, grid_out, 'nearest_s2d', reuse_weights=True)
            diff_dir[model]=regridder(diff)

    return diff_dir


# In[22]:


variables = ["nao_pattern_djf", "psl_spatialmean_djf", "tas_spatialmean_djf", "pr_spatialmean_djf", "nao_tas_regression_djf", "nao_pr_regression_djf"]
array = [1,2]

i=0

for var in variables:
    for exp in experiments:
        if exp != "piControl":
            diff_lig_nao_new = identify_ensemble_difference_field_specific(exp, var)
            n=0
            average=0

            for model in diff_lig_nao_new.keys():

                if len(list(diff_lig_nao_new[model])) != 0:
                    #print(diff_lig_nao[model])
                    average=(n*average+diff_lig_nao_new[model])/(n+1)
                    n=n+1

            time_period_np = list(diff_lig_nao_new.values())
            #time_period_mean = np.nanmean(time_period_np, axis=0)
            time_period_np.append(average)
            keys = list(diff_lig_nao_new.keys())
            keys.append(var+"_mean")
            #draw_models(time_period_np, keys, "regrid_diff_and_ensemble", lons, lats, 4, 5, array[i], array[i+1])
            i=i+1


# # ensemble mean differences

# In[146]:


def ensemble_difference_mean(var, experiments):
    
    exp_mean_dir = {}

    array = [1,2]

    for exp in experiments:
        if exp != "piControl":
            differences = identify_ensemble_difference_field_specific(exp, var)
            n=0
            average=0

            for model in differences.keys():

                if len(list(differences[model])) != 0:
                    #print(diff_lig_nao[model])
                    average=(n*average+differences[model])/(n+1)
                    n=n+1
            
            exp_mean_dir[exp]=average
            
    return exp_mean_dir


# # ENS MEAN FOR ALL VARIABLES FOR THE DISCUSSION # UPDATE # PERSE

# In[ ]:


# NAO

nao_pattern_difference = ensemble_difference_mean("nao_pattern_djf", experiments)
nao_pattern_lig = ensemble_difference_mean("nao_pattern_djf", ["lig127k"])
nao_pattern_lgm = ensemble_difference_mean("nao_pattern_djf", ["lgm"]) 
nao_pattern_mh = ensemble_difference_mean("nao_pattern_djf", ["midHolocene"]) 
nao_pattern_1pct = ensemble_difference_mean("nao_pattern_djf", ["1pctCO2"]) 


# In[ ]:


# TAS REGRESSION

tas_pattern_difference = ensemble_difference_mean("nao_tas_regression_djf", experiments)
tas_pattern_lig = ensemble_difference_mean("nao_tas_regression_djf", ["lig127k"])
tas_pattern_lgm = ensemble_difference_mean("nao_tas_regression_djf", ["lgm"])
tas_pattern_mh = ensemble_difference_mean("nao_tas_regression_djf", ["midHolocene"])
tas_pattern_1pct = ensemble_difference_mean("nao_tas_regression_djf", ["1pctCO2"])


# In[ ]:


# PR REGRESSION

pr_pattern_difference = ensemble_difference_mean("nao_pr_regression_djf", experiments) 
pr_pattern_lig = ensemble_difference_mean("nao_pr_regression_djf", ["lig127k"]) 
pr_pattern_lgm = ensemble_difference_mean("nao_pr_regression_djf", ["lgm"]) 
pr_pattern_mh = ensemble_difference_mean("nao_pr_regression_djf", ["midHolocene"]) 
pr_pattern_1pct = ensemble_difference_mean("nao_pr_regression_djf", ["1pctCO2"]) 


# In[214]:


# UPDATE-here - DO NOT RUN AGAIN IF NOT NECESSARY-TAKES AGES!! - MEAN STATE CLIMATOLOGY PLOTS - TAS, PR, PSL

#TAS
tas_mean_difference = ensemble_difference_mean("tas_spatialmean_djf", experiments)
tas_mean_diff_lig = ensemble_difference_mean("tas_spatialmean_djf", ["lig127k"])
tas_mean_diff_lgm = ensemble_difference_mean("tas_spatialmean_djf", ["lgm"])
tas_mean_diff_mh = ensemble_difference_mean("tas_spatialmean_djf", ["midHolocene"])
tas_mean_diff_1pct = ensemble_difference_mean("tas_spatialmean_djf", ["1pctCO2"])


# In[265]:


#PR
pr_mean_difference = ensemble_difference_mean("pr_spatialmean_djf", experiments)
pr_mean_diff_lig = ensemble_difference_mean("pr_spatialmean_djf", ["lig127k"])
pr_mean_diff_lgm = ensemble_difference_mean("pr_spatialmean_djf", ["lgm"])
pr_mean_diff_mh = ensemble_difference_mean("pr_spatialmean_djf", ["midHolocene"])
pr_mean_diff_1pct = ensemble_difference_mean("pr_spatialmean_djf", ["1pctCO2"])


# In[276]:


#PSL
psl_mean_difference = ensemble_difference_mean("psl_spatialmean_djf", experiments)
psl_mean_diff_lig = ensemble_difference_mean("psl_spatialmean_djf", ["lig127k"])
psl_mean_diff_lgm = ensemble_difference_mean("psl_spatialmean_djf", ["lgm"])
psl_mean_diff_mh = ensemble_difference_mean("psl_spatialmean_djf", ["midHolocene"])
psl_mean_diff_1pct = ensemble_difference_mean("psl_spatialmean_djf", ["1pctCO2"])


# In[294]:


# DRAW MEAN STATE CLIMATOLOGY

# NAO 
#draw_models(nao_pattern_difference.values(), nao_pattern_difference.keys(), "tas_diff_mean", lons, lats, 1,4,-11.7,11.7)
#draw_models_one(nao_pattern_lig.values(), nao_pattern_lig.keys(), "nao_diff_mean_lig", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(nao_pattern_lgm.values(), nao_pattern_lgm.keys(), "nao_diff_mean_lgm", lons, lats, 1,1,-34.7,34.7)
#draw_models_one(nao_pattern_mh.values(), nao_pattern_mh.keys(), "nao_diff_mean_mh", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(nao_pattern_1pct.values(), nao_pattern_1pct.keys(), "nao_diff_mean_1pct", lons, lats, 1,1,-14.7,14.7)


# TAS REGRESSION
#draw_models(tas_pattern_difference.values(), tas_pattern_difference.keys(), "tas_diff_mean", lons, lats, 1,4,-11.7,11.7)
#draw_models_one(tas_pattern_lig.values(), tas_pattern_lig.keys(), "tas_re_diff_mean_lig", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(tas_pattern_lgm.values(), tas_pattern_lgm.keys(), "tas_re_diff_mean_lgm", lons, lats, 1,1,-34.7,34.7)
#draw_models_one(tas_pattern_mh.values(), tas_pattern_mh.keys(), "tas_re_diff_mean_mh", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(tas_pattern_1pct.values(), tas_pattern_1pct.keys(), "tas_re_diff_mean_1pct", lons, lats, 1,1,-14.7,14.7)


# PR REGRESSION
#draw_models(pr_pattern_difference.values(), tas_mean_difference.keys(), "tas_diff_mean", lons, lats, 1,4,-11.7,11.7)
#draw_models_one(pr_pattern_lig.values(), pr_pattern_lig.keys(), "pr_re_diff_mean_lig", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(pr_pattern_lgm.values(), pr_pattern_lgm.keys(), "pr_re_diff_mean_lgm", lons, lats, 1,1,-34.7,34.7)
#draw_models_one(pr_pattern_mh.values(), pr_pattern_mh.keys(), "pr_re_diff_mean_mh", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(pr_pattern_1pct.values(), pr_pattern_1pct.keys(), "pr_re_diff_mean_1pct", lons, lats, 1,1,-14.7,14.7)



# TAS (keep coordinates for these, letters: a,d,g,j, color: bwr, %0f) 
#draw_models(tas_mean_difference.values(), tas_mean_difference.keys(), "tas_diff_mean", lons, lats, 1,4,-11.7,11.7)
#draw_models_one(tas_mean_diff_lig.values(), tas_mean_diff_lig.keys(), "tas_diff_mean_lig", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(tas_mean_diff_lgm.values(), tas_mean_diff_lgm.keys(), "tas_diff_mean_lgm", lons, lats, 1,1,-34.7,34.7)
#draw_models_one(tas_mean_diff_mh.values(), tas_mean_diff_mh.keys(), "tas_diff_mean_mh", lons, lats, 1,1,-5.7,5.7)
#draw_models_one(tas_mean_diff_1pct.values(), tas_mean_diff_1pct.keys(), "tas_diff_mean_1pct", lons, lats, 1,1,-14.7,14.7)


# PR (no coordinates, letters: b,e,h,k, color:, %1f)
#draw_models(pr_mean_difference.values(), pr_mean_difference.keys(), "pr_diff_mean", lons, lats, 1,4,-2,2)
#draw_models_one(pr_mean_diff_lig.values(), pr_mean_diff_lig.keys(), "pr_diff_mean_lig", lons, lats, 1,1,-2,2)
#draw_models_one(pr_mean_diff_lgm.values(), pr_mean_diff_lgm.keys(), "pr_diff_mean_lgm", lons, lats, 1,1,-2,2)
#draw_models_one(pr_mean_diff_mh.values(), pr_mean_diff_mh.keys(), "pr_diff_mean_mh", lons, lats, 1,1,-2,2)
#draw_models_one(pr_mean_diff_1pct.values(), pr_mean_diff_1pct.keys(), "pr_diff_mean_1pct", lons, lats, 1,1,-2,2)

# PSL (no coorfinates, letters: c,f,i,l, color, %0f)
#draw_models(psl_mean_difference.values(), psl_mean_difference.keys(), "psl_diff_mean", lons, lats, 1,4,-15,15)
#draw_models_one(psl_mean_diff_lig.values(), psl_mean_diff_lig.keys(), "psl_diff_mean_lig", lons, lats, 1,1,-5.8,5.8)
#draw_models_one(psl_mean_diff_lgm.values(), psl_mean_diff_lgm.keys(), "psl_diff_mean_lgm", lons, lats, 1,1,-35,35)
#draw_models_one(psl_mean_diff_mh.values(), psl_mean_diff_mh.keys(), "psl_diff_mean_mh", lons, lats, 1,1,-5.8,5.8)
draw_models_one(psl_mean_diff_1pct.values(), psl_mean_diff_1pct.keys(), "psl_diff_mean_1pct", lons, lats, 1,1,-5.8,5.8)


# # DIFFERENCE OF ALL EXPERIMENTS

# In[38]:


# LIG
diff_lig_nao = identify_ensemble_difference_field_specific("lig127", "nao_pattern_djf")
diff_lig_tas_nao = identify_ensemble_difference_field_specific("lig127", "nao_tas_regression_djf")
diff_lig_pr_nao = identify_ensemble_difference_field_specific("lig127", "nao_pr_regression_djf")


# In[57]:


# LGM
diff_lgm_nao = identify_ensemble_difference_field_specific("lgm", "nao_pattern_djf")
diff_lgm_tas_nao = identify_ensemble_difference_field_specific("lgm", "nao_tas_regression_djf")
diff_lgm_pr_nao = identify_ensemble_difference_field_specific("lgm", "nao_pr_regression_djf")


# In[27]:


#MH
diff_mh_nao = identify_ensemble_difference_field_specific("midHolocene", "nao_pattern_djf")
diff_mh_tas_nao = identify_ensemble_difference_field_specific("midHolocene", "nao_tas_regression_djf")
diff_mh_pr_nao = identify_ensemble_difference_field_specific("midHolocene", "nao_pr_regression_djf")


# In[38]:


#1PCT
diff_1pct_nao = identify_ensemble_difference_field_specific("1pctCO2", "nao_pattern_djf")
diff_1pct_tas_nao = identify_ensemble_difference_field_specific("1pctCO2", "nao_tas_regression_djf")
diff_1pct_pr_nao = identify_ensemble_difference_field_specific("1pctCO2", "nao_pr_regression_djf")


# # LIG127K

# In[105]:


list_lig = ["HadGEM3", "IPSL-CM6A", "AWI-CM", "Ensemble Mean"]


# In[42]:


### DIFF_LIG_NAO ###
n=0
average=0

list_lig = ["HadGEM3", "IPSL-CM6A", "AWI-CM", "Ensemble Mean"]

for model in diff_lig_nao.keys():
    
    if len(list(diff_lig_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lig_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lig_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lig_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, list_lig, "lig127k_regrid_diff_and_ensemble", lons, lats, 1, 4, -2, 2)


# In[111]:


### DIFF_LIG_TAS_NAO ###
n=0
average=0

for model in diff_lig_tas_nao.keys():
    
    if len(list(diff_lig_tas_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lig_tas_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lig_tas_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lig_tas_nao.keys())
keys.append("nao_pattern_djf"+"_mean")
draw_models(time_period_np, list_lig, "lig127k_nao_tas_diff_and_ensemble", lons, lats, 1, 4, -1, 1)


# In[109]:


### DIFF_LIG_PR_NAO ###
n=0
average=0

for model in diff_lig_pr_nao.keys():
    
    if len(list(diff_lig_pr_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lig_pr_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lig_pr_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lig_pr_nao.keys())
keys.append("nao_pattern_djf"+"_mean")
draw_models(time_period_np, list_lig, "lig127k_nao_pr_diff_and_ensemble", lons, lats, 1, 4, -0.5, 0.5)


# # LGM

# In[78]:


n=0
average=0

list_lgm = ["AWI-CM", "GISS-E2-R", "MIROC-ESM", "COSMOS", "FGOALS-g2", "MPI-ESM-P", "CNRM-CM5", "IPSL-CM5A", "MRI-CGCM3", "CCSM4", "Ensemble Mean"]

for model in diff_lgm_nao.keys():
    
    if len(list(diff_lgm_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lgm_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lgm_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lgm_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, list_lgm, "lgm_regrid_diff_and_ensemble", lons, lats, 3, 4, -5.8, 5.8)


# In[123]:


### DIFF_LGM_TAS_NAO ### CHANGE BACK TO ORIGINAL
n=0
average=0

for model in diff_lgm_tas_nao.keys():
    
    if len(list(diff_lgm_tas_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lgm_tas_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lgm_tas_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lgm_tas_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, list_lgm, "lgm_nao_tas_diff_and_ensemble", lons, lats, 3, 4, -3, 3)


# In[103]:


### DIFF_LGM_PR_NAO ###
n=0
average=0

for model in diff_lgm_pr_nao.keys():
    
    if len(list(diff_lgm_pr_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_lgm_pr_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_lgm_pr_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_lgm_pr_nao.keys())
keys.append("nao_pattern_djf"+"_mean")
draw_models(time_period_np, list_lgm, "lgm_nao_pr_diff_and_ensemble", lons, lats, 3, 4, -1.2, 1.2)


# # MIDHOLOCENE

# In[28]:


n=0
average=0

list_midHolocene = ["IPSL-CM5A","IPSL-CM6A",
                    "HadGEM3","GISS-E2-R",
                    "CNRM-CM5", "CCSM4","FGOALS-g2",
                    "CSIRO-Mk3", "MPI-ESM-P", "BCC-CSM1",
                    "AWI-CM","MRI-CGCM3","MIROC-ESM","Ensemble Mean"]

for model in diff_mh_nao.keys():
    
    if len(list(diff_mh_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_mh_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_mh_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_mh_nao.keys())
keys.append("Ensemble Mean")
print(keys)
draw_models(time_period_np, keys, "mh_regrid_diff_and_ensemble", lons, lats, 4, 4, -2, 2)


# In[31]:


n=0
average=0

for model in diff_mh_tas_nao.keys():
    
    if len(list(diff_mh_tas_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_mh_tas_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_mh_tas_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_mh_tas_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, keys, "mh_nao_tas_diff_and_ensemble", lons, lats, 4, 4, -1, 1)


# In[37]:


n=0
average=0

for model in diff_mh_pr_nao.keys():
    
    if len(list(diff_mh_pr_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_mh_pr_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_mh_pr_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_mh_pr_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, keys, "mh_nao_pr_diff_and_ensemble", lons, lats, 4, 4, -0.5, 0.5)


# ## 1PCTCO2

# In[47]:


n=0
average=0

for model in diff_1pct_nao.keys():
    
    if len(list(diff_1pct_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_1pct_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_1pct_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_1pct_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, keys, "1pct_regrid_diff_and_ensemble", lons, lats, 4, 5, -2, 2)


# In[45]:


# 1PCT NAO TAS DIFF #
n=0
average=0

for model in diff_1pct_tas_nao.keys():
    
    if len(list(diff_1pct_tas_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_1pct_tas_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_1pct_tas_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_1pct_tas_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, keys, "1pct_nao_tas_diff_and_ensemble", lons, lats, 4, 5, -2, 2)


# In[43]:


# 1PCT NAO PR DIFF #
n=0
average=0

for model in diff_1pct_pr_nao.keys():
    
    if len(list(diff_1pct_pr_nao[model])) != 0:
        #print(diff_lig_nao[model])
        average=(n*average+diff_1pct_pr_nao[model])/(n+1)
        n=n+1
        
time_period_np = list(diff_1pct_pr_nao.values())
#time_period_mean = np.nanmean(time_period_np, axis=0)
time_period_np.append(average)
keys = list(diff_1pct_pr_nao.keys())
keys.append("Ensemble Mean")
draw_models(time_period_np, keys, "1pct_nao_pr_diff_and_ensemble", lons, lats, 4, 5, -0.5, 0.5)


# # all ensembles difference fields in each experiment together

# In[104]:


"""identify difference field for ensemble members"""

def identify_ensemble_difference_field(means, keys):
    print(keys.index('piControl'))
    pi_index = keys.index('piControl')
    
    ens_difference = []
    
    for variable_list in means:
        var_diff = []
        for i in range(len(variable_list)):
            if i != pi_index:
                var_diff.append(variable_list[pi_index]-variable_list[i])
        ens_difference.append(var_diff)
        
    return ens_difference


# In[128]:


nao = []
psl = []
tas = []
pr = []

for a in experiments:
    print(a)
    variables_dir = identify_ensemble_members(a)

    ensembles_means = []
    
    print("MOIMOIMOMOIMOIMOIMOIM")

    for time_period in variables_dir:
        
        time_period_np = np.array(list(variables_dir[time_period].values()))
        time_period_mean = np.nanmean(time_period_np, axis=0)
        if "nao_pattern_djf" == time_period:
            nao.append(time_period_mean) 
        elif "psl_spatialmean_djf" == time_period:
            psl.append(time_period_mean)
        elif "tas_spatialmean_djf" == time_period:
            tas.append(time_period_mean)
        elif "pr_spatialmean_djf" == time_period:
            pr.append(time_period_mean)

var_list = [nao,psl,tas,pr]

a = identify_ensemble_difference_field(var_list, experiments) # change experiments to another list of names

draw_models(a[0], experiments, "nao_all", lons, lats, 1, 4, min_val=-1, max_val=1)
draw_models(a[1], experiments, "psl_all", lons, lats, 1, 4, min_val=-20,max_val=20)
draw_models(a[2], experiments, "tas_all", lons, lats, 1, 4, min_val=-15, max_val=15)
draw_models(a[3], experiments, "pr_all", lons, lats, 1, 4, min_val=-2,max_val=2)


# # scatter plots

# In[20]:


ipcc_dataset = xr.open_dataset('ipcc_region_array-2.nc')
ipcc = ipcc_dataset.variables["ipcc_region_array"]


#Populate some dictionaries for searching over
ts_dict={}
for i in range(ipcc.sizes['timeseries_id']):
        ts_dict[ipcc_dataset['timeseries'].data[i]]=ipcc_dataset['timeseries_id'].data[i]
gcm_dict={}
for i in range(ipcc.sizes['model_id']):
        gcm_dict[ipcc_dataset['models'].data[i]]=ipcc_dataset['model_id'].data[i]
expt_dict={}
for i in range(ipcc.sizes['expt_id']):
        expt_dict[ipcc_dataset['expts'].data[i]]=ipcc_dataset['expt_id'].data[i]
print(gcm_dict)


# In[21]:


print(ipcc[:].dims)


# In[40]:


#Try loading in some sample data
x_variable="ipcc_NEU_pr"
y_variable="nao_pattern_djf"
x_var_is_stddev=False
y_var_is_stddev=False

piCtl=ipcc[:,0,ts_dict[x_variable],x_var_is_stddev].data
others=ipcc[:,:,ts_dict[x_variable],x_var_is_stddev].data
print(others[:,2])
print(others[:,5])
x=others.T-piCtl.T
x.shape

piCtl=ipcc[:,0,ts_dict[y_variable],y_var_is_stddev].data
others=ipcc[:,:,ts_dict[y_variable],y_var_is_stddev].data
print(others[:,2])
print(others[:,5])
y=others.T-piCtl.T
y.shape


# In[41]:


colors = ['black','darkgreen','blue','lime','darkgreen','blue','lime','red','maroon','goldenrod','darkorange']
markers = ['o','+','+','+','o','o','o','o','o','o','o']

#Plot with lig for this one...
for expt in ipcc_dataset['expt_id'].data[1:]:
    this_label=ipcc_dataset['expts'].data[expt]
    #if this_label.startswith('lig') != True: 
    plt.scatter(x[expt,:], y[expt,:], color=colors[expt],marker=markers[expt],label=this_label)
    
plt.xlabel("mean "+x_variable)
plt.ylabel("Mean "+y_variable)
plt.legend(loc='best')


# In[42]:


modes_dataset = xr.open_dataset("modes_array.nc")
modes=modes_dataset.modes_array
print(modes_dataset)
#Populate some dictionaries for searching over
mode_dict={}
for i in range(modes_dataset.sizes['modes_id']):
        mode_dict[modes_dataset['modenames'].data[i]]=modes_dataset['modes_id'].data[i]
print(mode_dict)


# In[137]:


y_variable="nao_pattern_djf"
piCtl=modes[:,0,mode_dict[y_variable]].data
others=modes[:,:,mode_dict[y_variable]].data
print(others[:,2])
print(others[:,5])
y=others.T/piCtl.T
y.shape


# In[24]:


colors = ['black','darkgreen','blue','lime','darkgreen','blue','lime','red','maroon','goldenrod','darkorange']
markers = ['o','+','+','+','o','o','o','o','o','o','o']

#Plot with lig for this one...
for expt in ipcc_dataset['expt_id'].data[1:]:
    this_label=ipcc_dataset['expts'].data[expt]
    if this_label.startswith('lig') != True: 
        plt.scatter(x[expt,:], y[expt,:], color=colors[expt],marker=markers[expt],label=this_label)
    
plt.xlabel("mean "+x_variable)
plt.ylabel(y_variable+" Amplitude Change")
plt.legend(loc='best')


# # standard deviations for each nao_pattern in each model and each experiment - then, ensemble mean of all the std together

# In[228]:


for exp in experiments:
    print(exp)
    variables_dir = identify_ensemble_members(exp)
    
    for mod in variables_dir["nao_pattern_djf"].keys():
        print("\t",mod)
        nao = variables_dir["nao_pattern_djf"][mod][9:69,90:220] #[110:170,90:220]
        std = np.nanstd(nao)
        print("\t",std)


# In[231]:


variables_dir["nao_pattern_djf"]


# In[9]:


# ftp://wxmaps.org/pub/straus/CLIM_753/EOF.pdf
# https://www.nature.com/articles/nature12580#online-methods

for exp in experiments:
    print(exp)
    variables_dir = identify_ensemble_members(exp, "nao_pattern_djf")
    timeseries_dir = identify_ensemble_members(exp, "nao_timeseries_djf", True,False)
    
    for mod, mod2 in zip(variables_dir["nao_pattern_djf"].keys(),timeseries_dir["nao_timeseries_djf"].keys()):   # the leading spatial pattern is standardized by its respective spatial standard deviation
        print("\t",mod)
        nao = variables_dir["nao_pattern_djf"][mod][9:69,90:220] #[110:170,90:220] - remember changing this depending on the reverse_data
        nao_time = timeseries_dir["nao_timeseries_djf"][mod2]
        nao_std = np.nanstd(nao)
        nao_mean = np.nanmean(nao)
        new_EOF = (nao-nao_mean)/nao_std
        new_nao_time = nao_time*nao_std
        std_nao_time = np.nanstd(new_nao_time)
        print("\t",mod2)
        print("\t",std_nao_time)


# # make your own scatter plots

# In[83]:


def identify_model_name():

    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"
    
    models2 = []
        
    for r, d, f in os.walk(path):
        for file in f:
            if "DS_Store" not in file and "AWI-CM_midHolocene" not in file:
                pathtofile = os.path.join(r, file)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]

                index = pathtofile.rfind("_")

                model_name = pathtofile[:index]

                if model_name not in models2 and "C20-Reanalysi" not in model_name:
                    models2.append(model_name)
    
    return models2


# In[84]:


modelss = identify_model_name()
print(modelss)


# In[149]:


def scatter_plots(experiments, var):
    
    variable_list = ["tas_spatialmean_djf", "pr_spatialmean_djf"]
    models_times_dict = {}
    difference_models_times = {}
    regrid_models = {}
    
    models = identify_model_name() 
    
    grid_out = xe.util.grid_global(1,1)
    
    for model in models:
        models_times_dict[model]= {}
        difference_models_times[model] = {}
        regrid_models[model] = {}
        
    path = "/Users/venniarra/Desktop/jupyter-notebook/test2"

    for r, d, f in os.walk(path):
        for file in f:
            if "DS_Store" not in file and "C20-Reanalysi" not in file and "AWI-CM_midHolocene" not in file and "HadGEM3_midHolocene" not in file:
                pathtofile = os.path.join(r, file)

                open_file = xr.open_dataset(pathtofile, decode_times=False)

                pathtofile = re.search('\/[^\/]+\.cvdp',pathtofile).group(0)

                pathtofile = pathtofile[1:-5]

                index = pathtofile.rfind("_")

                model_name = pathtofile[:index]
                
                pathtofile = pathtofile[index+1:]

                models_times_dict[model_name][pathtofile] = open_file
                
                #pathtofile = pathtofile.replace("_", " ")
    
    for model in models:
        for exp in models_times_dict[model].keys():
            if "piControl" in exp:
                piControl_name = exp
        for exp in models_times_dict[model].keys():
            if "piControl" not in exp:
                diff = models_times_dict[model][exp].variables.get(var) - models_times_dict[model][piControl_name].variables.get(var)
                regridder = xe.Regridder(models_times_dict[model][exp], grid_out, 'nearest_s2d', reuse_weights=False)
                new_grid = regridder(diff.values)
                difference_models_times[model][exp] = new_grid
                new_grid = regridder(models_times_dict[model][exp].variables.get(var).values)
                regrid_models[model][exp] = new_grid
    
    x_axis_tas_neu = {}
    x_axis_pr_neu = {}
    x_axis_tas_med = {}
    x_axis_pr_med = {}
    y_axis = []
    y_dict = {}
    x_dict = {}
    exp_dict = {}

    
    for exp in experiments:
        exp_dict[exp] = []
        x_axis_tas_neu[exp] = []
        print(x_axis_tas_neu.keys())
        x_axis_pr_neu[exp] = []
        x_axis_tas_med[exp] = []
        x_axis_pr_med[exp] = []
        print("hello")
        
    for model in models:
        y_axis = []
        for exp in difference_models_times[model].keys():
            nao_area = regrid_models[model][exp][9:69,90:220]
            nao_std  = np.nanstd(nao_area)
            nao_time = models_times_dict[model][exp].variables.get("nao_timeseries_djf")
            new_nao_time = nao_time*nao_std
            std_nao_time = np.nanstd(new_nao_time)
            y_axis.append(std_nao_time)
            
            
            current_exp = models_times_dict[model][exp]
            regridder = xe.Regridder(current_exp, grid_out, 'nearest_s2d', reuse_weights=False)
            
            for index, item in enumerate(variable_list):
                
                print(model)
                print(exp)
                print(item)
                
                new_grid = regridder(current_exp.variables.get(item).values)
                
    
                print(models_times_dict[model].keys())
        
                if item == "tas_spatialmean_djf":
                    pi_grid = regridder(models_times_dict[model]["piControl"].variables.get(item).values)
                    x_tas_neu=(np.nanmean(new_grid[19:39,170:215] - pi_grid[19:39,170:215]))
                    print(x_axis_tas_neu.keys())
                    x_axis_tas_neu[exp].append(x_tas_neu)
                
                if item == "pr_spatialmean_djf":
                    pi_grid = regridder(models_times_dict[model]["piControl"].variables.get(item).values)
                    x_pr_neu=(np.nanmean(new_grid[19:39,170:215] - pi_grid[19:39,170:215]))
                    x_axis_pr_neu[exp].append(x_pr_neu)
                    
                if item == "tas_spatialmean_djf":
                    pi_grid = regridder(models_times_dict[model]["piControl"].variables.get(item).values)
                    x_tas_med=(np.nanmean(new_grid[39:59,170:215] - pi_grid[39:59,170:215]))
                    x_axis_tas_med[exp].append(x_tas_med)
                if item == "pr_spatialmean_djf":
                    pi_grid = regridder(models_times_dict[model]["piControl"].variables.get(item).values)
                    x_pr_med=(np.nanmean(new_grid[39:59,170:215] - pi_grid[39:59,170:215]))
                    x_axis_pr_med[exp].append(x_pr_med)
            
            exp_dict[exp].append(std_nao_time)
    
            
    x_axis = [exp_dict, x_axis_tas_neu, x_axis_pr_neu, x_axis_tas_med, x_axis_pr_med] # ordningen för a
   
   
    return x_axis
    


# In[ ]:


a=scatter_plots(experiments, "nao_pattern_djf")


# In[ ]:


print(a[0].keys())


# In[145]:


# then make scatter plot with values from above

#make function (changing x-axis values)

colors = ['black','darkgreen','blue','lime','darkgreen','blue','lime','red','maroon','goldenrod','darkorange']
markers = ['o','+','+','+','o','o','o','o','o','o','o']

"""#Plot with lig for this one...
for expt in ipcc_dataset['expt_id'].data[1:]:
this_label=ipcc_dataset['expts'].data[expt]
#if this_label.startswith('lig') != True"""
i=0
for b in a[0].keys():
    if "piControl" not in b:
        plt.scatter(a[1][b], a[0][b], color=colors[i], marker=markers[i], label=b)
        plt.xlabel("neu_tas")
        plt.ylabel("nao_amplitude")
        plt.legend(loc='best')
        i=i+1


# In[ ]:



        
    """ tas_NEU = regrid_models[model][exp].variables.get("tas_spatialmean_djf")[19:39,170:215]
            pr_NEU = regrid_models[model][exp].variables.get("pr_timeseries_djf")[19:39,170:215]
            tas_MED = regrid_models[model][exp].variables.get("tas_spatialmean_djf")[39:59,170:215]
            pr_MED = regrid_models[model][exp].variables.get("pr_timeseries_djf")[39:59,170:215]
            if var == "tas_spatialmean_djf"
            if var == "pr_spatialmean_djf"
            
        
        # Y-AXIS
        for mod3, mod4 in zip(pattern_dir["nao_pattern_djf"].keys(),time_dir["nao_timeseries_djf"].keys()):   # the leading spatial pattern is standardized by its respective spatial standard deviation
            nao = variables_dir["nao_pattern_djf"][mod3][9:69,90:220] #[110:170,90:220] - remember changing this depending on the reverse_data
            nao_time = timeseries_dir["nao_timeseries_djf"][mod4]
            nao_std = np.nanstd(nao)
            new_nao_time = nao_time*nao_std
            std_nao_time = np.nanstd(new_nao_time)
            
            # add: std for each experiment minus std pi and then append to directory keyed to variable (when calling variable, all models will appear)

        # X-AXIS (changeble)
        for mod5, mod6 in zip(climatology_tas_dir["tas_spatialmean_djf"].keys, climatology_pr_dir["pr_spatialmean_djf"].keys):
            
            ### NORTH ###
            temperature_NEU = climatology_tas_dir["tas_spatialmean_djf"][mod5][19:39,170:215]
            mean_temp_NEU = np.nanmean(temperature_NEU)
            # add: each experiment - pi and append to directory keyed to variable (when calling variable, all models will appear)
            #y_axis.append()

            precipitation_NEU = climatology_pr_dir["pr_spatialmean_djf"][mod5][19:39,170:215]
            mean_pr_NEU = np.nanmean(precipitation_NEU)
            # add: each experiment - pi and append to directory keyed to variable (when calling variable, all models will appear)
            #y_axis.append()

            ### SOUTH ###
            temperature_MED = climatology_tas_dir["tas_spatialmean_djf"][mod5][39:59,170:215]
            mean_temp_MED = np.nanmean(temperature_MED)
            # add: each experiment - pi and append to directory keyed to variable (when calling variable, all models will appear)
            #y_axis.append()

            precipitation_MED = climatology_pr_dir["pr_spatialmean_djf"][mod5][39:59,170:215]
            # add: each experiment - pi and append to directory keyed to variable (when calling variable, all models will appear)
            mean_pr_MED = np.nanmean(precipitation_MED)
            
            #y_axis.append()
        
    return   """


# In[ ]:


for b in range(len(a)):
    plt.scatter(a[b],a[-1], color=colors[b], marker=markers[b], label="a")
    plt.xlabel("neu_tas")
    plt.ylabel("nao_amplitude")
    plt.legend(loc='best')


# In[51]:


hello = scatter_plots(experiments, "nao_pattern_djf")


# In[ ]:


##  SCATTERS ##

# Dictionary for values above - amplitude values for each model in each time period
# Take experiment amplitude minus pi amplotude for each model and put these difference values into a DIFF dictionary

# Take average TAS/PR djf over three boxes in the north atlantic (NEU,SEU,MED) and make a dictionary with average values
# for TAS and PR for each model in each time period.
# Take time period value minus pi value for each model and put these difference values into a DIFF dictionary

# MAKE SCATTER PLOT WITH amplitude change in x-axis and TAS/PR change in y-axis


# In[55]:


# amplitude for 20C

#import file
C20 = xr.open_dataset('test2/C20-Reanalysis.cvdp_data.1871-2012_corr2.nc', decode_times=False)

# select variable
C20_nao_patternn = C20.variables["nao_pattern_djf"][::-1,:]
C20_nao_time = C20.variables["nao_timeseries_djf"][:]

# regrid
grid_out = xe.util.grid_global(1,1)

# create grid
regridder_M1 = xe.Regridder(C20, grid_out, 'nearest_s2d', reuse_weights=False)

# use new grid
C20_nao_10 = regridder_M1(C20_nao_patternn.values)[:] # [9:69,90:220]

#plt.imshow(C20_nao_10)

# std of pattern

C20_nao_std = np.nanstd(C20_nao_10)
C20_nao_mean = np.nanmean(C20_nao_10)
New_C20_nao = (C20_nao_10-C20_nao_mean)/C20_nao_std
new_C20_time = C20_nao_time*C20_nao_std
std_C20_time = np.nanstd(new_C20_time)

#print(std_C20_time)
draw_models_one([C20_nao_10], [" "], "DELETE", lons, lats, 1, 1, -5,5)


# In[230]:


timeseries_dir["nao_timeseries_djf"]
#nao_time = timeseries_dir["nao_timeseries_djf"][mod]
#print(nao_time)


# # % of positive NAO seasons in time series

# In[8]:


for exp in experiments:
    print(exp)
    timeseries_dir = identify_ensemble_members(exp, "nao_timeseries_djf", True,False)
    
    for mod in timeseries_dir["nao_timeseries_djf"].keys():
            print("\t",mod)
            nao_time = timeseries_dir["nao_timeseries_djf"][mod]
            positive = nao_time>0 # how many years are positive?
            percentage = positive.sum()/len(nao_time)*100
            print("\t",percentage)


# In[ ]:




