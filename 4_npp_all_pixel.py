from numpy.lib.type_check import nan_to_num
from osgeo import osr,ogr,gdal
from netCDF4 import Dataset
import numpy as np
import time
import os
# import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats
from sklearn.metrics import mean_squared_error,r2_score
import glob
import time
from functools import partial
from numba import jit
import csv
import time
import scipy
from scipy.stats import pearsonr
import math
from scipy.stats import gaussian_kde


def read_rasterfile(input_rasterfile):
    dataset=gdal.Open(input_rasterfile)
    im_width=dataset.RasterXSize
    im_height=dataset.RasterYSize
    im_bands=dataset.RasterCount
    im_geotrans=dataset.GetGeoTransform()
    im_proj=dataset.GetProjection()
    im_data=dataset.ReadAsArray(0,0,im_width,im_height)
    NoDataValue=dataset.GetRasterBand(1).GetNoDataValue()
    return [im_data,im_width,im_height,im_bands,im_geotrans,im_proj,NoDataValue]

def write_rasterfile(output_rasterfile,im_data,im_width,im_height,im_bands,im_geotrans,im_proj):
    driver=gdal.GetDriverByName("GTiff")
    datasetnew=driver.Create(output_rasterfile,im_width,im_height,im_bands,gdal.GDT_Float32)
    datasetnew.SetGeoTransform(im_geotrans)
    datasetnew.SetProjection(im_proj)
    datasetnew.GetRasterBand(1).SetNoDataValue(-9999.0)
    datasetnew.GetRasterBand(1).WriteArray(im_data)
    del datasetnew

def read_all_file(file_dir):
    final_file_list=[]
    final_file_list = glob.glob(os.path.join(file_dir,'*.tif'))
    final_file_list = sorted(final_file_list)
    return final_file_list

def read_rasterfile_setnan(input_rasterfile):
    dataset=gdal.Open(input_rasterfile)
    im_width=dataset.RasterXSize
    im_height=dataset.RasterYSize
    im_bands=dataset.RasterCount
    
    im_geotrans=dataset.GetGeoTransform()
    im_proj=dataset.GetProjection()
    im_data=dataset.ReadAsArray(0,0,im_width,im_height).astype(np.float)
    
    NoDataValue=dataset.GetRasterBand(1).GetNoDataValue()
    im_data[im_data==NoDataValue]=np.nan
    return [im_data,im_width,im_height,im_bands,im_geotrans,im_proj,NoDataValue]

def get_MSE_RMSE_R2_for_predicted(observations_array,predicted_array):
    mse=mean_squared_error(observations_array,predicted_array)
    rmse=np.sqrt(mse)
    r2=r2_score(observations_array,predicted_array)
    return mse,rmse,r2

def get_nparray_by_startNum_length_step(startNum,length,step):
    nparray=np.zeros(length)
    for i in range(length):
        nparray[i]=startNum+step*i
    return nparray

def deal_with_nc(ncfilepath,variable_name,variable_weidu):
    f = Dataset(ncfilepath)
    v = f.variables
    vndim=v[variable_name].ndim
    if vndim==3:
        rasterArray=v[variable_name][variable_weidu,:,:].data
    elif vndim==2:
        rasterArray=v[variable_name][:,:]
    elif vndim==1:
        rasterArray=v[variable_name][:]
    else:
        print ('data dimension wrong')
    return rasterArray

def pbias(estimate,observation):
    if observation==0:
        return 0
    else:
        return (estimate-observation)/observation

@jit
def process(data,new_data):
    n = 3
    for i in range(0,58):
        print(i)
        for j in range(0,(data.shape)[1]):
            for k in range(0,(data.shape)[2]):
                for ii in range(0,n):
                    for jj in range(0,n):
                        if j < (data.shape)[1] and k < (data.shape)[2] and np.isnan(data[i][j][k]) == False:
                            new_data [i][j*n+ii][k*n+jj] = data [i][j][k]
        print(i)


    print(new_data)
    return new_data[:]


@jit
def down_resample(data):
    new_data = np.full((568,1034),0.0)
    for i in range(0,1136,2):
        for j in range(0,2067,2):
            new_data[int(i/2),int(j/2)] = data[i,j]
    return new_data[:]


def data_resample(array,resample_parameter):
    new_array = np.full((int(math.ceil(array.shape[0]*1.0/resample_parameter)),int(math.ceil(array.shape[1]*1.0/resample_parameter))),-999.0)
    for i in range(0,array.shape[0],resample_parameter):
        for j in range(0,array.shape[1],resample_parameter):
            new_array[int(i/resample_parameter),int(j/resample_parameter)] = array[i,j]
    return new_array
    
gems_nc.close()
gems_nc = Dataset("C:/Users/hp/Downloads/wang20220316173017/ncproda_cgcm_yrb1km_rcp26.nc4")
crop_type = Dataset("F:/YRB_NC_files/management/mgmt_crop_type_yrb1km_1992_2015_dv02.nc4")
all_fips_filpath=u'G:/gems/GEMS_NPP_MODIS_VS/FIPS/yrb_FIPS_raster_1km.tif'
all_modis_npp_files = read_all_file('G:/gems/GEMS_NPP_MODIS_VS/YRB_MOD17A3_NPP_1km_kgC_m2')


fips=read_rasterfile_setnan(all_fips_filpath)[0]
fips = down_resample(fips)
all_fips=np.unique(fips[~np.isnan(fips)])
gems_all_data = gems_nc.variables["ncproda"][:].data
crop_type_data = crop_type.variables["crop"][:].data

gems_all_data[gems_all_data <= 10]=np.nan
year = [0,1,2,3,4,5,6,7]
all_del_gems_data = []
all_del_modis_data = []
all_ever_gems_data = []
all_ever_modis_data = []
all_grassland_gems_data = []
all_grassland_modis_data = []
all_shrubland_gems_data = []
all_shrubland_modis_data = []

    
for i in year:
    gems_year_data = gems_all_data[i+8]

    modis_year_data=read_rasterfile_setnan(all_modis_npp_files[i])[0]*1000.0
    modis_year_data = down_resample(modis_year_data)

    crop_type_year_data = crop_type_data[8+i]
    crop_type_year_data = down_resample(crop_type_year_data)
    
    for fip in all_fips:
        all_del_gems_data.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==8)]))
        all_del_modis_data.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==8) ]))
        
        all_ever_gems_data.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==9) ]))
        all_ever_modis_data.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==9)]))
        
        all_grassland_gems_data.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==11) ]))
        all_grassland_modis_data.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==11) ]))
        
        all_shrubland_gems_data.append( np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==12)]))
        all_shrubland_modis_data.append(np.nanmean(modis_year_data[(fips==fip)  & (crop_type_year_data==12)]))
    print(str(i+2000) + " down")
    print(len(all_del_gems_data),len(all_del_modis_data),len(all_ever_gems_data),len(all_ever_modis_data),len(all_grassland_gems_data),len(all_grassland_modis_data),len(all_shrubland_gems_data),len(all_shrubland_modis_data))


all_del_gems_data_1 = []
all_del_modis_data_1 = []
all_ever_gems_data_1 = []
all_ever_modis_data_1 = []
all_grassland_gems_data_1 = []
all_grassland_modis_data_1 = []
all_shrubland_gems_data_1 = []
all_shrubland_modis_data_1 = []

    
for i in range(8,16,1):
    gems_year_data = gems_all_data[i+8]

    modis_year_data=read_rasterfile_setnan(all_modis_npp_files[i])[0]*1000.0
    modis_year_data = down_resample(modis_year_data)

    crop_type_year_data = crop_type_data[8+i]
    crop_type_year_data = down_resample(crop_type_year_data)
    
    for fip in all_fips:
        all_del_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==8)]))
        all_del_modis_data_1.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==8) ]))
        
        all_ever_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==9) ]))
        all_ever_modis_data_1.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==9)]))
        
        all_grassland_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==11) ]))
        all_grassland_modis_data_1.append(np.nanmean(modis_year_data[(fips==fip) & (crop_type_year_data==11) ]))
        
        all_shrubland_gems_data_1.append( np.nanmean(gems_year_data[(fips==fip) & (crop_type_year_data==12)]))
        all_shrubland_modis_data_1.append(np.nanmean(modis_year_data[(fips==fip)  & (crop_type_year_data==12)]))
    print(str(i+2000) + " down")
    print(len(all_del_gems_data_1),len(all_del_modis_data_1),len(all_ever_gems_data_1),len(all_ever_modis_data_1),len(all_grassland_gems_data_1),len(all_grassland_modis_data_1),len(all_shrubland_gems_data_1),len(all_shrubland_modis_data_1))


modis_cal = ([all_del_modis_data,all_ever_modis_data,all_grassland_modis_data,all_shrubland_modis_data])
gems_cal  = ([all_del_gems_data ,all_ever_gems_data ,all_grassland_gems_data ,all_shrubland_gems_data ])
modis_val = ([all_del_modis_data_1,all_ever_modis_data_1,all_grassland_modis_data_1,all_shrubland_modis_data_1])
gems_val  = ([all_del_gems_data_1 ,all_ever_gems_data_1 ,all_grassland_gems_data_1 ,all_shrubland_gems_data_1])

fig,ax = plt.subplots(4,2,figsize=(6,10),sharex=False, sharey=False)
# fig, ax = plt.subplots(figsize=(30, 30))
fontdict={'fontsize': 8, 'fontweight': rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'family':'arial'}
plt.rcParams['font.size']=8
name = ['Deciduous','Evergreen','Grassland','Shrubland']
for i in range(4):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_cal[i])):
        if modis_cal[i][j]!=0 and gems_cal[i][j]!=0 and ~np.isnan(modis_cal[i][j]) and ~np.isnan(gems_cal[i][j]):
            soybeans_gems_final_data.append(gems_cal[i][j])
            soybeans_modis_final_data.append(modis_cal[i][j])  
            
    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)
    
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='gray',linestyle='--',zorder=0)
    ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
    ax[i][0].set(xlim=(0,1000))
    ax[i][0].set(ylim=(0,1000))
    ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
    ax[i][0].text(40,920,"RMSE = " + str(rmse), fontsize=8)
    ax[i][0].text(40,830,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][0].text(40,740,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    ax[i][0].tick_params(axis='both',which='both',direction='in')
    
ax[0][0].text(600,100,name[0],fontsize=10)
ax[1][0].text(600,100,name[1],fontsize=10)
ax[2][0].text(600,100,name[2],fontsize=10)
ax[3][0].text(600,100,name[3],fontsize=10)

# ax[0][0].spines['right'].set_visible(False)
# ax[0][0].spines['top'].set_visible(False)
    
    
for i in range(4):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_val[i])):
        if modis_val[i][j]!=0 and gems_val[i][j]!=0 and ~np.isnan(modis_val[i][j]) and ~np.isnan(gems_val[i][j]):
            soybeans_gems_final_data.append(gems_val[i][j])
            soybeans_modis_final_data.append(modis_val[i][j])  
            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)

    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='gray',linestyle='--',zorder=0)
    ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
#     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=1,marker='+',c='red')
    ax[i][1].set(xlim=(0,1000))
    ax[i][1].set(ylim=(0,1000))
    ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
    ax[i][1].text(40,920,"RMSE = " + str(rmse), fontsize=8)
    ax[i][1].text(40,830,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][1].text(40,740,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    ax[i][1].tick_params(axis='both',which='both',direction='in')
#     ax[i][1].text(615,20,name[i],fontsize=8,fontweight ='bold')
    
ax[0][1].text(600,100,name[0],fontsize=10)
ax[1][1].text(600,100,name[1],fontsize=10)
ax[2][1].text(600,100,name[2],fontsize=10)
ax[3][1].text(600,100,name[3],fontsize=10)

ax[0][0].set_title("Calibration",fontsize=10,fontweight ='bold',pad=10)
ax[0][0].xaxis.set_label_position('top') 
ax[0][1].set_title("Validation",fontsize=10,fontweight ='bold',pad=10)
ax[0][1].xaxis.set_label_position('top') 

ax[3][0].set(xlabel='MODIS Npp\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[3][1].set(xlabel='MODIS Npp\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[0][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[1][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[3][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')

    
plt.subplots_adjust(wspace=0.22,hspace=0.22)    
plt.savefig("C:/Users/hp/Desktop/npp1.png",dpi=1000)

fig,ax = plt.subplots(4,2,figsize=(6,10),sharex=False, sharey=False)
# fig, ax = plt.subplots(figsize=(30, 30))
fontdict={'fontsize': 8, 'fontweight': rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'family':'arial'}
plt.rcParams['font.size']=8
name = ['Deciduous','Evergreen','Grassland','Shrubland']
for i in range(4):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_cal[i])):
        if modis_cal[i][j]!=0 and gems_cal[i][j]!=0 and ~np.isnan(modis_cal[i][j]) and ~np.isnan(gems_cal[i][j]):
            soybeans_gems_final_data.append(gems_cal[i][j])
            soybeans_modis_final_data.append(modis_cal[i][j])  
            
    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)
    
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='gray',linestyle='--',zorder=0)
    ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
    ax[i][0].set(xlim=(0,1000))
    ax[i][0].set(ylim=(0,1000))
    ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
    ax[i][0].text(40,920,"RMSE = " + str(rmse), fontsize=8)
    ax[i][0].text(40,830,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][0].text(40,740,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    ax[i][0].tick_params(axis='both',which='both',direction='in')
    
ax[0][0].text(600,100,name[0],fontsize=10)
ax[1][0].text(600,100,name[1],fontsize=10)
ax[2][0].text(600,100,name[2],fontsize=10)
ax[3][0].text(600,100,name[3],fontsize=10)

# ax[0][0].spines['right'].set_visible(False)
# ax[0][0].spines['top'].set_visible(False)
    
    
for i in range(4):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_val[i])):
        if modis_val[i][j]!=0 and gems_val[i][j]!=0 and ~np.isnan(modis_val[i][j]) and ~np.isnan(gems_val[i][j]):
            soybeans_gems_final_data.append(gems_val[i][j])
            soybeans_modis_final_data.append(modis_val[i][j])  
            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)

    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='gray',linestyle='--',zorder=0)
    ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
#     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=1,marker='+',c='red')
    ax[i][1].set(xlim=(0,1000))
    ax[i][1].set(ylim=(0,1000))
    ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
    ax[i][1].text(40,920,"RMSE = " + str(rmse), fontsize=8)
    ax[i][1].text(40,830,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][1].text(40,740,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    ax[i][1].tick_params(axis='both',which='both',direction='in')
#     ax[i][1].text(615,20,name[i],fontsize=8,fontweight ='bold')
    
ax[0][1].text(600,100,name[0],fontsize=10)
ax[1][1].text(600,100,name[1],fontsize=10)
ax[2][1].text(600,100,name[2],fontsize=10)
ax[3][1].text(600,100,name[3],fontsize=10)

ax[0][0].set_title("Calibration",fontsize=10,fontweight ='bold',pad=10)
ax[0][0].xaxis.set_label_position('top') 
ax[0][1].set_title("Validation",fontsize=10,fontweight ='bold',pad=10)
ax[0][1].xaxis.set_label_position('top') 

ax[3][0].set(xlabel='MODIS Npp\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[3][1].set(xlabel='MODIS Npp\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[0][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[1][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[3][0].set(ylabel='GEMS NPP\n(g C m$^-$$^2$ yr$^-$$^1$)')

    
plt.subplots_adjust(wspace=0.22,hspace=0.22)    
plt.savefig("C:/Users/hp/Desktop/npp1.png",dpi=1000)

root = Dataset("G:/3.LP_crop_type/mgmt_crops_chn1km_1992_2015_tab_dv01.nc4")
crop_harvest = root.variables['crop_harvest'][:]

gems_nc = Dataset("C:/Users/hp/Desktop/gems/wang20220113144843/cgrain_cgcm_yrb1km_rcp26.nc4")
crop_type = Dataset("F:/YRB_NC_files/management/mgmt_crop_type_yrb1km_1992_2015_dv02.nc4")
crop_type_data = crop_type.variables["crop"][:].data

gems_all_data = gems_nc.variables['cgrain'][:].data
gems_all_data[gems_all_data<10]=np.nan

f = Dataset("G:/3.LP_crop_type/mgmt_crops_chn1km_1992_2015_tab_dv01.nc4")
v = f.variables
#获取每一个fips
china_fips=list(np.unique(v['crop_harvest'][:,0,0]['county_id']))
#构建字典，存储每个fips每一年的7种耕作方式的产量总和
fip_a_cgrain={}
for i in np.arange(len(china_fips)):
    #只要2000年到2015年的数据
    each_year_cgrain=np.full([16,29],-9999.0)
    for y in range(16):
        #取第一种耕作方式'Total Crop'
        each_year_cgrain[y]=v['crop_harvest'][i][y+8]['yield'][0]
    #存储到该fip字典
    fip_a_cgrain[china_fips[i]]=each_year_cgrain
    
all_corn_gems_data = []
all_corn_modis_data = []
all_soybeans_gems_data = []
all_soybeans_modis_data = []
all_wheat_gems_data = []
all_wheat_modis_data = []

for i in range(8):
    gems_year_data = gems_all_data[i+8]
    crop_type_year_data = crop_type_data[8+i]
    crop_type_year_data = down_resample(crop_type_year_data)
    for fip in all_fips:
        if fip in fip_a_cgrain:
            all_corn_gems_data.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data>=103)&(crop_type_year_data<=104)]))
            all_corn_modis_data.append(fip_a_cgrain[fip][i][3])

            all_soybeans_gems_data.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data==119)]))
            all_soybeans_modis_data.append(fip_a_cgrain[fip][i][19])

            all_wheat_gems_data.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data==128)]))
            all_wheat_modis_data.append(fip_a_cgrain[fip][i][28])
    print(str(i+2000) + " down")
    print(len(all_corn_gems_data),len(all_corn_modis_data),len(all_soybeans_gems_data),len(all_soybeans_modis_data),len(all_wheat_gems_data),len(all_wheat_modis_data))

all_corn_gems_data_1 = []
all_corn_modis_data_1 = []
all_soybeans_gems_data_1 = []
all_soybeans_modis_data_1 = []
all_wheat_gems_data_1 = []
all_wheat_modis_data_1 = []

for i in range(8,16,1):
    gems_year_data = gems_all_data[i+8]
    crop_type_year_data = crop_type_data[8+i]
    crop_type_year_data = down_resample(crop_type_year_data)
    for fip in all_fips:
        if fip in fip_a_cgrain:
            all_corn_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data>=103)&(crop_type_year_data<=104)]))
            all_corn_modis_data_1.append(fip_a_cgrain[fip][i][3])

            all_soybeans_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data==119)]))
            all_soybeans_modis_data_1.append(fip_a_cgrain[fip][i][19])

            all_wheat_gems_data_1.append(np.nanmean(gems_year_data[(fips==fip)&(crop_type_year_data==128)]))
            all_wheat_modis_data_1.append(fip_a_cgrain[fip][i][28])
    print(str(i+2000) + " down")
    print(len(all_corn_gems_data_1),len(all_corn_modis_data_1),len(all_soybeans_gems_data_1),len(all_soybeans_modis_data_1),len(all_wheat_gems_data_1),len(all_wheat_modis_data_1))

modis_cal = ([all_corn_modis_data,all_soybeans_modis_data,all_wheat_modis_data])
gems_cal   = ([all_corn_gems_data ,all_soybeans_gems_data ,all_wheat_gems_data])
modis_val = ([all_corn_modis_data_1,all_soybeans_modis_data_1,all_wheat_modis_data_1])
gems_val   = ([all_corn_gems_data_1 ,all_soybeans_gems_data_1 ,all_wheat_gems_data_1])

for i in range(len(modis_cal[2])):
    if modis_cal[2][i]>250:
        modis_cal[2][i]=0
for i in range(len(modis_val[2])):
    if modis_val[2][i]>250:
        modis_val[2][i]=0
        
for i in range(len(gems_cal[2])):
    if gems_cal[2][i]>250:
        gems_cal[2][i]=0
for i in range(len(gems_val[2])):
    if gems_val[2][i]>250:
        gems_val[2][i]=0

fig,ax = plt.subplots(3,2,figsize=(6,8),sharex=False, sharey=False)
# fig, ax = plt.subplots(figsize=(30, 30))
fontdict={'fontsize': 8, 'fontweight': rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'family':'arial'}
plt.rcParams['font.size']=8
name = ['corn','soybeans','wheat']

for i in range(3):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []

    
    for j in range(len(modis_val[i])):
        if modis_cal[i][j]!=0 and gems_cal[i][j]!=0 and ~np.isnan(modis_cal[i][j]) and ~np.isnan(gems_cal[i][j]):
            soybeans_gems_final_data.append(gems_cal[i][j])
            soybeans_modis_final_data.append(modis_cal[i][j])
    
    

    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)
    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)
    

    
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='black',linestyle='--',zorder=0)
#     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
    ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
    ax[i][0].set(xlim=(0,500))
    ax[i][0].set(ylim=(0,500))
    ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='--')
    ax[i][0].text(15,460,"RMSE = " + str(rmse), fontsize=8)
    ax[i][0].text(15,420,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][0].text(15,380,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)

ax[0][0].text(320,60,name[0],fontsize=10)
ax[1][0].text(320,60,name[1],fontsize=10)
ax[2][0].text(320,60,name[2],fontsize=10)
    
    
for i in range(3):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_val[i])):
        if modis_val[i][j]!=0 and gems_val[i][j]!=0 and ~np.isnan(modis_val[i][j]) and ~np.isnan(gems_val[i][j]):
            soybeans_gems_final_data.append(gems_val[i][j])
            soybeans_modis_final_data.append(modis_val[i][j])  

            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)

    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    
    ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='black',linestyle='--',zorder=0)
#     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='mediumturquoise',linewidth=0.7)
#     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=7,marker='o',c='salmon',linewidth=0.5)
    ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
    ax[i][1].set(xlim=(0,500))
    ax[i][1].set(ylim=(0,500))
    ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='--')
    ax[i][1].text(15,460,"RMSE = " + str(rmse), fontsize=8)
    ax[i][1].text(15,420,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
    ax[i][1].text(15,380,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    
ax[0][1].text(320,60,name[0],fontsize=10)
ax[1][1].text(320,60,name[1],fontsize=10)
ax[2][1].text(320,60,name[2],fontsize=10)
    

ax[0][0].set_title("Calibration",fontsize=10,fontweight ='bold',pad=10)
ax[0][0].xaxis.set_label_position('top') 

ax[0][1].set_title("Validation",fontsize=10,fontweight ='bold',pad=10)
ax[0][1].xaxis.set_label_position('top') 


ax[2][0].set(xlabel='Obs Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][1].set(xlabel='Obs Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')


ax[0][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[1][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')




plt.subplots_adjust(wspace=0.22,hspace=0.22)   
plt.savefig("C:/Users/hp/Desktop/cgrain.png",dpi=1000)

fig,ax = plt.subplots(3,2,figsize=(6,8),sharex=False, sharey=False)
# fig, ax = plt.subplots(figsize=(30, 30))
fontdict={'fontsize': 8, 'fontweight': rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'family':'arial'}
plt.rcParams['font.size']=8
name = ['corn','soybeans','wheat']

for i in range(3):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []

    
    for j in range(len(modis_val[i])):
        if modis_cal[i][j]!=0 and gems_cal[i][j]!=0 and ~np.isnan(modis_cal[i][j]) and ~np.isnan(gems_cal[i][j]):
            soybeans_gems_final_data.append(gems_cal[i][j])
            soybeans_modis_final_data.append(modis_cal[i][j])
    
    

    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)
    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)
    

    
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    if i == 0:
        ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
    #     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
        ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][0].set(xlim=(0,500))
        ax[i][0].set(ylim=(0,500))
        ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][0].text(20,450,"RMSE = " + str(rmse), fontsize=8)
        ax[i][0].text(20,400,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][0].text(20,350,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
        
    if i == 1:
        ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
    #     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
        ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][0].set(xlim=(0,300))
        ax[i][0].set(ylim=(0,300))
        ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][0].text(15,270,"RMSE = " + str(rmse), fontsize=8)
        ax[i][0].text(15,240,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][0].text(15,210,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8) 
        
    if i == 2:
        ax[i][0].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
    #     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
        ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][0].set(xlim=(0,300))
        ax[i][0].set(ylim=(0,300))
        ax[i][0].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][0].text(15,270,"RMSE = " + str(rmse), fontsize=8)
        ax[i][0].text(15,240,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][0].text(15,210,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    ax[i][0].tick_params(axis='both',which='both',direction='in')

ax[0][0].text(320,60,name[0],fontsize=10)
ax[1][0].text(180,35,name[1],fontsize=10)
ax[2][0].text(180,35,name[2],fontsize=10)
    
    
for i in range(3):
    soybeans_gems_final_data = []
    soybeans_modis_final_data = []
    
    for j in range(len(modis_val[i])):
        if modis_val[i][j]!=0 and gems_val[i][j]!=0 and ~np.isnan(modis_val[i][j]) and ~np.isnan(gems_val[i][j]):
            soybeans_gems_final_data.append(gems_val[i][j])
            soybeans_modis_final_data.append(modis_val[i][j])  

            
    soybeans_gems_final_data = np.array(soybeans_gems_final_data)
    soybeans_modis_final_data = np.array(soybeans_modis_final_data)

    print(soybeans_gems_final_data.shape)
    print(soybeans_modis_final_data.shape)

    xy = np.vstack([soybeans_gems_final_data,soybeans_modis_final_data])
    z = gaussian_kde(xy)(xy)
    
    r2=str(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[2])[:4]
    rmse=round(get_MSE_RMSE_R2_for_predicted(soybeans_gems_final_data,soybeans_modis_final_data)[1],2)
    pear = pearsonr(soybeans_gems_final_data,soybeans_modis_final_data)[0]
    ra = round(pear*pear,2)

    sum_gems = 0
    sum_modis = 0
    for jj in range(len(soybeans_gems_final_data)):
        sum_gems = sum_gems +  soybeans_gems_final_data[jj]
        sum_modis = sum_modis + soybeans_modis_final_data[jj]
    bias = pbias(sum_gems/len(soybeans_gems_final_data),sum_modis/len(soybeans_gems_final_data))

    A=np.vstack([soybeans_modis_final_data,np.ones(len(soybeans_modis_final_data))]).T
    m,c=np.linalg.lstsq(A,soybeans_gems_final_data,rcond=None)[0]


    if i==0:
        ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
        #     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='mediumturquoise',linewidth=0.7)
        #     ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=7,marker='o',c='salmon',linewidth=0.5)
        ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][1].set(xlim=(0,500))
        ax[i][1].set(ylim=(0,500))
        ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][1].text(20,450,"RMSE = " + str(rmse), fontsize=8)
        ax[i][1].text(20,400,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][1].text(20,350,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    if i == 1:
        ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
    #     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
        ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][1].set(xlim=(0,300))
        ax[i][1].set(ylim=(0,300))
        ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][1].text(15,270,"RMSE = " + str(rmse), fontsize=8)
        ax[i][1].text(15,240,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][1].text(15,210,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8) 
        
    if i == 2:
        ax[i][1].plot([0,2000],[0,2000],linewidth=1,color='grey',linestyle='--',zorder=0)
    #     ax[i][0].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=13,marker='o',c='white',edgecolors='dodgerblue',linewidth=0.7)
        ax[i][1].scatter(soybeans_modis_final_data,soybeans_gems_final_data,s=2,marker='o',c=z,cmap=plt.cm.jet,linewidth=0.5)
        ax[i][1].set(xlim=(0,300))
        ax[i][1].set(ylim=(0,300))
        ax[i][1].plot(np.arange(0,2000,1),np.arange(0,2000,1)*m+c,linewidth=1,color='red',linestyle='-')
        ax[i][1].text(15,270,"RMSE = " + str(rmse), fontsize=8)
        ax[i][1].text(15,240,'R'+'\u00b2'+ " = " + str(ra), fontsize=8)
        ax[i][1].text(15,210,"PB = " + str(round(float(bias)*100,2))+'% ', fontsize=8)
    
    
    
    ax[i][1].tick_params(axis='both',which='both',direction='in')
    
ax[0][1].text(320,60,name[0],fontsize=10)
ax[1][1].text(180,35,name[1],fontsize=10)
ax[2][1].text(180,35,name[2],fontsize=10)
    

ax[0][0].set_title("Calibration",fontsize=10,fontweight ='bold',pad=10)
ax[0][0].xaxis.set_label_position('top') 

ax[0][1].set_title("Validation",fontsize=10,fontweight ='bold',pad=10)
ax[0][1].xaxis.set_label_position('top') 


ax[2][0].set(xlabel='Obs Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][1].set(xlabel='Obs Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')


ax[0][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[1][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')
ax[2][0].set(ylabel='GEMS Grain\n(g C m$^-$$^2$ yr$^-$$^1$)')




plt.subplots_adjust(wspace=0.22,hspace=0.22)   
plt.savefig("C:/Users/hp/Desktop/cgrain1.png",dpi=1000)