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
import matplotlib as mpl


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
x = []
y = []
for i in range(1034):
    x.append(94.6071304991+i*(119.780064108-94.6071304991)/1034)
len(x)
for i in range(568):
    y.append(41.9862608605-i*(41.9862608605-31.2039974889)/568)
file_path = ("C:/Users/hp/Desktop/wang/")
#npp_nc = Dataset(file_path + "/ncproda_cgcm_yrb1km_rcp26.nc4")
npp_nc = Dataset(file_path + "/somtc_cgcm_yrb1km_rcp26.nc4")
#npp_data = npp_nc.variables['ncproda'][:].data
npp_data = npp_nc.variables['somtc'][:].data
npp_data[npp_data<10] = npp_data[0,0,0]
npp_data[npp_data==-9999] = np.nan
if os.path.exists(file_path+'SOMTC') != 1:
    os.mkdir(file_path+'SOMTC')
    
if os.path.exists(file_path+'SOMTC/2000_2015') != 1:
    os.mkdir(file_path+'SOMTC/2000_2015')
    
if os.path.exists(file_path+'SOMTC/2000_2050') != 1:
    os.mkdir(file_path+'SOMTC/2000_2050')
    
if os.path.exists(file_path+'SOMTC/2000_20501') != 1:
    os.mkdir(file_path+'SOMTC/2000_20501')    

if os.path.exists(file_path+'SOMTC/2035_2050') != 1:
    os.mkdir(file_path+'SOMTC/2035_2050')
    
fig,ax = plt.subplots(4,4,figsize=(5.5,4),sharex=False, sharey=False)
for i in range(4):
    for j in range(4):
        im = ax[i][j].contourf(x,y,npp_data[i*4+j+8],cmap = 'RdYlGn_r',levels=np.arange(0,7000.1,1000),extend='both')
        ax[i][j].axis('off')
        ax[i][j].text(102,41,str(1992+8+j+4*i),fontsize=6)
              
plt.subplots_adjust(wspace=-0.4,hspace=-0.11)
fig.colorbar(im, ax=[ax[0][0],ax[0][1],ax[0][2],ax[0][3],ax[1][0],ax[1][1],ax[1][2],ax[1][3],ax[2][0],ax[2][1],ax[2][2],ax[2][3],ax[3][0], ax[3][1],ax[3][2],ax[3][3]],fraction=0.047, pad = 0,orientation = 'horizontal',label='SOMTC (g C m$^-$$^2$)',extend='both')
plt.savefig(file_path + "SOMTC/2000_2015/2000_2015.jpg",dpi = 1000)

import matplotlib.colors as colors 

cmap = colors.ListedColormap(['green', 'limegreen','gold','orange','coral','orangered','red'])
#cmap = colors.ListedColormap([ 'CornflowerBlue','LightSkyBlue','gold','orange','coral','orangered','red'])
bounds=[0,5,10]
norm = colors.BoundaryNorm(bounds, cmap.N)
for i in range(0,59,5):
    plt.figure(figsize=(1034/100,568/100))
    im = plt.contourf(x,y,npp_data[i],cmap = cmap,levels=np.arange(0,7000.1,1000),extend='both')
    plt.text(102,40,str(1992+i),fontsize=16)
    plt.axis('off')
    plt.colorbar(im,fraction = 0.045, pad = - 0.12,orientation = 'horizontal',label='SOMTC (g C m$^-$$^2$)',extend='both')
#     plt.show()
    plt.savefig(file_path + "SOMTC/2000_2050/somtc_" + str(1992+i)+ ".jpg",dpi = 400,bbox_inches='tight')
    
from PIL import Image
import os
import os.path
import numpy as np
# import cv2
#指明被遍历的文件夹
rootdir = r'C:\Users\hp\Desktop\wang\SOMTC\2000_2050'
for parent, dirnames, filenames in os.walk(rootdir):#遍历每一张图片
    for filename in filenames:
        print('parent is :' + parent)
        print('filename is :' + filename)
        currentPath = os.path.join(parent, filename)
        print('the fulll name of the file is :' + currentPath)
   
        img = Image.open(currentPath)
        print (img.format, img.size, img.mode)
        #img.show()
        box1 = (800,100,2700,1977)#设置左、上、右、下的像素
        image1 = img.crop(box1) # 图像裁剪
        image1.save(r"C:\Users\hp\Desktop\wang\SOMTC\2000_2050n"+'\\'+filename) #存储裁剪得到的图像
        
import os
import imageio
 
def create_gif(image_list, gif_name):
 
    frames = []
    for image_name in image_list:
        if image_name.endswith('.jpg'):
            print(image_name)
            frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name, frames, 'GIF', duration = 0.5)
 
    return
 
def main():
 
    path = file_path + "SOMTC/2000_2050n/"
    a = [ path+img for img in  os.listdir(path)]
    image_list = []
    for name in a:
#         if (name[-5] == str(0))|(name[-5] == str(5)):
        image_list.append(name)
    gif_name = file_path + "SOMTC/2000_20501/somtc.gif"
    create_gif(image_list, gif_name)
    
if __name__ == "__main__":
    main()
    fig,ax = plt.subplots(4,4,figsize=(5.5,4),sharex=False, sharey=False)
    
    
for i in range(4):
    for j in range(4):
        im = ax[i][j].contourf(x,y,npp_data[i*4+j+8+35],cmap = 'RdYlGn_r',levels=np.arange(0,7000.1,1000),extend='both')
        ax[i][j].axis('off')
        ax[i][j].text(102,41,str(1992+8+j+4*i+35),fontsize=6)
              
plt.subplots_adjust(wspace=-0.4,hspace=-0.11)
fig.colorbar(im, ax=[ax[0][0],ax[0][1],ax[0][2],ax[0][3],ax[1][0],ax[1][1],ax[1][2],ax[1][3],ax[2][0],ax[2][1],ax[2][2],ax[2][3],ax[3][0], ax[3][1],ax[3][2],ax[3][3]],fraction=0.047, pad = 0,orientation = 'horizontal',label='SOMSC (g C m$^-$$^2$)',extend='both')
plt.savefig(file_path + "SOMTC/2035_2050/2035_2050.jpg",dpi = 1000)

npp_2050_2035 = np.full((npp_data.shape[1],npp_data.shape[2]),-1.0)
for i in range(43,59,1):
    npp_2050_2035 += npp_data[i]
npp_2050_2035 /= 16
plt.figure(figsize=(8,4))
im = plt.contourf(x,y,npp_2050_2035,cmap = 'RdYlGn_r',levels=np.arange(0,900.1,150),extend='both')
# plt.text(107,32,str(1992+8+i),fontsize=16)
plt.axis('off')
plt.colorbar(im,fraction = 0.05, pad = -0.01,orientation = 'horizontal',label='NPP (g C m$^-$$^2$ yr$^-$$^1$)',extend='both')
plt.savefig(file_path + "npp/2035_2050/even_35_50.jpg",dpi = 1000)
npp_2015_2000 = np.full((npp_data.shape[1],npp_data.shape[2]),-1.0)

for i in range(8,24,1):
    npp_2015_2000 += npp_data[i]
npp_2015_2000 /= 16
plt.figure(figsize=(8,4))
im = plt.contourf(x,y,npp_2015_2000,cmap = 'RdYlGn_r',levels=np.arange(0,900.1,150),extend='both')
plt.axis('off')
plt.colorbar(im,fraction = 0.05, pad = -0.01,orientation = 'horizontal',label='NPP (g C m$^-$$^2$ yr$^-$$^1$)',extend='both')
plt.savefig(file_path + "npp/2000_2015/even_00_15.jpg",dpi = 1000)
import matplotlib.colors as colors 
plt.figure(figsize=(8,4))

cmap = colors.ListedColormap(['green', 'limegreen','gold','orange','coral','orangered','red'])
bounds=[0,5,10]
norm = colors.BoundaryNorm(bounds, cmap.N)

# plt.hist2d(xvals, yvals, cmap=cmap)

npp_2050_2000 = plt.contourf(x,y,(npp_2050_2035-npp_2015_2000),levels=np.arange(-200,500.1,100),cmap = cmap,extend='both')
plt.colorbar(npp_2050_2000,fraction = 0.05, pad = -0.01,orientation = 'horizontal',label='NPP (g C m$^-$$^2$ yr$^-$$^1$)',extend='both')
plt.axis('off')
plt.savefig(file_path + "npp/npp_diff.jpg",dpi = 1000)