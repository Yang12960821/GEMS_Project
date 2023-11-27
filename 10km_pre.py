import calendar
from numpy.lib.type_check import nan_to_num
from osgeo import osr,ogr,gdal
from netCDF4 import Dataset
import numpy as np
import time
import os
import matplotlib.pyplot as plt
import scipy.stats as stats
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
    im_data=dataset.ReadAsArray(0,0,im_width,im_height).astype(np.float32)
    
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
    new_data = np.full((114,207),0.0,dtype=np.float32)
    for i in range(0,1136,10):
        for j in range(0,2067,10):
            new_data[int(i/10),int(j/10)] = data[i,j]
    return new_data[:]


def data_resample(array,resample_parameter):
    new_array = np.full((int(math.ceil(array.shape[0]*1.0/resample_parameter)),int(math.ceil(array.shape[1]*1.0/resample_parameter))),-999.0)
    for i in range(0,array.shape[0],resample_parameter):
        for j in range(0,array.shape[1],resample_parameter):
            new_array[int(i/resample_parameter),int(j/resample_parameter)] = array[i,j]
    return new_array

def extract_var_from_nc(file_path,var_name):
    file = Dataset(file_path)
    data = file.variables[var_name][:]
    file.close()
    return data
    
era5_pre_day = []
for year in range(1992,2016,1):
    #year_file = read_all_file("D:/zmx/LP_data/LP//pre/"+str(year))
    year_file = read_all_file("D:/zmx/LP_data/LP//pre/"+str(year))
    for name in year_file:
        era5_pre_day.append(read_rasterfile_setnan(name)[0])
    print(year)
    
era5_pre_day_np = np.array(era5_pre_day)
era5_pre_day_np.shape
era5_pre_day = []
for i in range(8766):
    era5_pre_day.append(down_resample(era5_pre_day_np[i]))
#     print(i)
era5_pre_day = np.array(era5_pre_day)
era5_pre_day.shape
root = Dataset("D:/zmx/LP_data/climate_ppt_yrb1km_1992_2015_dv01.nc4")
x = np.arange(8766,dtype = np.uint32)
mydataset.close()

mydataset = Dataset("D:/zmx/2.daya_nc/pre.nc4",'w')
mydataset.createDimension('x',era5_pre_day.shape[2])
mydataset.createDimension('y',era5_pre_day.shape[1])
mydataset.createDimension('time',era5_pre_day.shape[0])

lons = mydataset.createVariable('x',np.float32,('x',))
lons.long_name = "x distance on the projection plane from the UL corner"
lons.axis = "X"
lons.standard_name = "projection_x_coordinate"
lons.units = "meters"
lons[:] = root.variables['x'][0:2067:10]
lats = mydataset.createVariable('y',np.float32,('y',))
lats.long_name = "y distance on the projection plane from the UL corner"
lats.axis = "Y"
lats.standard_name = "projection_y_coordinate"
lats.units = "meter"
lats[:] = root.variables['y'][0:1136:10]
times=mydataset.createVariable('time','uint32',('time',))
times.standard_name='time'
times.long_time='time'
times.axis='T'
times.units = "days since 1992"
times.calendar = "gregorian"
times[:] = x[:]
Krasovsky_1940_Albers=mydataset.createVariable('albers_conical_equal_area','double')
Krasovsky_1940_Albers.grid_mapping_name = "albers_conical_equal_area"
Krasovsky_1940_Albers.standard_parallel = 25.0, 47.0
Krasovsky_1940_Albers.longitude_of_central_meridian = 105.0
Krasovsky_1940_Albers.latitude_of_projection_origin = 0
Krasovsky_1940_Albers.false_easting = 0.0
Krasovsky_1940_Albers.false_northing = 0.0
Krasovsky_1940_Albers.datum = "D_Krasovsky_1940"
Krasovsky_1940_Albers.GeoTransform = -852356.7064127066, 10000.0, 0.0, 4523979.996925818, 0.0, -10000.0
Krasovsky_1940_Albers._CoordinateTransformType = "Projection"
Krasovsky_1940_Albers._CoordinateSystems = "ProjectionCoordinateSystem"
Krasovsky_1940_Albers._CoordinateAxisTypes = "GeoX GeoY"

data_variable = mydataset.createVariable('prec','float',('time','y','x'),)
data_variable.long_time='precipitation'                                                
data_variable.units = "mm"
data_variable.grid_mapping = "albers_conical_equal_area"
data_variable[:] = era5_pre_day[:]
mydataset.close()