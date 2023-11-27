# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 22:31:05 2022
@author: Binb
"""

import os
import netCDF4 as nc
import numpy as np
from osgeo import gdal, osr, ogr
import glob

data = r'G:\GEMSupleft\reference_fips_chn1km_dv01.nc4'
nc_data_obj = nc.Dataset(data)
print(nc_data_obj.variables)


def nc2tif(data, Output_folder):
    tmp_data = nc.Dataset(data)  # 利用.Dataset()方法读取nc数据
    Lat_data = tmp_data.variables['lat'][:]
    Lon_data = tmp_data.variables['lon'][:]
    re_data = tmp_data.variables['SMsurf'][:]
    tmp_arr = np.asarray(tmp_data.variables['SMsurf'])
    # tmp_arr = np.asarray(tmp_data.variables['time'])

    # 影像的左上角&右下角坐标
    Lonmin, Latmax, Lonmax, Latmin = [Lon_data.min(), Lat_data.max(), Lon_data.max(), Lat_data.min()]
    print(Lonmin, Latmax, Lonmax, Latmin)

    # 分辨率计算
    Num_lat = len(Lat_data)  # 5146
    Num_lon = len(Lon_data)  # 7849
    Lat_res = (Latmax - Latmin) / (float(Num_lat) - 1)
    Lon_res = (Lonmax - Lonmin) / (float(Num_lon) - 1)
    # print(Num_lat, Num_lon)
    # print(Lat_res, Lon_res)
    for i in range(len(tmp_arr[:])):
        # i=0,1,2,3,4,5,6,7,8,9,...
        # 创建tif文件
        driver = gdal.GetDriverByName('GTiff')
        out_tif_name = Output_folder + '/' + data.split('/')[-1].split('.')[0] + '_' + str(i + 1) + '.tif'
        out_tif = driver.Create(out_tif_name, Num_lon, Num_lat, 1, gdal.GDT_Int16)
        # 设置影像的显示范围
        # Lat_re前需要添加负号
        geotransform = (Lonmin, Lon_res, 0.0, Latmax, 0.0, -Lat_res)
        out_tif.SetGeoTransform(geotransform)
        # 定义投影
        prj = osr.SpatialReference()
        prj.ImportFromEPSG(4326)  # WGS84
        out_tif.SetProjection(prj.ExportToWkt())
        # 数据导出
        out_tif.GetRasterBand(1).WriteArray(tmp_arr[i])  # 将数据写入内存，此时没有写入到硬盘
        out_tif.FlushCache()  # 将数据写入到硬盘
        out_tif = None  # 关闭tif文件


def main():
    # Input_folder = r"D:\\1234\\"
    Output_folder = r"E:/data/GLEAM/SM_SURF/"
    # 输出数据
    # data_list = glob.glob(os.path.join(Input_folder + '*.nc'))
    # print(data_list)
    # for i in range(len(data_list)):
    #     data = data_list[i]
    #     nc2tif(data, Output_folder)
    #     print('转tif成功')
    data = r'E:/data/GLEAM/nc/SMsurf_1980-2021_GLEAM_v3.6a_YR.nc'
    nc2tif(data, Output_folder)
    print('转tif成功')


main()