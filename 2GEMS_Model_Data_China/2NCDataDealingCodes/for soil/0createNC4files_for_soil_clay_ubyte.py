# -*- coding: utf-8 -*-
"""
create NC4files
"""

import os
import time
import numpy as np
from osgeo import osr,ogr,gdal
from netCDF4 import Dataset


def deal_with_nc(ncfilepath,variable_name,variable_weidu):
	#print ncfilepath
	ncfile=ncfilepath.encode('GB2312')  #避免出现中文不支持的问题
	f = Dataset(ncfile)
	v = f.variables
	vndim=v[variable_name].ndim
	if vndim==3:
		rasterArray=v[variable_name][variable_weidu,:,:]	
	elif vndim==2:
		rasterArray=v[variable_name][:,:]
	elif vndim==1:
		rasterArray=v[variable_name][:]
	else:
		print 'data dimension wrong'
	return rasterArray


#递归读取文件夹下所有的特定格式的文件
def read_all_files_inside_filepath(filepath,filetype,final_file_list):
	pathDir=os.listdir(filepath)
	for each in pathDir:
		newDir=os.path.join(filepath,each)
		if os.path.isfile(newDir):
			if os.path.splitext(newDir)[1]==filetype:
				final_file_list.append(newDir.replace("\\","/"))				
		else:
			read_all_files_inside_filepath(newDir,filetype,final_file_list)


def write_rasterfile(output_rasterfile,im_data,im_width,im_height,im_bands,im_geotrans,im_proj):	
	driver=gdal.GetDriverByName("GTiff")
	datasetnew=driver.Create(output_rasterfile,im_width,im_height,im_bands,gdal.GDT_Float32)
	datasetnew.SetGeoTransform(im_geotrans)
	datasetnew.SetProjection(im_proj)
	datasetnew.GetRasterBand(1).SetNoDataValue(65535.0)
	datasetnew.GetRasterBand(1).WriteArray(im_data)
	del datasetnew



def convert_nclistfiles_to_raster(nc_files_list,nc_variables,output_rasterfile_path,im_width,im_height,im_bands,im_geotrans,im_proj):
	totalstart=time.time()
	print u"转换开始..."
	for each in nc_files_list:		
		#print each.replace("\\","/")
		print u"第"+str(nc_files_list.index(each)+1)+u"开始..."
		starttime=time.time()
		output_rasterfile=output_rasterfile_path+"/"+"stand_age.tif"
		im_data=deal_with_nc(each.replace("\\","/"),nc_variables,7)
		write_rasterfile(output_rasterfile,im_data,im_width,im_height,im_bands,im_geotrans,im_proj)		
		endtime=time.time()
		print str(round((endtime-starttime)/60.0,2))+"  min"+"\n"
	totalend=time.time()
	print "All done, total time:"+str(round((totalend-totalstart)/3600.0,2))+"  h"+"\n"


def convert_nclistfiles_to_raster(nc_files_list,nc_variables,output_rasterfile_path,im_width,im_height,im_bands,im_geotrans,im_proj,variable_weidu,start_time):
	for each in nc_files_list:
		output_rasterfile=output_rasterfile_path+"/"+nc_variables+"_"+str(start_time+int(variable_weidu))+".tif"
		im_data=deal_with_nc(each.replace("\\","/"),nc_variables,variable_weidu)
		write_rasterfile(output_rasterfile,im_data,im_width,im_height,im_bands,im_geotrans,im_proj)


def read_rasterfile(input_rasterfile):
	dataset=gdal.Open(input_rasterfile)
	#读取栅格数据基本属性
	im_width=dataset.RasterXSize
	im_height=dataset.RasterYSize
	im_bands=dataset.RasterCount
	#读取栅格数据坐标信息和投影信息
	im_geotrans=dataset.GetGeoTransform()
	im_proj=dataset.GetProjection()
	#读取栅格数据的具体值
	im_data=dataset.ReadAsArray(0,0,im_width,im_height)	
	#获取空值
	NoDataValue=dataset.GetRasterBand(1).GetNoDataValue()
	return [im_data,im_width,im_height,im_bands,im_geotrans,im_proj,NoDataValue]

#根据起始值，数据量和步长构建array
def get_nparray_by_startNum_length_step(startNum,length,step):
	nparray=np.zeros(length)
	for i in range(length):
		nparray[i]=startNum+step*i
	return nparray



#对这样的数据进行排序
#XXXX/precip.1979_1.img
def list_sort_by_data(list_data):
	return int(list_data.split('/')[-1].split('.')[1][:4])*100+int(list_data.split('/')[-1].split('.')[1].split('_')[-1])

#进度条展示
def print_current_process(process):
	currentbarlength=int(process*20)
	alllength=""
	for i in range(currentbarlength):
		alllength=alllength+"*"
	print alllength+"..."+str(process*100.0)+"%"

#input_rasterfiles_type='.img'
#data_variable_name='CO2_ppm'
#all_vars_units=[times,lats,lons,data_variable]
#times_units='months since 1970'
def create_NC_by_raster(input_rasterfiles_path,input_rasterfiles_type,output_NC_file,data_variable_name,data_variable_longname,all_vars_units,Scale_factor,begin_end_date,filetitle):
	final_file_list=[]
	read_all_files_inside_filepath(input_rasterfiles_path,input_rasterfiles_type,final_file_list)
	#降雨、温度需要排序
	#final_file_list.sort(key=list_sort_by_data)
	if len(final_file_list)>0:
		rasterinfo=read_rasterfile(final_file_list[0])
		im_width,im_height,im_bands,im_geotrans,im_proj,NoDataValue=rasterinfo[1],rasterinfo[2],rasterinfo[3],rasterinfo[4],rasterinfo[5],rasterinfo[6]
	else:
		print 'no \''+input_rasterfiles_type+'\' type of data in the path: '+input_rasterfiles_path
		exit()
	#Fill_value=NoDataValue
	mydataset=Dataset(output_NC_file,'w',format='NETCDF4')
	#构建维度
	mydataset.createDimension('soil_layers_top',len(final_file_list))
	mydataset.createDimension('soil_layers_bottom',len(final_file_list))
	mydataset.createDimension('y',im_height)
	mydataset.createDimension('x',im_width)
	#char类型投影参数变量
	#mydataset.createDimension('Krasovsky_1940_Albers',1)
	#构建变量
	soil_layers_top=mydataset.createVariable('soil_layers_top','u2',('soil_layers_top',))
	soil_layers_bottom=mydataset.createVariable('soil_layers_bottom','u2',('soil_layers_bottom',))
	lats=mydataset.createVariable('y',np.float32,('y',))
	lons=mydataset.createVariable('x',np.float32,('x',))
	#投影变量
	Krasovsky_1940_Albers=mydataset.createVariable('albers_conical_equal_area','S1')
	#再存储数据的时候，需要先存放纬度，再存放精度
	if len(soil_layers_top)>1:
		data_variable=mydataset.createVariable(data_variable_name,'u1',('soil_layers_top','y','x'),fill_value=255)
	else:
		data_variable=mydataset.createVariable(data_variable_name,'u1',('y','x'),fill_value=255)
	#对变量的单位进行设定
	soil_layers_top.units=all_vars_units[0]
	soil_layers_bottom.units=all_vars_units[0]
	lats.units=all_vars_units[1]
	lons.units=all_vars_units[2]
	data_variable.units=all_vars_units[3]
	#变量的相关参数
	#soil_layers_top.standard_name='soil_layers_top'
	soil_layers_top.long_time="top depth of soil layers"
	soil_layers_top.axis='soil_layers_top'
	soil_layers_bottom.long_time="lower depth of soil layers"
	soil_layers_bottom.axis='soil_layers_bottom'
	data_variable.long_name=data_variable_longname
	lats.long_name="y distance on the projection plane from the UL corner"
	lats.standard_name="projection_y_coordinate"
	lats.axis = "Y"
	lons.long_name="x distance on the projection plane from the UL corner"
	lons.standard_name="projection_x_coordinate"
	lons.axis = "X"
	Krasovsky_1940_Albers.grid_mapping_name = "albers_conical_equal_area"
	Krasovsky_1940_Albers.standard_parallel = (25.0,47.0)
	Krasovsky_1940_Albers.longitude_of_central_meridian=105.0
	#是参数latitude_of_projection_origin的问题，设置为0就正确了，不是经纬度范围的均值
	Krasovsky_1940_Albers.latitude_of_projection_origin=0
	Krasovsky_1940_Albers.false_easting = 0.0
	Krasovsky_1940_Albers.false_northing = 0.0
	Krasovsky_1940_Albers.datum = "D_Krasovsky_1940"
	Krasovsky_1940_Albers.spatial_ref=im_proj
	Krasovsky_1940_Albers.GeoTransform=im_geotrans
	Krasovsky_1940_Albers._CoordinateTransformType = "Projection";
	Krasovsky_1940_Albers._CoordinateSystems = "ProjectionCoordinateSystem"
	Krasovsky_1940_Albers._CoordinateAxisTypes = "GeoX GeoY"
	data_variable.grid_mapping = 'albers_conical_equal_area'
	#然后将数据存入该文件
	#左侧经度开始；im_geotrans[0]，长度：im_width
	lon=get_nparray_by_startNum_length_step(im_geotrans[0],len(lons),im_geotrans[1])
	lat=get_nparray_by_startNum_length_step(im_geotrans[3],len(lats),im_geotrans[-1])
	lats[:]=lat
	lons[:]=lon
	soil_layers_top[:]=np.array([0,5,20,50,100,150])
	soil_layers_bottom[:]=np.array([5,20,50,100,150,999])
	#soil_layers[:]=np.array([0,5,20,50,100,150])
	#soil_layers[:]=np.arange(len(soil_layers))
	#soil_layers[:]=np.array([0,5,20,50,100,150])
	for i in range(len(soil_layers_top)):
		if data_variable.ndim==2:
			rasterdata=read_rasterfile(final_file_list[i])[0]
			rasterdata[rasterdata==NoDataValue]=255
			data_variable[:]=rasterdata
		elif data_variable.ndim>2:
			rasterdata=read_rasterfile(final_file_list[i])[0]
			rasterdata[rasterdata==NoDataValue]=255
			data_variable[i]=rasterdata
		print_current_process(round((i+1)*1.0/len(soil_layers_top),2))
	#data_variable.scale_factor=Scale_factor
	#添加全局信息
	mydataset.title=filetitle
	mydataset.source='Ecohydrology Lab'
	mydataset.references='Ecohydrology Lab'
	mydataset.institution='Institude of Global environment change'
	mydataset.history= 'Created ' + time.ctime(time.time())
	mydataset.setncattr('begin_date_(yyyy-mm)','1992-01')
	mydataset.setncattr('end_date_(yyyy-mm)','2050-12')
	mydataset.setncattr('begin_date_(yyyy)', begin_end_date[0])
	mydataset.setncattr('end_date_(yyyy)', begin_end_date[1])
	#mydataset.begin_date= "2016"
 	#mydataset.end_date = begin_end_date[1]
 	#mydataset.SOUTH_WEST_CORNER_LAT=36.0
 	#mydataset.SOUTH_WEST_CORNER_LON = 105.0
	mydataset.setncattr('upleft(x,y) center',[-852356.7,4523980.0])
	mydataset.setncattr('pixsize(x,y)',[1000.0,1000.0])
	mydataset.close()



input_rasterfiles_path=u"F:/2GEMS/1GEMS_Model_Data_YRB/1YRB_data/raster_data/soil/soil_clay/0_updated_clay/6层数据/"
input_rasterfiles_type='.tif'
output_NC_file=u"F:/2GEMS/1GEMS_Model_Data_YRB/1YRB_data/nc4_data/2.soil/3.soil_clay_updated/soil_clay_yrb1km_dv01.nc4"
data_variable_name='clay_ttl'
data_variable_longname="total clay"
all_vars_units=['cm','meters','meters','percent']
#begin_end_date=["0cm","200cm"]
begin_end_date=["1992","2015"]
#Fill_value=-3.40282346639e+038
filetitle="Gridded Total Clay for 6 Soil Layers"
Scale_factor=1
create_NC_by_raster(input_rasterfiles_path,input_rasterfiles_type,output_NC_file,data_variable_name,data_variable_longname,all_vars_units,Scale_factor,begin_end_date,filetitle)