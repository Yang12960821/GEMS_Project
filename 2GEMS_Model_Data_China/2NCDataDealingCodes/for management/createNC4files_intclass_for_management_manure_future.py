# -*- coding: utf-8 -*-
"""
create NC4files
nc4文件数据类型
_nptonctype  = {'S1' : NC_CHAR,
                'i1' : NC_BYTE,
                'u1' : NC_UBYTE,
                'i2' : NC_SHORT,
                'u2' : NC_USHORT,
                'i4' : NC_INT,
                'u4' : NC_UINT,
                'i8' : NC_INT64,
                'u8' : NC_UINT64,
                'f4' : NC_FLOAT,
                'f8' : NC_DOUBLE}
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
	mydataset.createDimension('time',len(final_file_list))
	mydataset.createDimension('y',im_height)
	mydataset.createDimension('x',im_width)
	#char类型投影参数变量
	#mydataset.createDimension('Krasovsky_1940_Albers',1)
	#构建变量
	#times=mydataset.createVariable('time',np.int16,('time',))
	times=mydataset.createVariable('time','u4',('time',))
	lats=mydataset.createVariable('y',np.float32,('y',))
	lons=mydataset.createVariable('x',np.float32,('x',))
	#投影变量
	Krasovsky_1940_Albers=mydataset.createVariable('albers_conical_equal_area','S1')
	#再存储数据的时候，需要先存放纬度，再存放精度
	if len(times)>1:
		#data_variable=mydataset.createVariable(data_variable_name,np.int16,('time','y','x'),fill_value=np.ubyte(Fill_value))
		data_variable=mydataset.createVariable(data_variable_name,'u2',('time','y','x'),fill_value=-9999)
	else:
		#data_variable=mydataset.createVariable(data_variable_name,np.int16,('y','x'),fill_value=np.ubyte(Fill_value))
		data_variable=mydataset.createVariable(data_variable_name,'u2',('y','x'),fill_value=-9999)
	#对变量的单位进行设定
	data_variable.long_name=data_variable_longname
	data_variable.units=all_vars_units[3]
	#变量的相关参数
	times.standard_name='time'
	times.long_time='time'
	times.units=all_vars_units[0]
	times.axis='T'
	lats.long_name="y distance on the projection plane from the UL corner"
	lats.axis = "Y"	
	lats.standard_name="projection_y_coordinate"
	lats.units=all_vars_units[1]
	lons.long_name="x distance on the projection plane from the UL corner"
	lons.axis = "X"
	lons.standard_name="projection_x_coordinate"
	lons.units=all_vars_units[2]
	Krasovsky_1940_Albers.grid_mapping_name = "albers_conical_equal_area"
	Krasovsky_1940_Albers.standard_parallel = (25.0,47.0)
	Krasovsky_1940_Albers.longitude_of_central_meridian=105.0
	#是参数latitude_of_projection_origin的问题，设置为0就正确了，不是经纬度范围的均值
	Krasovsky_1940_Albers.latitude_of_projection_origin=0
	Krasovsky_1940_Albers.false_easting = 0.0
	Krasovsky_1940_Albers.false_northing = 0.0
	Krasovsky_1940_Albers.semi_major_axis = 6378245.0
	Krasovsky_1940_Albers.semi_minor_axis = 6356863.018773047
	Krasovsky_1940_Albers.datum = "D_Krasovsky_1940"
	#Krasovsky_1940_Albers.spatial_ref=im_proj
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
	times[:]=np.arange(len(times))
	for i in range(len(times)):
		if data_variable.ndim==2:
			imdata=read_rasterfile(final_file_list[i])[0]
			imdata[imdata==NoDataValue]=-9999
			data_variable[:]=imdata
		elif data_variable.ndim>2:
			imdata=read_rasterfile(final_file_list[i])[0]
			imdata[imdata==NoDataValue]=-9999
			data_variable[i]=imdata
		print_current_process(round((i+1)*1.0/len(times),2))
	#data_variable[data_variable<=0]=255
	#data_variable.flag_values=np.arange(1,18,1).astype(np.ubyte)
	#data_variable.flag_values=np.arange(1,18,1)
	#data_variable.flag_meanings='water developed mechanically_disturbed_national_forest mechanically_disturbed_other_public mechanically_disturbed_private mining barren deciduous_forest evergreen_forest mixed_forest grassland shrubland agriculture hay_pasture herbaceous_wetland woody_wetland perennial_ice_snow'
	#data_variable[data_variable<0]=0
	#print data_variable.dtype
	#添加scale_factor
	#默认set_auto_scale是True
	#data_variable.set_auto_scale(True)
	#data_variable.scale_factor=Scale_factor
	#help(data_variable)
	#添加全局信息
	mydataset.title=filetitle
	mydataset.source='Ecohydrology Lab'
	mydataset.references='Ecohydrology Lab'
	mydataset.institution='Institude of Global environment change'
	mydataset.history= 'Created ' + time.ctime(time.time())
	mydataset.setncattr('begin_date (yyyy)',begin_end_date[0])
	mydataset.setncattr('end_date (yyyy)', begin_end_date[1])
	#mydataset.begin_date= "2016"
 	#mydataset.end_date = begin_end_date[1]
 	#mydataset.SOUTH_WEST_CORNER_LAT=36.0
 	#mydataset.SOUTH_WEST_CORNER_LON = 105.0
	mydataset.setncattr('upleft(x,y) center',[-852356.7,4523980.0])
	mydataset.setncattr('pixsize(x,y)',[1000.0,1000.0])
	mydataset.close()


#1个NC文件构建2个变量
def create_NC_by_raster_2_variables(input_rasterfiles_1_path,input_rasterfiles_1_type,input_rasterfiles_2_path,input_rasterfiles_2_type,output_NC_file,data_variable1_name,data_variable2_name,data_variable_1_longname,data_variable_2_longname,all_vars_2vars_units,Scale_factor1,Scale_factor2,begin_end_date,filetitle):
	final_file1_list=[]
	read_all_files_inside_filepath(input_rasterfiles_1_path,input_rasterfiles_1_type,final_file1_list)
	#final_file1_list=[final_file1_list1[-1]]*35
	#降雨、温度需要排序
	#final_file_list.sort(key=list_sort_by_data)
	if len(final_file1_list)>0:
		rasterinfo=read_rasterfile(final_file1_list[0])
		im_width,im_height,im_bands,im_geotrans,im_proj,NoDataValue1=rasterinfo[1],rasterinfo[2],rasterinfo[3],rasterinfo[4],rasterinfo[5],rasterinfo[6]
	else:
		print 'no \''+input_rasterfiles_1_type+'\' type of data in the path: '+input_rasterfiles_1_path
		exit()
	#Fill_value1=NoDataValue1
	#变量2参数，需要保证数据的基本参数都要一致
	final_file2_list=[]
	read_all_files_inside_filepath(input_rasterfiles_2_path,input_rasterfiles_2_type,final_file2_list)
	#final_file2_list=[final_file2_list1[-1]]*35
	#降雨、温度需要排序
	#final_file_list.sort(key=list_sort_by_data)
	if len(final_file2_list)>0:
		rasterinfo2=read_rasterfile(final_file2_list[0])
		im_width2,im_height2,im_bands2,im_geotrans2,im_proj2,NoDataValue2=rasterinfo2[1],rasterinfo2[2],rasterinfo2[3],rasterinfo2[4],rasterinfo2[5],rasterinfo2[6]
	else:
		print 'no \''+input_rasterfiles_2_type+'\' type of data in the path: '+input_rasterfiles_2_path
		exit()
	#Fill_value2=NoDataValue2
	mydataset=Dataset(output_NC_file,'w',format='NETCDF4')
	#构建维度
	mydataset.createDimension('time',len(final_file1_list))
	mydataset.createDimension('y',im_height)
	mydataset.createDimension('x',im_width)
	#char类型投影参数变量
	#mydataset.createDimension('Krasovsky_1940_Albers',1)
	#构建变量
	times=mydataset.createVariable('time','u4',('time',))
	lats=mydataset.createVariable('y',np.float32,('y',))
	lons=mydataset.createVariable('x',np.float32,('x',))
	#投影变量
	Krasovsky_1940_Albers=mydataset.createVariable('albers_conical_equal_area','S1')
	#再存储数据的时候，需要先存放纬度，再存放精度
	if len(times)>1:
		data_variable1=mydataset.createVariable(data_variable1_name,'u4',('time','y','x'),fill_value=9999)
		data_variable2=mydataset.createVariable(data_variable2_name,'u4',('time','y','x'),fill_value=9999)
	else:
		data_variable1=mydataset.createVariable(data_variable1_name,'u4',('y','x'),fill_value=9999)
		data_variable2=mydataset.createVariable(data_variable2_name,'u4',('y','x'),fill_value=9999)
	#对变量的单位进行设定
	times.units=all_vars_2vars_units[0]
	lats.units=all_vars_2vars_units[1]
	lons.units=all_vars_2vars_units[2]
	data_variable1.units=all_vars_2vars_units[3]
	data_variable2.units=all_vars_2vars_units[4]
	#变量的相关参数
	times.standard_name='time'
	times.long_time='time'
	times.axis='T'
	data_variable1.long_name=data_variable_1_longname
	data_variable2.long_name=data_variable_2_longname
	lats.standard_name="projection_y_coordinate"
	lats.axis = "Y"
	lons.standard_name="projection_x_coordinate"
	lons.axis = "X"
	Krasovsky_1940_Albers.grid_mapping_name = "albers_conical_equal_area"
	Krasovsky_1940_Albers.standard_parallel = (25.0,47.0)
	Krasovsky_1940_Albers.longitude_of_central_meridian=105.0
	#是参数latitude_of_projection_origin的问题，设置为0就正确了，不是经纬度范围的均值
	Krasovsky_1940_Albers.latitude_of_projection_origin=0
	Krasovsky_1940_Albers.false_easting = 0.0
	Krasovsky_1940_Albers.false_northing = 0.0
	Krasovsky_1940_Albers.semi_major_axis = 6378245.0
	Krasovsky_1940_Albers.semi_minor_axis = 6356863.018773047
	Krasovsky_1940_Albers.datum = "D_Krasovsky_1940"
	Krasovsky_1940_Albers.spatial_ref=im_proj
	Krasovsky_1940_Albers.GeoTransform=im_geotrans
	Krasovsky_1940_Albers._CoordinateTransformType = "Projection";
	Krasovsky_1940_Albers._CoordinateSystems = "ProjectionCoordinateSystem"
	Krasovsky_1940_Albers._CoordinateAxisTypes = "GeoX GeoY"
	data_variable1.grid_mapping = 'albers_conical_equal_area'
	data_variable2.grid_mapping = 'albers_conical_equal_area'
	#然后将数据存入该文件
	#左侧经度开始；im_geotrans[0]，长度：im_width
	lon=get_nparray_by_startNum_length_step(im_geotrans[0],len(lons),im_geotrans[1])
	lat=get_nparray_by_startNum_length_step(im_geotrans[3],len(lats),im_geotrans[-1])
	lats[:]=lat
	lons[:]=lon
	times[:]=np.arange(len(times))
	for i in range(len(times)):
		if data_variable1.ndim==2:
			imdata1=read_rasterfile(final_file1_list[i])[0]
			imdata2=read_rasterfile(final_file2_list[i])[0]
			imdata1[imdata1<0]=0
			imdata2[imdata2<0]=0
			data_variable1[:]=imdata1
			data_variable2[:]=imdata2
		elif data_variable1.ndim>2:
			imdata1=read_rasterfile(final_file1_list[i])[0]
			imdata2=read_rasterfile(final_file2_list[i])[0]
			imdata1[imdata1<0]=0
			imdata2[imdata2<0]=0
			data_variable1[i]=imdata1
			data_variable2[i]=imdata2
		print_current_process(round((i+1)*1.0/len(times),2))
	#添加scale_factor
	#默认set_auto_scale是True
	#data_variable.set_auto_scale(True)
	#data_variable1.scale_factor=Scale_factor1
	#data_variable2.scale_factor=Scale_factor2
	#help(data_variable)
	#添加全局信息
	mydataset.title=filetitle
	mydataset.source='Ecohydrology Lab'
	mydataset.references='Ecohydrology Lab'
	mydataset.institution='Institude of Global environment change'
	mydataset.history= 'Created ' + time.ctime(time.time())
	mydataset.setncattr('begin_date (yyyy)',begin_end_date[0])
	mydataset.setncattr('end_date (yyyy)', begin_end_date[1])
	#mydataset.begin_date= "2016"
 	#mydataset.end_date = begin_end_date[1]
 	#mydataset.SOUTH_WEST_CORNER_LAT=36.0
 	#mydataset.SOUTH_WEST_CORNER_LON = 105.0
	mydataset.setncattr('upleft(x,y) center',[-852356.7,4523980.0])
	mydataset.setncattr('pixsize(x,y)',[1000.0,1000.0])
	mydataset.close()


'''
input_rasterfiles_path=u"H:/GEMS_Model_data/my_studyarea_data/raster_data/management/fertilizer/YRB-Fertilizer_gN_m2_1km/"
input_rasterfiles_type='.tif'
output_NC_file=u"H:/GEMS_Model_data/my_studyarea_data/nc4_data/4.management/6.fertilizer/mgmt_fertilizer_yrb1km_1992_2015_dv01.nc4"
data_variable_name='fertilizer'
data_variable_longname='total fertilizer'
all_vars_units=['years since 1992-12','meters','meters','gN/m2/yr']
begin_end_date=["1992","2015"]
#Fill_value=0
filetitle='Land Carbon Derived Total Fertilizer'
Scale_factor=1
create_NC_by_raster(input_rasterfiles_path,input_rasterfiles_type,output_NC_file,data_variable_name,data_variable_longname,all_vars_units,Scale_factor,begin_end_date,filetitle)

'''


#构建2个变量的NC文件
#未来的这个数据全用历史数据的最后一个数据代替
input_rasterfiles_1_path=u"F:/2GEMS/1GEMS_Model_Data_YRB/1YRB_data/raster_data/management/fertilizer/0YRB_Manure_1km_update/Future/C_rate"
input_rasterfiles_1_type='.tif'
input_rasterfiles_2_path=u"F:/2GEMS/1GEMS_Model_Data_YRB/1YRB_data/raster_data/management/fertilizer/0YRB_Manure_1km_update/Future/N_rate"
input_rasterfiles_2_type='.tif'
output_NC_file=u"F:/2GEMS/1GEMS_Model_Data_YRB/1YRB_data/nc4_data/4.management/7.mgmt_manure/mgmt_manure_yrb1km_CGCM_rcp26_2016_2050_dv01.nc4"
data_variable1_name='C_rate'
data_variable2_name='N_rate'
#Fill_value1=-3.40282346639e+038
#Fill_value2=-3.40282346639e+038
data_variable_1_longname="rate of manure C applied"
data_variable_2_longname="rate of manure N applied"
all_vars_2vars_units=['years since 2016','meters','meters','gC/m^2','gN/m^2']
Scale_factor1=1
Scale_factor2=1
begin_end_date=["2016","2050"]
filetitle='Land Carbon Derived Total Fertilizer'
create_NC_by_raster_2_variables(input_rasterfiles_1_path,input_rasterfiles_1_type,input_rasterfiles_2_path,input_rasterfiles_2_type,output_NC_file,data_variable1_name,data_variable2_name,data_variable_1_longname,data_variable_2_longname,all_vars_2vars_units,Scale_factor1,Scale_factor2,begin_end_date,filetitle)