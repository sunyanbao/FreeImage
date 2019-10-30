#include  "stdafx.h"
#include "read_tiff.h"
#include "gdal\include\gdal_priv.h"
#include "gdal\include\ogr_spatialref.h"
#include "gdal\include\ogrsf_frmts.h"
#include "kernel_service\service\COM_Service.h"
/*#include "gdal\include\shapefil.h"  */ 
#include "shapefil.h"  

#include "FIIO_Mem.h"

#pragma comment(lib,"gdal_i.lib")
#pragma comment(lib,"shapelib_x64.lib")

#ifdef _DEBUG  
#pragma comment(lib, "FreeImaged.lib")  
#else  
#pragma comment(lib, "FreeImage.lib")  
#endif 


HandleTif::HandleTif()
{
	//tiff_path_ = "C:\\Users\\admin\\Desktop\\黄山ioc\\环保样例数据\\PM10.MOD021KM.A2019263_DAY_AVERAGE.tif";
	shp_path_ = "D:/river_shp/river.shp";
	//color_tex_path_ = "D:\\DeepEyev3\\localgis\\hsioc\\map\\color_index.png";
	//pic_path_ = "D:\\DeepEyev3\\localgis\\hsioc\\map\\heatmap.png";
	lla_map_.clear();
	gray_value_vec_.clear();
	tiff_wid_pixel_ = 0;
	tiff_hei_pixel_ = 0;
	min_gray_value_ = 0;
	max_gray_value_ = 0;
	min_lon_ = max_lon_=min_lat_=max_lat_=0;
}

HandleTif::~HandleTif()
{
}

bool HandleTif::clear()
{
	return true;
}


bool HandleTif::ReadTif(CString& tiff_path)
{
	//tif文件读取
	const char *charName = tiff_path.GetBuffer();
	//注册
	GDALAllRegister();
	//以防中文名不能正常读取
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	GDALDataset *m_pDataset = (GDALDataset*)GDALOpen(charName, GA_ReadOnly);
	if(m_pDataset==NULL)
	{
		AfxMessageBox("读取图像失败！");
		return FALSE;
	}

	tiff_wid_pixel_ = m_pDataset->GetRasterXSize(); //获取图像宽
	tiff_hei_pixel_ = m_pDataset->GetRasterYSize(); //获取图像高，像素单位，左上角像素是0，0
	//GT[0],GT[3]是做上角坐标(投影坐标)，GT[1],GT[5]是图像总向横向分辨率，就是每个像素代表的长度，GT[2],GT[4]是旋转相关，一般是0，不用管
	double adfGeoTransform[6];
	m_pDataset->GetGeoTransform(adfGeoTransform);
	//arcmap中对应米的单位
	double rbx,rby,ltx,lty;
	rbx = adfGeoTransform[0] + tiff_wid_pixel_ * adfGeoTransform[1];//右下角坐标
	rby = adfGeoTransform[3] + tiff_hei_pixel_ * adfGeoTransform[5];
	ltx = adfGeoTransform[0];
	lty = adfGeoTransform[3];//左上角坐标

	//投影坐标转经纬度
	OGRSpatialReference fRef,tRef;
	char *tmp = NULL;
	const char* projRef = m_pDataset->GetProjectionRef();
	//获得projRef的一份拷贝
	//由于projRef是const char*,下面的一个函数不接受，所以需要转换成非const
	tmp = (char *)malloc(strlen(projRef) + 1);
	strcpy_s(tmp, strlen(projRef)+1, projRef);

	//设置原始的坐标参数，和test.tif一致
	fRef.importFromWkt(&tmp);
	//设置转换后的坐标
	tRef.SetWellKnownGeogCS("WGS84");

	OGRCoordinateTransformation *coordTrans;
	coordTrans = OGRCreateCoordinateTransformation(&fRef, &tRef);
	//对应arcmap的十进制度
	coordTrans->Transform(1, &rbx, &rby);
	coordTrans->Transform(1, &ltx, &lty);
	max_lon_ = rbx;
	min_lat_ = rby;
	min_lon_ = ltx;
	max_lat_ = lty;

	//获取第1波段的波段指针，参数就是表示第几波段的意思  
	//读取数据
	GDALRasterBand  *m_pRasterBand = m_pDataset->GetRasterBand(1);
	double adfMinMax[2];
	int bGotMin,bGotMax;
	adfMinMax [0] = m_pRasterBand-> GetMinimum(&bGotMin);
	adfMinMax [1] = m_pRasterBand-> GetMaximum(&bGotMax);
	if(!(bGotMin && bGotMax))
		GDALComputeRasterMinMax((GDALRasterBandH)m_pRasterBand,TRUE,adfMinMax);
	min_gray_value_ = adfMinMax [0];
	max_gray_value_ = adfMinMax [1];

 	float *inBuf;
 	inBuf = new float[tiff_wid_pixel_*tiff_hei_pixel_];
 	ZeroMemory(inBuf,sizeof(float)*tiff_wid_pixel_*tiff_hei_pixel_);
 	CPLErr err;
 	//向inBuf读入数据
 	err=m_pRasterBand->RasterIO(GF_Read,0,0,tiff_wid_pixel_,tiff_hei_pixel_,inBuf,tiff_wid_pixel_,tiff_hei_pixel_,GDT_Float32,0,0);
 
 	if(err==CE_Failure)
 	{
 		AfxMessageBox("读取输入图像数据失败！");
 		return FALSE;
 	}
 	for (int j=0;j<tiff_hei_pixel_;j++)
 	{
 		for (int i=0;i<tiff_wid_pixel_;i++)
 		{
			float item = inBuf[j*tiff_wid_pixel_+i];
 			//if( _isnan(item)==0)//筛选有效值，去除黑色区域的点
 			gray_value_vec_.push_back(item);//存储灰度值
 		}
 	}

	delete inBuf;
	inBuf = NULL;
	GDALClose(m_pDataset);
	return true;
}

bool HandleTif::ReadShp()
{
	GDALAllRegister();
	GDALDataset *m_pDataset;
	CPLSetConfigOption("SHAPE_ENCODING","");  //解决中文乱码问题
	//读取shp文件
	m_pDataset = (GDALDataset*) GDALOpenEx(shp_path_.c_str(), GDAL_OF_VECTOR, NULL, NULL, NULL );

	if( m_pDataset == NULL )
	{
		AfxMessageBox("打开shp文件失败！");
		return FALSE;
	}

	OGRLayer  *m_pLayer;
	m_pLayer = m_pDataset->GetLayer(0); //读取层
	OGRFeature *m_pFeature;

	m_pLayer->ResetReading();
	double x,y,area;
	OGRPoint *point;
	OGREnvelope *envelope;
	OGRGeometry *m_pGeometry;
	point = new OGRPoint;
	envelope = new OGREnvelope;
	while( (m_pFeature = m_pLayer->GetNextFeature()) != NULL )
	{
		m_pGeometry = m_pFeature->GetGeometryRef();
		OGRwkbGeometryType pGeoType=m_pGeometry->getGeometryType();
		if( m_pGeometry != NULL && pGeoType == wkbPolygon )
		{
			OGRPolygon  *m_pPolygon = (OGRPolygon *) m_pGeometry;
			m_pPolygon->getEnvelope(envelope);
			double minX = envelope->MinX;
			double minY = envelope->MinY;
			double maxX = envelope->MaxX;
			double maxY = envelope->MaxY;
			for(double i=minX;i<maxX;i=i+0.0001){
				for(double j=minY;j<maxY;j=j+0.0001){
 					point->setX(i);
 					point->setY(j);
 					point->setM(0);
 					point->setZ(0);
					if(m_pPolygon->IsPointOnSurface(point)){
						lla_map_.insert(make_pair(i,j));
					}
				}
			}
// 			OGRLinearRing *pOGRLinearRing = m_pPolygon->getExteriorRing();
// 			int pointNum = pOGRLinearRing->getNumPoints();
// 			for(int i=0;i<pointNum;i++){
// 				pOGRLinearRing->getPoint(i,&point);//获取多边形内每个点的坐标
// 				i++;
// 				mapGeo.insert(pair<double,double>(point.getX(),point.getY()));
//  			}
		}
		else if( m_pGeometry != NULL && pGeoType == wkbMultiPolygon){
			OGRMultiPolygon *pMulPolygon=(OGRMultiPolygon*)m_pGeometry;
			OGRPolygon *m_pPolygon=NULL;
			for(int i=0;i<pMulPolygon->getNumGeometries();i++)
			{
				m_pPolygon=(OGRPolygon*)pMulPolygon->getGeometryRef(i);
// 				OGRLinearRing *pOGRLinearRing = m_pPolygon->getExteriorRing();
// 				int pointNum = pOGRLinearRing->getNumPoints();
// 				for(int i=0;i<pointNum;i++){
// 					pOGRLinearRing->getPoint(i,&point);//获取多边形内每个点的坐标
// 					i++;
// 					mapGeo.insert(pair<double,double>(point.getX(),point.getY()));
// 						}
				m_pPolygon->getEnvelope(envelope);
				double minX = envelope->MinX;
				double minY = envelope->MinY;
				double maxX = envelope->MaxX;
				double maxY = envelope->MaxY;
				for(double i=minX;i<maxX;i=i+0.0001){
					for(double j=minY;j<maxY;j=j+0.0001){
						point->setX(i);
						point->setY(j);
						point->setM(0);
						point->setZ(0);
						if(m_pPolygon->IsPointOnSurface(point)){
							lla_map_.insert(make_pair(i,j));
						}
					}
				}
			}
		OGRFeature::DestroyFeature( m_pFeature );
		}
	}
	delete envelope;
	envelope = NULL;
	GDALClose( m_pDataset );
	return true;
}

bool HandleTif::ConvertTiffToImage(CString& tiff_path,CString& image_path,CString& color_index_path)
{
	if(!ReadTif(tiff_path)){
		return false;
	}
// 	if(!ReadTif()||!ImportShpFile()){
// 		return false;
// 	}
	map<double,double>::iterator it;
	it = lla_map_.begin();
	int x,y;
	float lon,lat,gray_value;
	CString strContainer;
	CString strPara;
	float lon_per_pixel=0,lat_per_pixel=0;
	lon_per_pixel = (max_lon_-min_lon_)/tiff_wid_pixel_;
	lat_per_pixel = (max_lat_-min_lat_)/tiff_hei_pixel_;
	FIBITMAP* color_bitmap=FreeImage_Load(FIF_PNG,color_index_path.GetBuffer(),0);//加载图片
	int width = FreeImage_GetWidth(color_bitmap);
	int height = FreeImage_GetHeight(color_bitmap);
	RGBQUAD rgbvalue;
	rgbvalue.rgbBlue = 0;
	rgbvalue.rgbGreen = 0;
	rgbvalue.rgbRed = 255;
	rgbvalue.rgbReserved = 255;
	float ratio = 0;
	int delta_pixel = 1;
	FIBITMAP* bitmap=FreeImage_Allocate(tiff_wid_pixel_/delta_pixel,tiff_hei_pixel_/delta_pixel,32);//创建新图片;
	int index = 0;
	for (int j=0;j<tiff_hei_pixel_;j+=delta_pixel)
	{
		for (int i=0;i<tiff_wid_pixel_;i+=delta_pixel)
		{
// 			lon = lon_per_pixel*i+min_lon_;
// 			lat = -lat_per_pixel*j+max_lat_;
			gray_value = gray_value_vec_[j*tiff_wid_pixel_+i];
			ratio = (gray_value-min_gray_value_)/(max_gray_value_-min_gray_value_);
			index =  (width - 1)*(1-ratio);
			if (index<0)
			{
				index = 0;
			}
			else if (index>width)
			{
				index = width;
			}
			if (gray_value>min_gray_value_)//
			{
				FreeImage_GetPixelColor(color_bitmap,index,0,&rgbvalue);//获取像素颜色
				rgbvalue.rgbReserved = 255;
				FreeImage_SetPixelColor(bitmap,i,tiff_hei_pixel_-1-j,&rgbvalue);//写入像素颜色
			}
			//FreeImage_SetPixelColor(re,i,j,&color);//写入像素颜色
// 			strPara.Format("[%f,%f,%f],",lon,lat,weight);
// 			if( _isnan(weight)==0){
// 				strContainer += strPara;
// 			}
		}
	}
	FreeImage_Save(FIF_PNG,bitmap,image_path.GetBuffer(),0);//保存图片

// 	strContainer = strContainer.Left(strContainer.GetLength()-1);//去掉最后一个逗号
// 	CString strService;
// 	strService.Format("{\"service_id\":\"smart_city_effect.hot_area.create_real_hot_area\", \"diagonal_pos\":\
// 					  [\
// 					  [ %f, %f, 100 ],\
// 					  [ %f, %f, 100 ]\
// 					  ],\
// 					  \"hot_area_data\" :[%s],\
// 					  \"min_view_height\" : 10000\
// 					  }",min_lon_,max_lat_,max_lon_,min_lat_,strContainer);
// 	service::call(strService);
	return true;
}

bool HandleTif::ImportShpFile()
{
	SHPHandle hShpHandle = NULL;  // shp文件
	DBFHandle hDbfHandl  = NULL;  // dbf文件
	int     iFieldCount  = 0;     // 属性字段个数

	// 全部以只读方式打开文件
	hShpHandle = SHPOpen( shp_path_.c_str(), "rb");
	if ( !hShpHandle )
	{
		CString strMsg;
		strMsg.Format( _T("shp文件：%s，打开失败！请检查！"), shp_path_ );
		AfxMessageBox(strMsg);
		return false;
	}

	hDbfHandl = DBFOpen( shp_path_.c_str(), "rb" );
	if ( hDbfHandl )
	{
		// 读取属性字段
		iFieldCount = DBFGetFieldCount(hDbfHandl);
		DBFFieldType type;
		int iWidth = 0, iDecimals;
		iFieldCount = iFieldCount+1;
		int i=0;
	}

	// 读取数据
	int    nEntities = -1,
		   nShapeType= -1;
	SHPObject* pObj = NULL;
	SHPGetInfo( hShpHandle, &nEntities, &nShapeType, NULL, NULL );

	// 空间数据
	for (int i=0; i!=nEntities; ++i)
	{
		pObj = SHPReadObject( hShpHandle, i );
		if ( pObj==NULL ) 
			continue;

		CString sGraphType("");
		// 面
		if ( nShapeType==SHPT_POLYGON || nShapeType==SHPT_POLYGONZ || nShapeType==SHPT_POLYGONM )
		{
			sGraphType = "3";
			int nParts = pObj->nParts;
			int nPartCount = 0;
			for (int j=0; j!=nParts; ++j)
			{
				//vector<Vector3d> vecPos;
				if ( j==nParts-1 ) 
					nPartCount = pObj->nVertices;
				else
					nPartCount = pObj->panPartStart[j+1];
				int k = pObj->panPartStart[j];
				for (; k!=nPartCount; ++k)
				{
					lla_map_.insert(make_pair(pObj->padfX[k],pObj->padfY[k]));
					//Vector3d vPos;
// 					vPos.x = pObj->padfX[k];
// 					vPos.y = pObj->padfY[k];
// 					vPos.z = 9999.0;
// 					if ( nShapeType==SHPT_POLYGONZ )
// 						vPos.z = pObj->padfZ[k];
					//vecPos.push_back(vPos);
				}
			}
		}
		SHPDestroyObject(pObj);
	}
	return true;
}