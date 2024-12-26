#jupyter notebook中书写并运行

import matplotlib.pyplot as plt
import geopandas as gpd

input_file_name = 'my_name_tra.geojson'
gdf = gpd.read_file(input_file_name)

gdf.plot()
plt.show()
#可以先借助matplotlib库可视化姓名，检查geojson文件是否正确
#或者直接在geojson.io网页中输入geojson文件，在线在地图上显示姓名

#--------------------------------我是分界线------------------------------

import ipyleaflet as iplf
import geojson

default_center = [31.757,115.39]
default_zoom = 10

map = iplf.Map(center = default_center,zoom = default_zoom)
#申请成为天地图开发者，获得密钥，填写下面的img_w和cia_w
img_w = 'https://t0.tianditu.gov.cn/cia_w/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0&LAYER=cia&STYLE=default&TILEMATRIXSET=w&FORMAT=tiles&TILECOL={x}&TILEROW={y}&TILEMATRIX={z}&tk=您的密钥'
#卫星影像地图图层
cia_w = 'https://t0.tianditu.gov.cn/cia_w/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0&LAYER=cia&STYLE=default&TILEMATRIXSET=w&FORMAT=tiles&TILECOL={x}&TILEROW={y}&TILEMATRIX={z}&tk=您的密钥'
#区域属性字段的标记图层

img_w_layer = iplf.TileLayer(
    url = img_w,
    attribution = 'Map data © <a href="http://www.tianditu.gov.cn/">天地图</a>',
    name = '天地图卫星影像'
)
map.add_layer(img_w_layer)

cia_w_layer = iplf.TileLayer(
    url = cia_w,
    attribution = 'Map data © <a href="http://www.tianditu.gov.cn/">天地图</a>',
    name = '天地图卫星影像'
)
map.add_layer(cia_w_layer)

#要注意load和loads的区别哦
def read_geojson_obj(filename):
    with open(filename,'r',encoding = 'utf-8') as f:
        geojson_obj = geojson.load(f)
        return geojson_obj

input_file_name = 'my_name_tra.geojson'
geojson_data = read_geojson_obj(input_file_name)

geo_layer = iplf.GeoJSON(data = geojson_data)
map.add_layer(geo_layer)

map
