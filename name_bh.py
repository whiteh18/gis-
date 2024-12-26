#左斜
#虽然但是我觉得我这个左斜有点不太对，没有达到姓名下部分保持不变，只使姓名上部分左斜的效果
import json
def calculate_center(coordinates):
    min_x = float('inf')
    min_y = float('inf')
    max_x = float('-inf')
    max_y = float('-inf')
    for feature in coordinates:
        if feature['geometry']['type'] == 'MultiLineString':
            for line in feature['geometry']['coordinates']:
                for point in line:
                    x, y = point[0], point[1]
                    min_x = min(min_x, x)
                    max_x = max(max_x, x)
                    min_y = min(min_y, y)
                    max_y = max(max_y, y)
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    return center_x, center_y

#左斜,暂时考虑错切变换
b = 0.2
zuoxie_matrix = [
    [1,b,0],
    [0,1,0],
    [0,0,1]]

with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)
# 计算姓名中心坐标
center_x, center_y = calculate_center(geojson_obj['features'])  

#为了使姓名位置不变,以姓名的中心作为原点
#对应原理书上的相对于(Xf,Yf)的坐标变换

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * zuoxie_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('zuoxie.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------
#其实其他的变换只是改变了平面坐标变换矩阵

#耸肩
a = 1.2
b = 0.5
songjian_matrix = [
    [1,0,0],
    [0,a,b],
    [0,0,1]]
import json
with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * songjian_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('songjian.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#比例
Sx = 0.5
Sy = 0.5
scale_matrix = [
    [Sx,0,0],
    [0,Sy,0],
    [0,0,1]]
import json
with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * scale_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('scale.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#逆时针旋转30°
import math
cosa = math.cos(math.pi / 6)
sina = math.sin(math.pi / 6)
rotate_matrix = [
    [cosa,sina,0],
    [-sina,cosa,0],
    [0,0,1]]
import json
with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * rotate_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('ni_rotate.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)
#-----------------------------------------------分界线--------------------------------------------------

#顺时针旋转15°
import math
cosa = math.cos(math.pi / 12)
sina = math.sin(math.pi / 12)
rotate2_matrix = [
    [cosa,-sina,0],
    [sina,cosa,0],
    [0,0,1]]
import json
with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * rotate2_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('shun_rotate.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#平移
Tx = 0.5
Ty = 0.5
move_matrix = [
    [1,0,0],
    [0,1,0],
    [Tx,Ty,1]]
import json
with open('my_name_tra.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                # 先平移坐标，使其以中心为原点
                x, y = line[j][0], line[j][1]
                translated_x = x - center_x
                translated_y = y - center_y
                line[j] = [translated_x, translated_y, 1]

                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * move_matrix[n][m])
                result.pop()

                # 再将坐标平移回原来的位置
                result[0] += center_x
                result[1] += center_y
                line[j] = result[:2]

# 将变换后的 GeoJSON 保存为新文件
with open('move.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------
#最简单是使用matplotlib库进行变化前后对比可视化显示
#或者在geojson.io中输入geojson文件在线显示

#----------------显示在两幅图中---------
import matplotlib.pyplot as plt
import geopandas as gpd

gdf1 = gpd.read_file('my_name_tra.geojson')
gdf2 = gpd.read_file('zuoxie.geojson')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

gdf1.plot(ax = ax1, color='blue')
ax1.set_title('origin')

gdf2.plot(ax = ax2, color='green')
ax2.set_title('zuoxie')

plt.show()

#--------------显示在一幅图中----------------
import matplotlib.pyplot as plt
import geopandas as gpd

gdf1 = gpd.read_file('my_name_tra.geojson')
gdf2 = gpd.read_file('zuoxie.geojson')

ax = gdf1.plot(figsize=(10, 10), alpha=0.5, color='blue')
gdf2.plot(ax=ax, alpha=0.5, color='green')
plt.show()

#----------------还可以结合ipyleaflet对比显示-----------------------

import ipyleaflet as iplf
import geojson

default_center = [31.757,115.39]
default_zoom = 10

map = iplf.Map(center = default_center,zoom = default_zoom)
img_w = 'https://t0.tianditu.gov.cn/img_w/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0&LAYER=img&STYLE=default&TILEMATRIXSET=w&FORMAT=tiles&TILECOL={x}&TILEROW={y}&TILEMATRIX={z}&tk=61cd74f9e0a2c4cd8c750024caa46748'
cia_w = 'https://t0.tianditu.gov.cn/cia_w/wmts?SERVICE=WMTS&REQUEST=GetTile&VERSION=1.0.0&LAYER=cia&STYLE=default&TILEMATRIXSET=w&FORMAT=tiles&TILECOL={x}&TILEROW={y}&TILEMATRIX={z}&tk=61cd74f9e0a2c4cd8c750024caa46748'

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

def read_geojson_obj(filename):
    with open(filename,'r',encoding = 'utf-8') as f:
        geojson_obj = geojson.load(f)
        return geojson_obj

input_file_name = 'my_name_tra.geojson'
geojson_data = read_geojson_obj(input_file_name)

geo_layer = iplf.GeoJSON(data = geojson_data,style = {'color':'blue'})
map.add_layer(geo_layer)

input_file_name2 = 'zuoxie.geojson'
geojson_data2 = read_geojson_obj(input_file_name2)

geo_layer2 = iplf.GeoJSON(data = geojson_data2,style = {'color':'green'})
map.add_layer(geo_layer2)

map
