#变换矩阵中的关键值需要再考量

#左斜
b = 0.2
zuoxie_matrix = [
    [1,b,0],
    [0,1,0],
    [0,0,1]]
import json
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与变换矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * zuoxie_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_zuoxie.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#耸肩
a = 1.2
b = 0.5
songjian_matrix = [
    [1,0,0],
    [0,a,b],
    [0,0,1]]
import json
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与变换矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * songjian_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_songjian.geojson', 'w',encoding = 'utf-8') as f:
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
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征将坐标与矩阵相乘
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * scale_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_scale.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#旋转
#逆时针旋转30°
import math
cosa = math.cos(math.pi / 6)
sina = math.sin(math.pi / 6)
rotate_matrix = [
    [cosa,sina,0],
    [-sina,cosa,0],
    [0,0,1]]
import json
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征并对多线坐标进行旋转
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * rotate_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_ni_rotate.geojson', 'w',encoding = 'utf-8') as f:
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
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征并对多线坐标进行旋转
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * rotate2_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_shun_rotate.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------

#平移
Tx = 1
Ty = 1
move_matrix = [
    [1,0,0],
    [0,1,0],
    [Tx,Ty,1]]
import json
with open('my_name.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征并对多线坐标进行平移
for feature in geojson_obj['features']:
    if feature['geometry']['type'] == 'MultiLineString':
        for line in feature['geometry']['coordinates']:
            length = len(line)
            for j in range(length):
                line[j].append(1)
                result = [0, 0, 0]
                for m in range(3):
                    for n in range(3):
                        result[m] += (line[j][n] * move_matrix[n][m])
                result.pop()        
                line[j] = result

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_move.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)

#-----------------------------------------------分界线--------------------------------------------------
#最简单是使用matplotlib库进行变化前后对比可视化显示
#不显示子啊地图上，是因为变换矩阵中的值选取的都比较大，在地图上变化较大，不易观察
#后续可调整

import matplotlib.pyplot as plt
import geopandas as gpd

gdf1 = gpd.read_file('my_name.geojson')
gdf2 = gpd.read_file('my_name_zuoxie.geojson')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

gdf1.plot(ax = ax1, color='blue')
ax1.set_title('origin')

gdf2.plot(ax = ax2, color='green')
ax2.set_title('zuoxie')

plt.show()
