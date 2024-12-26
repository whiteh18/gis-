#--------------1.1 + 经纬度格网-----------------

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
data = gpd.read_file('china.shp')

# 创建一个图形和轴对象
fig, ax = plt.subplots()

# 绘制经纬网格
lat_range = np.arange(15, 60, 5)
lon_range = np.arange(70, 145, 5)

for lat in lat_range:
    plt.plot([70, 140], [lat, lat], color='gray', linestyle='--', linewidth=0.5)
for lon in lon_range:
    plt.plot([lon, lon], [15, 55], color='gray', linestyle='--', linewidth=0.5)
# 设置轴标签
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
# 显示图形
# 绘制中国版图
data.plot(ax=ax)
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------

#---------------------1.2 + 经纬度格网-----------------------------

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon

# 读取北京54坐标系下的中国版图数据
data = gpd.read_file('china.shp')

# 兰伯特投影参数
lat_1 = 20  # 标准纬线1
lat_2 = 40  # 标准纬线2
lat_0 = 0  # 中心纬线
lon_0 = 105  # 中心经线

# 参考椭球体 克拉索夫斯基椭球体
a = 6378245  # 长半轴
f = 1 / 298.3  # 扁率
e = np.sqrt(f * (2 - f))  # 第一偏心率

#根据投影变换的原理公式书写函数
def lambert_conversion(lon, lat):
    # 将经纬度转换为弧度
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat_1_rad = np.radians(lat_1)
    lat_2_rad = np.radians(lat_2)
    lat_0_rad = np.radians(lat_0)
    lon_0_rad = np.radians(lon_0)

    # 计算投影公式中的常数
    n = (np.sin(lat_1_rad) + np.sin(lat_2_rad)) / 2
    F = (np.cos(lat_1_rad) + np.cos(lat_2_rad)) / 2
    rho_0 = a * F / n
    rho = a * n / (1 + n * np.sin(lat_rad) / np.cos(lat_rad))
    theta = n * (lon_rad - lon_0_rad)

    # 计算投影转换后的坐标
    x = rho * np.sin(theta)
    y = rho_0 - rho * np.cos(theta)

    return x, y


# 转换china.shp数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            lon, lat = point
            x, y = lambert_conversion(lon, lat)
            new_geom_coords.append((x, y))
        new_coords.append([new_geom_coords])  
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom.geoms:
            new_poly_coords = []
            for point in poly.exterior.coords:
                lon, lat = point
                x, y = lambert_conversion(lon, lat)
                new_poly_coords.append((x, y))
            new_geom_coords.append(new_poly_coords)
        new_coords.append(new_geom_coords)

# 创建新的GeoDataFrame
new_geometries = []
for coords in new_coords:
    if isinstance(coords[0][0], (list, tuple)) and all(isinstance(c, (list, tuple)) for c in coords[0]):
        new_geometries.append(MultiPolygon([Polygon(c) for c in coords]))
    else:
        new_geometries.append(Polygon(coords[0]))

new_data = gpd.GeoDataFrame(data.drop('geometry', axis=1), geometry=new_geometries)

# 绘制转换后的地图
fig, ax = plt.subplots(figsize=(10, 10))
new_data.plot(ax=ax)

# 添加经纬网格绘制
# 获取数据的经纬度范围
min_lon, min_lat, max_lon, max_lat = 70, 15, 140, 55
lon_values = np.arange(min_lon, max_lon + 5, 5)  # 每隔5度生成经度值
lat_values = np.arange(min_lat, max_lat + 5, 5)  # 每隔5度生成纬度值

# 将经纬度转换为兰伯特投影坐标并绘制经线
for lon in lon_values:
    x_coords = []
    y_coords = []
    for lat in lat_values:
        x, y = lambert_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

# 将经纬度转换为兰伯特投影坐标并绘制纬线
for lat in lat_values:
    x_coords = []
    y_coords = []
    for lon in lon_values:
        x, y = lambert_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

plt.title('Lambert Conformal Conic Projection')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------
#墨卡托投影转换同理

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon

# 读取北京54坐标系下的中国版图数据
data = gpd.read_file('china.shp')

# 墨卡托投影参数
lon_0 = 0  # 中心经线
lat_0 = 0  # 中心纬线

# 参考椭球体相关参数
a = 6378137.0  # 长半轴 WGS84椭球体长半轴
f = 1 / 298.257223563  # 扁率 WGS84椭球体扁率
e = np.sqrt(f * (2 - f))  # 第一偏心率

#定义投影函数
def mercator_conversion(lon, lat):
    #转为弧度
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    lon_0_rad = np.radians(lon_0)
    lat_0_rad = np.radians(lat_0)

    # 计算墨卡托投影坐标的x轴坐标（经度方向）
    x = a * (lon_rad - lon_0_rad)

    # 计算墨卡托投影坐标的y轴坐标（纬度方向）
    y = a * np.log(np.tan(np.pi / 4 + lat_rad / 2) * ((1 - e * np.sin(lat_rad)) / (1 + e * np.sin(lat_rad))) ** (e / 2))

    return x, y

# 转换china.shp数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            lon, lat = point
            x, y = mercator_conversion(lon, lat)
            new_geom_coords.append((x, y))
        new_coords.append([new_geom_coords])  
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom.geoms:
            new_poly_coords = []
            for point in poly.exterior.coords:
                lon, lat = point
                x, y = mercator_conversion(lon, lat)
                new_poly_coords.append((x, y))
            new_geom_coords.append(new_poly_coords)
        new_coords.append(new_geom_coords)

# 创建新的GeoDataFrame
new_geometries = []
for coords in new_coords:
    if isinstance(coords[0][0], (list, tuple)) and all(isinstance(c, (list, tuple)) for c in coords[0]):
        new_geometries.append(MultiPolygon([Polygon(c) for c in coords]))
    else:
        new_geometries.append(Polygon(coords[0]))

new_data = gpd.GeoDataFrame(data.drop('geometry', axis=1), geometry=new_geometries)

# 绘制转换后的地图
fig, ax = plt.subplots(figsize=(10, 10))
new_data.plot(ax=ax)

# 添加经纬网格绘制
# 获取数据的经纬度范围
min_lon, min_lat, max_lon, max_lat = 70, 15, 140, 55
lon_values = np.arange(min_lon, max_lon + 5, 5)  # 每隔5度生成经度值
lat_values = np.arange(min_lat, max_lat + 5, 5)  # 每隔5度生成纬度值

# 将经纬度转换为墨卡托投影坐标并绘制经线
for lon in lon_values:
    x_coords = []
    y_coords = []
    for lat in lat_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

# 将经纬度转换为墨卡托投影坐标并绘制纬线
for lat in lat_values:
    x_coords = []
    y_coords = []
    for lon in lon_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

plt.title('Mercator Projection')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------
#-----------------------2.1 + 经纬度格网-------------------------

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
data = gpd.read_file('world.shp')

# 创建一个图形和轴对象
fig, ax = plt.subplots()

# 绘制经纬网格
lat_range = np.arange(-90, 95, 5)
lon_range = np.arange(-180, 185, 5)

for lat in lat_range:
    plt.plot([-180, 180], [lat, lat], color='gray', linestyle='-', linewidth=0.5)
for lon in lon_range:
    plt.plot([lon, lon], [-90, 90], color='gray', linestyle='-', linewidth=0.5)
# 设置轴标签
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
# 显示图形
# 绘制中国版图
data.plot(ax=ax)
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------
#----------------------------------2.2 + 经纬度格网-------------------------
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon
#import cartopy.crs as ccrs

# 读取北京54坐标系下的中国版图数据
data = gpd.read_file('world.shp')

# 墨卡托投影参数
lon_0 = 0  # 中心经线
lat_0 = 0  # 中心纬线

def mercator_conversion(lon, lat):
    valid_lat = np.clip(lat, -85, 85)  # 限制纬度范围在[-85°, 85°] 因为tan的原因
    # 将经纬度转换为弧度
    lat_rad = np.radians(valid_lat)
    lon_rad = np.radians(lon)
    lon_0_rad = np.radians(lon_0)

    # 计算投影后坐标
    x = lon_rad - lon_0_rad
    y = np.log(np.tan(np.pi / 4 + lat_rad / 2))

    return x, y

# 转换world.shp数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            lon, lat = point
            x, y = mercator_conversion(lon, lat)
            new_geom_coords.append((x, y))
        new_coords.append([new_geom_coords])  
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom.geoms:
            new_poly_coords = []
            for point in poly.exterior.coords:
                lon, lat = point
                x, y = mercator_conversion(lon, lat)
                new_poly_coords.append((x, y))
            new_geom_coords.append(new_poly_coords)
        new_coords.append(new_geom_coords)

# 创建新的GeoDataFrame
new_geometries = []
for coords in new_coords:
    if isinstance(coords[0][0], (list, tuple)) and all(isinstance(c, (list, tuple)) for c in coords[0]):
        new_geometries.append(MultiPolygon([Polygon(c) for c in coords]))
    else:
        new_geometries.append(Polygon(coords[0]))

new_data = gpd.GeoDataFrame(data.drop('geometry', axis=1), geometry=new_geometries)
# 绘制转换后的地图
fig, ax = plt.subplots(figsize=(10, 10))
new_data.plot(ax=ax)

# 添加经纬网格绘制
# 获取数据的经纬度范围
min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
lon_values = np.arange(min_lon, max_lon + 5, 5)  # 每隔5度生成经度值
lat_values = np.arange(min_lat, max_lat + 5, 5)  # 每隔5度生成纬度值

# 将经纬度转换为墨卡托投影坐标并绘制经线
for lon in lon_values:
    x_coords = []
    y_coords = []
    for lat in lat_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

# 将经纬度转换为墨卡托投影坐标并绘制纬线
for lat in lat_values:
    x_coords = []
    y_coords = []
    for lon in lon_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

plt.title('World Mercator Projection')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')

plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------
#------------------------------------2.3 + 经纬度格网--------------------------------

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon
import cartopy.crs as ccrs
import math

# 读取北京54坐标系下的世界版图数据
data = gpd.read_file('world.shp')

# 墨卡托投影参数
lon_0 = 0  # 中心经线
lat_0 = 0  # 中心纬线

def mercator_conversion(lon, lat):
    valid_lat = np.clip(lat, -85, 85)  # 限制纬度范围在[-85°, 85°] 因为tan的原因
    # 将经纬度转换为弧度
    lat_rad = np.radians(valid_lat)
    lon_rad = np.radians(lon)
    lon_0_rad = np.radians(lon_0)

    # 计算投影后坐标
    x = lon_rad - lon_0_rad
    y = np.log(np.tan(np.pi / 4 + lat_rad / 2))

    return x, y

# 转换数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            lon, lat = point
            x, y = mercator_conversion(lon, lat)
            new_geom_coords.append((x, y))
        new_coords.append([new_geom_coords])  
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom.geoms:
            new_poly_coords = []
            for point in poly.exterior.coords:
                lon, lat = point
                x, y = mercator_conversion(lon, lat)
                new_poly_coords.append((x, y))
            new_geom_coords.append(new_poly_coords)
        new_coords.append(new_geom_coords)

# 创建新的 GeoDataFrame
new_geometries = []
for coords in new_coords:
    if isinstance(coords[0][0], (list, tuple)) and all(isinstance(c, (list, tuple)) for c in coords[0]):
        new_geometries.append(MultiPolygon([Polygon(c) for c in coords]))
    else:
        new_geometries.append(Polygon(coords[0]))

new_data = gpd.GeoDataFrame(data.drop('geometry', axis=1), geometry=new_geometries)
# 绘制转换后的地图
fig, ax = plt.subplots(figsize=(10, 10))
new_data.plot(ax=ax)

# 添加经纬网格绘制
# 获取数据的经纬度范围
min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
lon_values = np.arange(min_lon, max_lon + 5, 5)  # 每隔5度生成经度值
lat_values = np.arange(min_lat, max_lat + 5, 5)  # 每隔5度生成纬度值

# 将经纬度转换为墨卡托投影坐标并绘制经线
for lon in lon_values:
    x_coords = []
    y_coords = []
    for lat in lat_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)

# 将经纬度转换为墨卡托投影坐标并绘制纬线
for lat in lat_values:
    x_coords = []
    y_coords = []
    for lon in lon_values:
        x, y = mercator_conversion(lon, lat)
        x_coords.append(x)
        y_coords.append(y)
    ax.plot(x_coords, y_coords, color='gray', linestyle='--', linewidth=0.5)


# 北京和巴黎的经纬度
beijing_lat = 39.9333  # 北纬39°56′，转换为小数形式
beijing_lon = 116.3833  # 东经116°20′，转换为小数形式
paris_lat = 48.85  # 北纬48°51′，转换为小数形式
paris_lon = 2.3333  # 东经2°20′，转换为小数形式

# 将经纬度转换为弧度
beijing_lat_rad = np.radians(beijing_lat)
beijing_lon_rad = np.radians(beijing_lon)
paris_lat_rad = np.radians(paris_lat)
paris_lon_rad = np.radians(paris_lon)

# 计算大圆航线（基于球面三角学原理）
# 计算两点间的夹角（角距离）
sigma = np.arccos(np.sin(beijing_lat_rad) * np.sin(paris_lat_rad) +
                  np.cos(beijing_lat_rad) * np.cos(paris_lat_rad) * np.cos(paris_lon_rad - beijing_lon_rad))

# 设定路径点数量
num_points = 100
path_points = []
for i in range(num_points):
    fraction = i / (num_points - 1)
    A = np.sin((1 - fraction) * sigma) / np.sin(sigma)
    B = np.sin(fraction * sigma) / np.sin(sigma)
    x = A * np.cos(beijing_lat_rad) * np.cos(beijing_lon_rad) + B * np.cos(paris_lat_rad) * np.cos(paris_lon_rad)
    y = A * np.cos(beijing_lat_rad) * np.sin(beijing_lon_rad) + B * np.cos(paris_lat_rad) * np.sin(paris_lon_rad)
    z = A * np.sin(beijing_lat_rad) + B * np.sin(paris_lat_rad)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, np.sqrt(x ** 2 + y ** 2))
    path_points.append((np.degrees(lon), np.degrees(lat)))

# 转换大圆航线路径点的坐标到墨卡托投影坐标系
path_points_projected = [mercator_conversion(lon, lat) for lon, lat in path_points]


# 绘制大圆航线
x_path, y_path = zip(*path_points_projected)
ax.plot(x_path, y_path, color='red', linewidth=2)

plt.title('World Mercator Projection with Great Circle Route from Beijing to Paris')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')

plt.show()
