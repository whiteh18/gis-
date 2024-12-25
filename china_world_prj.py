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
    plt.plot([70, 140], [lat, lat], color='gray', linestyle='-', linewidth=0.5)
for lon in lon_range:
    plt.plot([lon, lon], [15, 55], color='gray', linestyle='-', linewidth=0.5)
# 设置轴标签，但是还是中文问题
#ax.set_xlabel('经度')
#ax.set_ylabel('纬度')
# 显示图形
# 绘制中国版图
data.plot(ax=ax)
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------

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


def lambert_conversion(lon, lat):
    # 将经纬度转换为弧度
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    lat_1_rad = np.radians(lat_1)
    lat_2_rad = np.radians(lat_2)
    lat_0_rad = np.radians(lat_0)
    lon_0_rad = np.radians(lon_0)

    # 计算常数
    n = (np.sin(lat_1_rad) + np.sin(lat_2_rad)) / 2
    F = (np.cos(lat_1_rad) + np.cos(lat_2_rad)) / 2
    rho_0 = a * F / n
    rho = a * n / (1 + n * np.sin(lat_rad) / np.cos(lat_rad))
    theta = n * (lon_rad - lon_0_rad)

    # 计算坐标
    x = rho * np.sin(theta)
    y = rho_0 - rho * np.cos(theta)

    return x, y


# 转换数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            lon, lat = point
            x, y = lambert_conversion(lon, lat)
            new_geom_coords.append((x, y))
        new_coords.append([new_geom_coords])  # Ensure it's a list of polygons
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
# 获取数据的大致经纬度范围（假设数据覆盖了一定的范围，这里简单示意，可以根据实际更精确调整）
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

# 参考椭球体相关参数（常用的地球参考椭球体参数，可根据实际需求调整）
a = 6378137.0  # 长半轴（WGS84椭球体长半轴，可根据实际参考椭球更换）
f = 1 / 298.257223563  # 扁率（WGS84椭球体扁率，可根据实际参考椭球更换）
e = np.sqrt(f * (2 - f))  # 第一偏心率

def mercator_conversion(lon, lat):
    """
    将经纬度坐标转换为墨卡托投影坐标
    """
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    lon_0_rad = np.radians(lon_0)
    lat_0_rad = np.radians(lat_0)

    # 计算墨卡托投影坐标的x轴坐标（经度方向）
    x = a * (lon_rad - lon_0_rad)

    # 计算墨卡托投影坐标的y轴坐标（纬度方向），用到了正切函数和对数函数的变换
    y = a * np.log(np.tan(np.pi / 4 + lat_rad / 2) * ((1 - e * np.sin(lat_rad)) / (1 + e * np.sin(lat_rad))) ** (e / 2))

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
        new_coords.append([new_geom_coords])  # Ensure it's a list of polygons
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
# 获取数据的大致经纬度范围（假设数据覆盖了一定的范围，这里简单示意，可以根据实际更精确调整）
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

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
data = gpd.read_file('world.shp')

# 创建一个图形和轴对象
fig, ax = plt.subplots()

# 绘制经纬网格
lat_range = np.arange(-90, 91, 5)
lon_range = np.arange(-180, 180, 5)

for lat in lat_range:
    plt.plot([-180, 180], [lat, lat], color='gray', linestyle='-', linewidth=0.5)
for lon in lon_range:
    plt.plot([lon, lon], [-90, 91], color='gray', linestyle='-', linewidth=0.5)
# 设置轴标签
#ax.set_xlabel('经度')
#ax.set_ylabel('纬度')
# 显示图形
# 绘制中国版图
data.plot(ax=ax)
plt.show()

#-------------------------------------------------------分界线--------------------------------------------------------

