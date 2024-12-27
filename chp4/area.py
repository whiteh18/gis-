#空间度量之面积

#---------------------------1----------------------------------
#待完善可视化图像的细节！！
import geopandas as gpd
import matplotlib.pyplot as plt

data = gpd.read_file('Area.shp')
areas = data['geometry'].area

#print(data.head())

fig, ax = plt.subplots(figsize=(10, 10))
data.plot(ax=ax, column='AREA', legend=True)

plt.title('Area of each city in Jiangsu Province')
plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')

#---------------------------2----------------------------------

import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon,MultiPolygon


data = gpd.read_file('Area.shp')

# 地球半径
R = 6378137
#R = 6371000

def inverse_mercator_conversion(x, y):
    lon = np.degrees(x / R)
    lat = np.degrees(2 * np.arctan(np.exp(y / R)) - np.pi / 2)
    return lon, lat

# 转换数据的坐标
new_coords = []
for geom in data.geometry:
    new_geom_coords = []
    if geom.geom_type == 'Polygon':
        for point in geom.exterior.coords:
            x, y = point
            lon, lat = inverse_mercator_conversion(x, y)
            new_geom_coords.append((lon, lat))
        new_coords.append(new_geom_coords)
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom:
            new_poly_coords = []
            for point in poly.exterior.coords:
                x, y = point
                lon, lat = inverse_mercator_conversion(x, y)
                new_poly_coords.append((lon, lat))
            new_coords.append(new_poly_coords)

# 创建新的 GeoDataFrame
new_data = gpd.GeoDataFrame(data.drop('geometry', axis=1), geometry=[Polygon(coords) for coords in new_coords])
#保存为新文件 便于3的计算
new_data.to_file('inverse_mercator_js.shp')
# 绘制转换后的地图
fig, ax = plt.subplots()
new_data.plot(ax=ax)

plt.title('Inverse Mercator')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()

#---------------------------3----------------------------------

import geopandas as gpd
from math import radians, sin, cos
import matplotlib.pyplot as plt

def calculate_area_on_ellipsoid(polygon):
    area = 0
    for i in range(len(polygon)):
        j = (i + 1) % len(polygon)
        lat1, lon1 = polygon[i]
        lat2, lon2 = polygon[j]
        area += (radians(lon2) - radians(lon1)) * (2 + sin(radians(lat1)) + sin(radians(lat2)))
    return abs(area * 6371000.0 ** 2 / 2)

# 读取江苏省经纬度数据
data = gpd.read_file('inverse_mercator_js.shp')

# 计算面积并添加到新的列
for index, row in data.iterrows():
    geometry = row['geometry']
    if geometry.geom_type == 'Polygon':
        exterior_coords = list(geometry.exterior.coords)
        area = round(calculate_area_on_ellipsoid(exterior_coords)/1000000,2)
    elif geometry.geom_type == 'MultiPolygon':
        total_area = 0
        for poly in geometry:
            exterior_coords = list(poly.exterior.coords)
            total_area += calculate_area_on_ellipsoid(exterior_coords)
        area = round(total_area/1000000,2)
    data.at[index, 'Area'] = area

# 将更新后的数据保存为新的 shp 文件
#data.to_file('new_Area.shp')

# 生成面积专题地图，并添加面积数据标注
fig, ax = plt.subplots(figsize=(10, 8))
data.plot(column='Area', ax=ax, legend=True, cmap='YlOrRd', scheme='quantiles')

# 遍历每个几何对象，添加面积数据标注
for index, row in data.iterrows():
    geometry = row['geometry']
    if geometry.geom_type == 'Polygon':
        centroid = geometry.centroid
        ax.annotate(text=f"{row['Area']:.2f} km²", xy=(centroid.x, centroid.y),
                    ha='center', va='center', fontsize=8, color='black')
    elif geometry.geom_type == 'MultiPolygon':
        for poly in geometry:
            centroid = poly.centroid
            ax.annotate(text=f"{row['Area']:.2f} km²", xy=(centroid.x, centroid.y),
                        ha='center', va='center', fontsize=8, color='black')

plt.title('Area of Regions in Jiangsu Province')
plt.show()
