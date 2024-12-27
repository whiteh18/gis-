#道格拉斯-普克压缩算法

#---------------------暂时的暂时的，要重新写-----------------------------------

#递归实现算法
import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import numpy as np
def doglas_pk(points,yuzhi):
    if len(points) < 3:
        return points
    dmax = 0.0
    index = 0
    for i in range(1,len(points) - 1):
        d = pl_distance(points[i],points[0],points[-1])
        if d > dmax:
            dmax = d
            index = i
    if dmax >= yuzhi:
        return doglas_pk(points[:index+1],yuzhi)[:-1] + doglas_pk(points[index:],yuzhi)
    else:
        return [points[0],points[-1]]
def pl_distance(point,line_start,line_end):
    x,y = point
    x1,y1 = line_start
    x2,y2 = line_end
    fz = abs((y2 - y1) * x + (x1 - x2) * y + x2 * y1 - x1 * y2)
    fm = np.sqrt((y2 - y1) ** 2 + (x1 - x2)**2)
    if fm != 0:
        return fz/fm
    else:
        return (y-y2)

# 读取兰伯特投影后的中国版图数据
data = gpd.read_file('china_L.shp')


original_coords = list(data.iloc[0]['geometry'].exterior.coords)


yuzhi = 2000

compressed_geometries = []
for index, row in data.iterrows():
    geometry = row['geometry']
    if geometry.geom_type == 'MultiPolygon':
        polygons = list(geometry.geoms)
        for polygon in polygons:
            original_coords = list(polygon.exterior.coords)
            compressed_coords = doglas_pk(original_coords, yuzhi)
            if len(compressed_coords) >= 3:
                compressed_geometries.append(Polygon(compressed_coords))
    else:
        original_coords = list(geometry.exterior.coords)
        compressed_coords = doglas_pk(original_coords, yuzhi)
        if len(compressed_coords) >= 3:
            compressed_geometries.append(Polygon(compressed_coords))

compressed_data = gpd.GeoDataFrame(geometry=compressed_geometries)
original_data = data.copy()

# 创建原始和压缩后的 GeoDataFrame
#original_data = gpd.GeoDataFrame(geometry=[gpd.Polygon(original_coords)])
#compressed_data = gpd.GeoDataFrame(geometry=[gpd.Polygon(compressed_coords)])

# 绘制原始和压缩后的地图
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
original_data.plot(ax=axes[0])
axes[0].set_title('Original Map')
compressed_data.plot(ax=axes[1])
axes[1].set_title('Compressed Map')
plt.show()
