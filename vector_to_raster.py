#矢量数据转为栅格数据

#-------------------点转栅格------------------------------
import geopandas as gpd
import rasterio 
import numpy as np

#读取矢量点数据
try:
    points = gpd.read_file('./shp_file/points.shp')
except Exception as e:
    print(f"读取矢量数据失败：{e}")
    exit()

#定义栅格参数
pixel_size = 10
bounds = points.total_bounds
#最左边的x坐标和最上边的y坐标
origin_x,origin_y = bounds[0],bounds[3]
#栅格的宽度和高度
width = int((bounds[2] - bounds[0]) / pixel_size)
height = int((bounds[3] - bounds[1]) / pixel_size)
#创建一个空栅格数组
raster = np.zeros((height, width), dtype=np.uint8)
#计算行列号并赋值
def coord_to_row_col(x,y,origin_x,origin_y,pixel_size):
    col = int((x - origin_x) / pixel_size) #+ 1
    row = int((origin_y - y) / pixel_size) #+ 1
    return row,col
for index,point in points.iterrows():
    x,y = point.geometry.x,point.geometry.y
    row,col = coord_to_row_col(x,y,origin_x,origin_y,pixel_size)
    #确保行列号在栅格范围内
    if 0 <= row < height and 0 <= col < width:
        raster[row,col] = 1

#创建栅格文件的仿射变换
transform = rasterio.transform.from_origin(origin_x, origin_y, pixel_size, pixel_size)
#保存栅格数据到文件
try:
    with rasterio.open('points.tif', 'w', driver='GTiff',
                   height=height, width=width, count=1, dtype=raster.dtype,
                   crs=points.crs, transform=transform) as dst:
        dst.write(raster, 1)
        print("栅格文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败：{e}")

#----------------------------------线转栅格1----------------------------------------------

import numpy as np
import geopandas as gpd
import rasterio

# 定义栅格化函数
def rasterize_line(points, pixel_size):
    # 获取点的行列坐标
    def coord_to_row_col(x, y, origin_x, origin_y, pixel_size):
        row = int((origin_y - y) / pixel_size)
        col = int((x - origin_x) / pixel_size)
        return row, col
    
    # 确定行和列的范围
    bounds = points.total_bounds
    origin_x, origin_y = bounds[0], bounds[3]
    width = int((bounds[2] - bounds[0]) / pixel_size)
    height = int((bounds[3] - bounds[1]) / pixel_size)
    
    # 创建一个空栅格数组
    raster = np.zeros((height, width), dtype=np.uint8)
    
    for index, line in points.iterrows():
        x1, y1 = line.geometry.coords[0]  # 线段的起点
        x2, y2 = line.geometry.coords[1]  # 线段的终点
        
        i1, j1 = coord_to_row_col(x1, y1, origin_x, origin_y, pixel_size)
        i2, j2 = coord_to_row_col(x2, y2, origin_x, origin_y, pixel_size)
        
        # 将端点涂黑
        if 0 <= i1 < height and 0 <= j1 < width:
            raster[i1, j1] = 1
        if 0 <= i2 < height and 0 <= j2 < width:
            raster[i2, j2] = 1

        # 计算行数差和列数差
        row_diff = abs(i2 - i1)
        col_diff = abs(j2 - j1)

        if row_diff > col_diff:
            # 按行划分
            for i in range(min(i1, i2), max(i1, i2) + 1):
                # 计算当前行的列坐标
                if j1 == j2:
                    j = j1  # 纵坐标相同，直接用j1或j2
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j] = 1
                else:
                    # 计算与直线的交点
                    t = (i - i1) / (i2 - i1) if (i2 != i1) else 0  # 避免除以零
                    x_intersect = x1 + t * (x2 - x1)
                    j = int((x_intersect - origin_x) / pixel_size)
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j] = 1
        else:
            # 按列划分
            for j in range(min(j1, j2), max(j1, j2) + 1):
                # 计算当前列的行坐标
                if i1 == i2:
                    i = i1  # 横坐标相同，直接用i1或i2
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j] = 1
                else:
                    # 计算与直线的交点
                    t = (j - j1) / (j2 - j1) if (j2 != j1) else 0  # 避免除以零
                    y_intersect = y1 + t * (y2 - y1)
                    i = int((origin_y - y_intersect) / pixel_size)
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j] = 1
    
    return raster


try:
    lines = gpd.read_file('./shp_file/Pybuild1.shp')
except Exception as e:
    print(f"读取矢量数据失败: {e}")
    exit()

# 定义栅格参数
pixel_size = 10

# 调用栅格化函数
raster_data = rasterize_line(lines, pixel_size)

# 保存栅格数据到文件
transform = rasterio.transform.from_origin(0, lines.total_bounds[3], pixel_size, pixel_size)
try:
    with rasterio.open('Pybuild1.tif', 'w', driver='GTiff',
                       height=raster_data.shape[0], width=raster_data.shape[1], count=1, dtype=raster_data.dtype,
                       crs=lines.crs, transform=transform) as dst:
        dst.write(raster_data, 1)
        print("线段栅格文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败: {e}")
#-----------------------------------线转栅格2-----------------------------------------

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import rasterize
from shapely.geometry import LineString

def vector_to_raster(vector_path, raster_path, pixel_size):
    # 读取矢量线数据
    gdf = gpd.read_file(vector_path)
    
    # 获取线的边界
    bounds = gdf.total_bounds
    min_x, min_y, max_x, max_y = bounds

    # 计算栅格的宽和高
    width = int((max_x - min_x) / pixel_size) + 1
    height = int((max_y - min_y) / pixel_size) + 1

    # 创建一个空栅格，值为0（背景）
    raster = np.zeros((height, width), dtype=np.uint8)

    # 将线段转化为栅格
    for geom in gdf.geometry:
        if isinstance(geom, LineString):
            # 获取线段的坐标
            x_coords, y_coords = geom.xy

            # 遍历线段上的每个点
            for i in range(len(x_coords) - 1):
                start = (x_coords[i], y_coords[i])
                end = (x_coords[i + 1], y_coords[i + 1])
                
                # 获取线段在栅格中的起点和终点索引
                start_row = int((max_y - start[1]) / pixel_size)
                start_col = int((start[0] - min_x) / pixel_size)
                end_row = int((max_y - end[1]) / pixel_size)
                end_col = int((end[0] - min_x) / pixel_size)

                # 使用 Bresenham 算法进行栅格化
                rasterize_line(raster, start_row, start_col, end_row, end_col)

    # 定义转化为栅格的变换
    transform = rasterio.transform.from_bounds(min_x, min_y, max_x, max_y, raster.shape[1], raster.shape[0])

    # 写入栅格文件
    with rasterio.open(raster_path, 'w', driver='GTiff', height=raster.shape[0], width=raster.shape[1],
                       count=1, dtype='uint8', crs=gdf.crs, transform=transform) as dst:
        dst.write(raster, 1)

def rasterize_line(raster, start_row, start_col, end_row, end_col):
    """实现 Bresenham 算法来填充线段的路径"""
    dx = end_col - start_col
    dy = end_row - start_row
    sx = 1 if dx > 0 else -1
    sy = 1 if dy > 0 else -1
    dx = abs(dx)
    dy = abs(dy)

    if dx > dy:
        err = dx / 2.0
        while start_col != end_col:
            if 0 <= start_row < raster.shape[0] and 0 <= start_col < raster.shape[1]:
                raster[start_row, start_col] = 1  # 填充当前点
            err -= dy
            if err < 0:
                start_row += sy
                err += dx
            start_col += sx
    else:
        err = dy / 2.0
        while start_row != end_row:
            if 0 <= start_row < raster.shape[0] and 0 <= start_col < raster.shape[1]:
                raster[start_row, start_col] = 1  # 填充当前点
            err -= dx
            if err < 0:
                start_col += sx
                err += dy
            start_row += sy

# 使用示例
vector_file_path = './shp_file/sh_prj6.shp'  
raster_file_path = './raster_file/sh_prj6.tif'  
pixel_size = 10  # 设置栅格的像素大小

vector_to_raster(vector_file_path, raster_file_path, pixel_size)
#-------------------------------------面转栅格---------------------------------------------
#不对不对 非常不对
import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import Polygon

def boundary_fill_raster(vector_path, raster_path, pixel_size, a):
    # 读取矢量多边形数据
    gdf = gpd.read_file(vector_path)

    # 获取多边形的总边界
    boundaries = gdf.boundary
    if boundaries.is_empty.any():
        raise ValueError("The provided geometries have no boundaries.")

    # 计算栅格的边界
    min_x, min_y, max_x, max_y = boundaries.total_bounds

    # 计算栅格的宽和高
    width = int((max_x - min_x) / pixel_size) + 1
    height = int((max_y - min_y) / pixel_size) + 1

    # 初始化栅格_array，值为0
    raster_array = np.zeros((height, width), dtype=np.int32)

    # 将边界坐标转换为行列索引
    def to_indices(x, y):
        col = int((x - min_x) / pixel_size)
        row = int((max_y - y) / pixel_size)
        return row, col

    # 遍历每个多边形
    for polygon in gdf.geometry:
        if isinstance(polygon, Polygon):
            exterior_coords = np.array(polygon.exterior.coords)

            # 处理多边形的边界
            for i in range(len(exterior_coords) - 1):
                start = exterior_coords[i]
                end = exterior_coords[i + 1]

                start_row, start_col = to_indices(*start)
                end_row, end_col = to_indices(*end)

                if start_row == end_row:  # 同一行
                    if end_col > start_col:  # 从左向右
                        for col in range(start_col, end_col + 1):
                            raster_array[start_row, col] += a
                    else:  # 从右向左
                        for col in range(end_col, start_col + 1):
                            raster_array[start_row, col] -= a

                else:  # 不在同一行
                    # 采用简单的线性插值来确定直线走向
                    dx = end_col - start_col
                    dy = end_row - start_row
                    steps = max(abs(dx), abs(dy))

                    x_inc = dx / steps
                    y_inc = dy / steps

                    x, y = start_col, start_row
                    for _ in range(steps + 1):
                        row_idx, col_idx = to_indices(x, y)
                        if row_idx >= 0 and row_idx < height and col_idx >= 0 and col_idx < width:
                            if dy >= 0:  # 下行，左侧加a
                                raster_array[row_idx, col_idx] += a
                            else:  # 上行，左侧减a
                                raster_array[row_idx, col_idx] -= a

                        x += x_inc
                        y += y_inc

    # 处理栅格数据保存
    transform = rasterio.transform.from_bounds(min_x, min_y, max_x, max_y, width, height)
    with rasterio.open(raster_path, 'w', driver='GTiff', height=raster_array.shape[0],
                       width=raster_array.shape[1], count=1, dtype='int32',
                       crs=gdf.crs, transform=transform) as dst:
        dst.write(raster_array, 1)

# 使用示例
vector_file_path = 'path/to/your/polygon_file.shp'  # 替换为你的多边形矢量文件路径
raster_file_path = 'path/to/your/output_raster.tif'  # 替换为你希望保存的栅格文件路径
pixel_size = 10  # 设置栅格的像素大小
a = 1  # 边界加权值

boundary_fill_raster(vector_file_path, raster_file_path, pixel_size, a)
