#矢量数据转为栅格数据
#老规矩.gen格式文件已经先转换为.shp格式文件

#-------------------矢量点转栅格------------------------------
import geopandas as gpd
import rasterio 
import numpy as np

#读取矢量点数据
try:
    points = gpd.read_file('./shp/point.shp')
except Exception as e:
    print(f"读取矢量数据失败：{e}")
    exit()

#定义栅格参数
pixel_size = 1
bounds = points.total_bounds
#最左边的x坐标和最上边的y坐标
origin_x,origin_y = bounds[0],bounds[3]
#栅格的宽度和高度
width = int((bounds[2] - bounds[0]) / pixel_size) + 1
height = int((bounds[3] - bounds[1]) / pixel_size) + 1
#创建一个空栅格数组
raster = np.zeros((height, width), dtype=np.uint8)
#计算行列号并赋值
def coord_to_row_col(x,y,origin_x,origin_y,pixel_size):
    col = int((x - origin_x) / pixel_size) 
    row = int((origin_y - y) / pixel_size) 
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
    with rasterio.open('point.tif', 'w', driver='GTiff',
                   height=height, width=width, count=1, dtype=raster.dtype,
                   crs=points.crs, transform=transform) as dst:
        dst.write(raster, 1)
        print("栅格文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败：{e}")

#----------------------------------线转栅格1(八方向栅格化)----------------------------------------------

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
    width = int((bounds[2] - bounds[0]) / pixel_size) + 1
    height = int((bounds[3] - bounds[1]) / pixel_size) + 1
    
    # 创建一个空栅格数组
    raster = np.zeros((height, width), dtype=np.uint8)
    
    for index, line in points.iterrows():
        x1, y1 = line.geometry.coords[0]  # 线段的起点
        x2, y2 = line.geometry.coords[-1]  # 线段的终点
        
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
            # 当行数大于列数
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
            # 当列数大于行数
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
    lines = gpd.read_file('./shp/line.shp')
except Exception as e:
    print(f"读取矢量数据失败: {e}")
    exit()

# 定义栅格参数
pixel_size = 10

# 调用栅格化函数
raster_data = rasterize_line(lines, pixel_size)

# 保存栅格数据到文件
transform = rasterio.transform.from_origin(line.total_bounds[0], lines.total_bounds[3], pixel_size, pixel_size)
try:
    with rasterio.open('line.tif', 'w', driver='GTiff',
                       height=raster_data.shape[0], width=raster_data.shape[1], count=1, dtype=raster_data.dtype,
                       crs=lines.crs, transform=transform) as dst:
        dst.write(raster_data, 1)
        print("线栅格文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败: {e}")
#-----------------------------------线转栅格2（全路径栅格化）-----------------------------------------

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
    width = int((bounds[2] - bounds[0]) / pixel_size) + 1
    height = int((bounds[3] - bounds[1]) / pixel_size) + 1
    
    # 创建一个空栅格数组
    raster = np.zeros((height, width), dtype=np.uint8)
    
    for index, line in points.iterrows():
        x1, y1 = line.geometry.coords[0]  # 线段的起点
        x2, y2 = line.geometry.coords[-1]  # 线段的终点
        
        i1, j1 = coord_to_row_col(x1, y1, origin_x, origin_y, pixel_size)
        i2, j2 = coord_to_row_col(x2, y2, origin_x, origin_y, pixel_size)
        
        # 将端点涂黑
        if 0 <= i1 < height and 0 <= j1 < width:
            raster[i1, j1 ] = 1
        if 0 <= i2 < height and 0 <= j2 < width:
            raster[i2, j2 ] = 1

        # 计算行数差和列数差
        row_diff = abs(i2 - i1)
        col_diff = abs(j2 - j1)

        if row_diff > col_diff:
            ix = i1
            # 当行数大于列数
            # 按行分带，计算行号
            for i in range(min(i1, i2), max(i1, i2) + 1):
                if y2 > y1:
                    x1,y1,x2,y2 = x2,y2,x1,y1  # 交换起点终点
                if j1 == j2:
                    j = j1  # 纵坐标相同，直接用j1或j2
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j ] = 1
                else:
                    # 计算该列的起末和该直线交点的行号
                    k = (y2 - y1) / (x2 - x1)
                    for j in range(min(j1,j2), max(j1,j2) + 1):
                        
                        ia = int((-(j - 1) * pixel_size + (x1 - origin_x)) * k + (- y1 + origin_y) / pixel_size) 
                        ie = int((- j * pixel_size + (x1 - origin_x)) * k + (- y1 + origin_y) / pixel_size) 
                        if 0 <= j < width and 0 <= ia < height and 0 <= ie < height:
                            raster[ia - 1:ie + 1, j] = 1
                        if i == i1 if i1 < i2 else i2:
                            ix = ia if ia < ie else ie
            raster[ix,j1] = 1

        else:
            # 当列数大于行数
            # 按列分带 计算各个行起始列号和结尾列号
            jx = j1
            for j in range(min(j1, j2), max(j1, j2) + 1):
                if y2 > y1:
                    x1,y1,x2,y2 = x2,y2,x1,y1  # 交换起点终点

                if i1 == i2:
                    i = i1  # 横坐标相同，直接用i1或i2
                    if 0 <= i < height and 0 <= j < width:
                        raster[i, j ] = 1
                else:
                    k = (y2 - y1) / (x2 - x1)  # 斜率
                    for i in range(min(i1,i2), max(i1,i2) + 1):
                        ja = int(((origin_y - (i -1) * pixel_size - y1) / k + x1 - origin_x) / pixel_size) - 1
                        je = int(((origin_y - i * pixel_size - y1) / k + x1 - origin_x) / pixel_size) - 1
                        
                        if 0 <= i < height and 0 <= ja < width and 0 <= je < width:
                            if ja < je:
                                raster[i, ja - 1:je + 1] = 1
                            else:
                                raster[i, je - 1:ja + 1] = 1
                        if i == i1 if i1 < i2 else i2:
                            jx = ja if ja < je else je
            raster[i1,jx] = 1
    
    return raster


try:
    lines = gpd.read_file('./shp/line.shp')
except Exception as e:
    print(f"读取矢量数据失败: {e}")
    exit()

# 定义栅格参数
pixel_size = 1

# 调用栅格化函数
raster_data = rasterize_line(lines, pixel_size)

# 保存栅格数据到文件
transform = rasterio.transform.from_origin(lines.total_bounds[0],lines.total_bounds[3], pixel_size, pixel_size)
try:
    with rasterio.open('line2.tif', 'w', driver='GTiff',
                       height=raster_data.shape[0], width=raster_data.shape[1], count=1, dtype=raster_data.dtype,
                       crs=lines.crs, transform=transform) as dst:
        dst.write(raster_data, 1)
        print("线栅格文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败: {e}")
#-------------------------------------面转栅格---------------------------------------------

import numpy as np
import geopandas as gpd
import rasterio

# 坐标转换函数，将地理坐标转成行、列坐标
def geo_to_grid_coord(x, y, origin_x, origin_y, pixel_size):
    row = int((origin_y - y) / pixel_size)
    col = int((x - origin_x) / pixel_size)
    return row, col

# 边界代数法实现矢量面转栅格面
def vector_polygon_to_raster(polygons, pixel_size):
    bounds = polygons.total_bounds
    origin_x = bounds[0]
    origin_y = bounds[3]
    width = int((bounds[2] - bounds[0]) / pixel_size) + 1
    height = int((bounds[3] - bounds[1]) / pixel_size) + 1

    raster = np.zeros((height, width), dtype=np.int32)

    for index, polygon in polygons.iterrows():
        exterior_ring = polygon.geometry.exterior.coords
        num_vertices = len(exterior_ring)

        for i in range(num_vertices - 1):
            x1, y1 = exterior_ring[i]
            x2, y2 = exterior_ring[i + 1]

            row1, col1 = geo_to_grid_coord(x1, y1, origin_x, origin_y, pixel_size)
            row2, col2 = geo_to_grid_coord(x2, y2, origin_x, origin_y, pixel_size)

            # 确保行列索引大于零
            if not (0 <= row1 < height and 0 <= col1 < width and 0 <= row2 < height and 0 <= col2 < width):
                continue  # 如果坐标超出范围则跳过

            y_diff = abs(row2 - row1)
            x_diff = abs(col2 - col1)
            steep = y_diff > x_diff

            if steep:
                # 沿 y 方向绘制
                if row1 > row2:
                    row1, row2 = row2, row1
                    col1, col2 = col2, col1
                dx = col2 - col1
                dy = row2 - row1
                derr = abs(dy) * 2
                error = 0
                col = col1

                for row in range(row1, row2 + 1):
                    if 0 <= row < height and 0 <= col < width:
                        raster[row, col] = 1
                    error += derr
                    if error > dx:
                        col += 1 if col1 < col2 else -1  # 确保列的更新方向
                        error -= dx * 2
            else:
                # 沿 x 方向绘制
                if col1 > col2:
                    row1, row2 = row2, row1
                    col1, col2 = col2, col1
                dx = col2 - col1
                dy = row2 - row1
                derr = abs(dy) * 2
                error = 0
                row = row1

                for col in range(col1, col2 + 1):
                    if 0 <= row < height and 0 <= col < width:
                        raster[row, col] = 1
                    error += derr
                    if error > dx:
                        row += 1 if row1 < row2 else -1  # 确保行的更新方向
                        error -= dx * 2

    return raster


try:
    polygons = gpd.read_file('./shp/polygon.shp')
except Exception as e:
    print(f"读取矢量数据失败: {e}")
    exit()

pixel_size = 1
raster_data = vector_polygon_to_raster(polygons, pixel_size)

# 保存栅格数据到文件
transform = rasterio.transform.from_origin(polygons.total_bounds[0], polygons.total_bounds[3], pixel_size, pixel_size)
try:
    with rasterio.open('polygon2.tif', 'w', driver='GTiff',
                       height=raster_data.shape[0], width=raster_data.shape[1], count=1, dtype=raster_data.dtype,
                       crs=polygons.crs, transform=transform) as dst:
        dst.write(raster_data, 1)
        print("栅格面文件已成功生成！")
except Exception as e:
    print(f"保存栅格文件失败: {e}")
