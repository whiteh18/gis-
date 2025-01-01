import numpy as np
import matplotlib.pyplot as plt

# 兰伯特投影参数
lat_1 = 20  # 标准纬线1
lat_2 = 40  # 标准纬线2
lat_0 = 0  # 中心纬线
lon_0 = 105  # 中心经线

# 参考椭球体为克拉索夫斯基椭球体
a = 6378245  # 长半轴
f = 1 / 298.3  # 扁率
e = np.sqrt(f * (2 - f))  # 第一偏心率

# 清空画布，可无
def clear_canvas(ax):
    ax.clear()

# 计算两点之间的距离，数学公式
def distance(a, b):
    return np.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)

# 计算点到直线的距离 a和b为线上两点，c为要计算到直线距离的点，使用数学公式计算
def point_line_distance(a, b, c):
    A = (b[1] - a[1]) / np.sqrt((b[1] - a[1]) ** 2 + (b[0] - a[0]) ** 2)
    B = (a[0] - b[0]) / np.sqrt((b[1] - a[1]) ** 2 + (b[0] - a[0]) ** 2)
    C = (a[1] * b[0] - a[0] * b[1]) / np.sqrt((b[1] - a[1]) ** 2 + (b[0] - a[0]) ** 2)
    d = np.abs(A * c[0] + B * c[1] + C)
    return d

# 利用道格拉斯普克算法的函数进行矢量数据压缩
# data和原始.gen格式数据形式几乎一样但是包含序号和经过兰伯特投影变换的坐标的列表
# [[序号],[[经度],[纬度]],[END]......[END][END]]
# threshold为阈值，可以表示压缩率的大小
def compress_map(data, threshold):
    #兰伯特下
    city = []
    #压缩后
    compress_data = []
    #记录当前压缩处理过的点的数量
    cur_number = 0
    for i in range(len(data) - 1):
        line = data[i]
        if len(line) == 0:
            continue
        if line == "END" and data[i - 1] == "END":
            #读到连续的第两个END，则说明数据处理完毕
            break
        #读到END的时候，这个序号下的坐标添加结束，对这个面中的边轮廓的点要素进行分析处理
        if line == "END":
            max_dist = 0
            min_dist = distance(city[0], city[1])
            #用来记录
            k = 0
            for j in range(len(city)):
                d = distance(city[0], city[j])
                if d > max_dist:
                    max_dist = d
                    k = j
                if d < min_dist:
                    min_dist = d
            if city[0][0] == city[-1][0] or city[0][1] == city[-1][1]:
                city1 = city[:k]
                city2 = city[k:]
                douglas_peucker(city1, threshold, compress_data, cur_number)
                douglas_peucker(city2, threshold, compress_data, cur_number)
            else:
                douglas_peucker(city, threshold, compress_data, cur_number)
            city = []
            compress_data.append("END")
            continue
        if len(line) == 1:
            compress_data.append(line)
            #序号不做处理，直接添加到compress_data列表中
            continue
        #添加序号后,就在city列表中添加属于该序号下的众多点的坐标
        city.append(line)
    return compress_data

# Douglas-Peucker算法 递归
def douglas_peucker(data, threshold, compress_data, cur_number):
    #若只剩下两个点，结束递归
    if len(data) == 2:
        compress_data.append(data[0])
        cur_number += 1
        compress_data.append(data[1])
        cur_number += 1
        return
    #只剩下一个点
    if len(data) == 1:
        compress_data.append(data[0])
        cur_number += 1
        return
    #为空
    if len(data) == 0:
        return
    #用来记录离特定线最远点的下标
    k = 0
    max_dist = 0
    for i in range(1, len(data) - 1):
        d = point_line_distance(data[0], data[-1], data[i])
        if d > max_dist:
            max_dist = d
            k = i
    #最大距离小于等于阈值则这条线之间的点都舍去
    if max_dist <= threshold:
        compress_data.append(data[0])
        cur_number += 1
        compress_data.append(data[-1])
        cur_number += 1
        return
    #最大距离大于阈值，从该点分两条路，递归进行判断
    data1 = data[:k]
    data2 = data[k:]
    douglas_peucker(data1, threshold, compress_data, cur_number)
    douglas_peucker(data2, threshold, compress_data, cur_number)

# 显示数据
def show_data(data, ax):
    x_points = []
    y_points = []
    count = 0
    for line in data:
        if line == "":
            continue
        if line == "END":
            ax.plot(x_points, y_points, color='black')
            x_points.clear()
            y_points.clear()
            count = 0
            continue
        if len(line) < 2:
            count = 1
            continue
        if count == 1:
            #t = get_point(line[0], line[1])
            #x_points.append(t[0])
            #y_points.append(t[1])
            x_points.append(line[0])
            y_points.append(line[1])
            count += 1
        else:
            #t = get_point(line[0], line[1])
            #x_points.append(t[0])
            #y_points.append(t[1])
            x_points.append(line[0])
            y_points.append(line[1])

# 转换为兰伯特投影
def to_lambert(data):
    #存放投影转换后的数据
    new_data = []
    #记录成功转换的点数
    pre_number = 0
    data = data.split('\n')  # 使用换行度分割数据，保存在列表中
    
    for line in data:
        line = line.strip()  # 遍历每一行数据，并去除行首和行尾的空白字符和换行符
        if line == "END":
            new_data.append("END")
            #如果当前行使"END",则将直接添加到new_data中，并跳过后续处理
            continue
        point = line.split(',')
        if len(point) < 2:
            new_data.append(point)
            #使标记为一个面的序号不参与下面后续的坐标处理
            continue
        try:
            # 尝试兰伯特投影转换点的坐标
            t = lambert_point(float(point[0]), float(point[1]))
            pre_number += 1
            new_data.append(t)
        except ValueError as e:
            print(f"无法转换数据: {line}, 错误: {e}")
    return new_data

# 兰伯特投影转换函数针对点进行坐标变换
def lambert_point(lon, lat):
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

    return [x, y]

# 获取点,对于要将地理坐标转换为画布坐标的
#但是这样还要修改兰伯特投影中的内容将y的计算公式中整体加负号
def get_point(x, y):
    x0 = (x - 5000000) / 10100 + 1100
    y0 = 570 - (y - 1000000) / 10100 - 55
    return [x0, y0]

# 显示地图
def show_map(file, threshold, ax):
    with open(file, 'r') as f:
        data = f.read()
    clear_canvas(ax)#可无
    new_data = to_lambert(data)
    if threshold > 0:
        compress_data = compress_map(new_data, threshold)
        show_data(compress_data, ax)
    else:
        show_data(new_data, ax)

# 主函数
def main():
    file_path = 'CHINA_Arc.gen' 
    threshold = 1800  # 设置压缩率为50%
    
    fig, ax = plt.subplots(figsize=(10, 10))
    show_map(file_path, threshold, ax)
    plt.show()

if __name__ == "__main__":
    main()
