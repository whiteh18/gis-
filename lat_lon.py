import geojson

with open('name.geojson','r',encoding = 'utf-8') as file:
    geojson_content = file.read()
    geojson_obj = geojson.loads(geojson_content)


min_x_1 = 999
max_y_1 = 0
min_y_1 = 999
for i in geojson_obj['features'][0]['geometry']['coordinates']:
    length = len(i)
    for j in range(length):
        if float(i[j][0]) < min_x_1:
            min_x_1 = float(i[j][0])
        if float(i[j][1]) > max_y_1:
            max_y_1 = float(i[j][1])
        if float(i[j][1]) < min_y_1:
            min_y_1 = float(i[j][1])
        print(i[j])
print(min_x_1)
print(max_y_1)
print(min_y_1)
#找最小x和最大y
#用于计算缩放比例和平移量



max_y_2 = 0
min_y_2 = 999
for i in geojson_obj['features'][1]['geometry']['coordinates']:
    length = len(i)
    for j in range(length):
        if float(i[j][1]) > max_y_2:
            max_y_2 = float(i[j][1])
        if float(i[j][1]) < min_y_2:
            min_y_2 = float(i[j][1])
        print(i[j])
print(max_y_2)
print(min_y_2)
#同理第二个字



max_x = 0
max_y_3 = 0
min_y_3 = 999
for i in geojson_obj['features'][2]['geometry']['coordinates']:
    length = len(i)
    for j in range(length):
        if float(i[j][0]) > max_x:
            max_x = i[j][0]
        if int(i[j][1]) > max_y_3:
            max_y_3 = float(i[j][1])
        if float(i[j][1]) < min_y_3:
            min_y_3 = float(i[j][1])
        print(i[j])
print(max_x)
print(max_y_3)
print(min_y_3)
#同理第三个字

#上面三步我应该是写复杂了，可能因为当时是边探索数据边写的T_T
#可以直接将三者结合起来


Sx = (115.5 - 115.3) / (max_x - min_x_1)
Sy = (31.9 - 31.8) / (max(max_y_1,max_y_2,max_y_3) - min(min_y_1,min_y_2,min_y_3))
print(Sx,Sy)
#比例变换中的Sx和Sy
#115.5,115.3,31.9,31.8对应家乡经纬度范围


i = geojson_obj['features'][1]['geometry']['coordinates'][6]
Tx = 115.3 - i[1][0] * Sx
Ty = 31.9 - i[1][1] * Sy
print(Tx,Ty)
#平移变换中的Tx和Ty
#i[0][0]和i[0][1]代表姓名中间那个点，好人工的确定方法



scale_matrix = [
    [Sx,0,0],
    [0,Sy,0],
    [0,0,1]]
move_matrix = [
    [1,0,0],
    [0,1,0],
    [Tx,Ty,1]]
point = [3,11,1]#检验点
result = [0,0,0]
for i in range(3):
    for j in range(3):
        result[i] += (point[j] * scale_matrix[j][i])        
result1 = [0,0,0]
for i in range(3):
    for j in range(3):
        result1[i] += (result[j] * move_matrix[j][i]) 
#用其中的一个点进行先缩放再平移的实验,观察结果是否合理
#result1 的值为 [115.19743589743588, 31.92325581395349, 1.0] ，后续还要将1去除

#-------------------------------------------------jupyter notebook写的,可能有些内容有重复----------------------------------------
import json
with open('my_name_exp.geojson', 'r',encoding = 'utf-8') as f:
    geojson_obj = json.load(f)

# 遍历每个特征并对多线的各个坐标进行变换
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
                result1 = [0, 0, 0]
                for x in range(3):
                    for y in range(3):
                        result1[x] += (result[y] * move_matrix[y][x])
                line[j] = result1
                line[j].pop()

# 将变换后的 GeoJSON 保存为新文件
with open('my_name_tra.geojson', 'w',encoding = 'utf-8') as f:
    json.dump(geojson_obj, f)
