# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:59:48 2024

@author: vicky
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from decimal import Decimal, ROUND_HALF_UP

# 讀取shapefiles
shapefile = r"C:\Users\AE133\Desktop\三維道路\test\AREA_RD_97221071_G.shp"
ouput_file_path = r'C:\Users\AE133\Desktop\三維道路\test\output'
#測試資料 A[0, 0, 2, 2] ,  B[6, 5, 4]
gdf = gpd.read_file(shapefile)

float_acc = 4 # If there are round off problem, using "1"
round_off = lambda num : float(Decimal(str(num)).quantize(Decimal('.{}'.format('0' * float_acc)), ROUND_HALF_UP))

#計算平均Z值
def avg_z(polygon):
    polygon_coords = list(polygon.exterior.coords)
    all_z = 0
    for i in range(len(polygon_coords) - 1):
        coord = polygon_coords[i]
        z = coord[2]
        all_z += z
    avg = all_z / (len(polygon_coords) - 1)
    return avg

if 'Diff' not in gdf.columns:
    gdf['Diff'] = np.nan

# 迭代每個polygon
for idx, row in gdf.iterrows():
    polygon = row.geometry
    for other_idx, other_row in gdf.iterrows():
        if idx < other_idx:
            other_polygon = other_row.geometry
            #判斷是否重疊
            if polygon.overlaps(other_polygon):
                overlap_polygon1 = polygon.intersection(other_polygon)
                overlap_polygon2 = other_polygon.intersection(polygon)
                if overlap_polygon1.geom_type == 'GeometryCollection':
                    for geom in list(overlap_polygon1.geoms):
                        if isinstance(geom, Polygon):
                            exterior1_z = 0
                            for coord in geom.exterior.coords:
                                exterior1_z += round_off(coord[2])
                            exterior1_avg_z = exterior1_z / len(geom.exterior.coords)
                            print("ex1 = ", exterior1_avg_z)
                    for geom in list(overlap_polygon2.geoms):
                        if isinstance(geom, Polygon):
                            exterior2_z = 0
                            for coord in geom.exterior.coords:
                                exterior2_z += round_off(coord[2])
                            exterior2_avg_z = exterior2_z / len(geom.exterior.coords)
                            print("ex2 = ", exterior2_avg_z)
                    difference = abs(exterior1_avg_z-exterior2_avg_z)
                    print("dif = ", difference)
                            # exterior_coords = [(round_off(coord[0]), round_off(coord[1]), round_off(coord[2])) 
                            # for coord in geom.exterior.coords]
                            # print("多邊形的外部邊界座標", exterior_coords)
                else:
                    difference = abs(avg_z(overlap_polygon1)-avg_z(overlap_polygon2))
                    # print(difference)
                
                #寫入gdf新欄位
                if np.isnan(gdf.at[idx, 'Diff']) or difference > gdf.at[idx, 'Diff']:
                    gdf.at[idx, 'Diff'] = difference
                if np.isnan(gdf.at[other_idx, 'Diff']) or difference > gdf.at[other_idx, 'Diff']:
                    gdf.at[other_idx, 'Diff'] = difference

output_shapefile = ouput_file_path + '\\Updated_Plygns3.shp'
gdf.to_file(filename=output_shapefile, driver='ESRI Shapefile', encoding='utf-8')



