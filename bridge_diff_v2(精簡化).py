# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 10:05:24 2024

@author: AE133
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from decimal import Decimal, ROUND_HALF_UP

# 讀取shapefiles
shapefile = r"D:\Users\ae133\Desktop\三維道路\test\AREA_RD_97221071_G.shp"
ouput_file_path = r"D:\Users\ae133\Desktop\三維道路\test\output"
#測試資料 A[0, 0, 2, 2] ,  B[6, 5, 4]
gdf = gpd.read_file(shapefile)

float_acc = 4 # If there are round off problem, using "1"
round_off = lambda num : float(Decimal(str(num)).quantize(Decimal('.{}'.format('0' * float_acc)), ROUND_HALF_UP))

#計算平均Z值
def avg_z(polygon):
    #解決'GeometryCollection'問題
    if polygon.geom_type == 'GeometryCollection':
        for geom in list(polygon.geoms):
            if isinstance(geom, Polygon):
                exterior_z = sum((coord[2])for coord in geom.exterior.coords)
                exterior_avg_z = exterior_z / len(geom.exterior.coords)
                return exterior_avg_z
    else:
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
                print(idx, other_idx)
                overlap_polygon1 = polygon.intersection(other_polygon)
                overlap_polygon2 = other_polygon.intersection(polygon)
                
                difference = abs(avg_z(overlap_polygon1)-avg_z(overlap_polygon2))
                print("po1 = ", avg_z(overlap_polygon1))
                print("po2 = ", avg_z(overlap_polygon2))
                print(difference)
                print()
                
                #寫入gdf新欄位
                if np.isnan(gdf.loc[idx, 'Diff']) or difference > gdf.loc[idx, 'Diff']:
                    gdf.loc[idx, 'Diff'] = difference
                if np.isnan(gdf.at[other_idx, 'Diff']) or difference > gdf.loc[other_idx, 'Diff']:
                    gdf.loc[other_idx, 'Diff'] = difference


output_shapefile = ouput_file_path + '\\Updated_Plygns4.shp'
gdf.to_file(filename=output_shapefile, driver='ESRI Shapefile', encoding='utf-8')