# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 11:59:48 2024

@author: vicky
"""

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

# 讀取shapefiles
shapefile = r"C:\Users\AE133\Desktop\三維道路\test\AREA_RD_97221071_G.shp"
#測試資料 A[0, 0, 2, 2] ,  B[6, 5, 4]
gdf = gpd.read_file(shapefile)

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


# 迭代每個polygon
for idx, row in gdf.iterrows():
    polygon = row.geometry
    for other_idx, other_row in gdf.iterrows():
        if idx < other_idx:
            other_polygon = other_row.geometry
            #判斷是否重疊
            if polygon.overlaps(other_polygon):
                # print(avg_z(polygon))
                # print(avg_z(other_polygon))
                difference = abs(avg_z(polygon)-avg_z(other_polygon))
                print(difference)
                #寫入gdf新欄位
                gdf.at[idx, 'Diff'] = difference
                gdf.at[other_idx, 'Diff'] = difference

gdf.to_file('new_AREA_RD_97221071_G.shp')



