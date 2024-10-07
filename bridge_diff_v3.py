# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 11:08:08 2024

@author: AE133
"""
    
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

# 讀取shapefiles
shapefile = r"P:\13018-113年及114年三維道路模型建置案\04_Python\檢核上下層高度\test\AREA_RD_97221071_G.shp"
output_file_path = r"P:\13018-113年及114年三維道路模型建置案\04_Python\檢核上下層高度\output"
gdf = gpd.read_file(shapefile)


def avg_z(overlap_coords):
    return sum(overlap_coords) / len(overlap_coords)

if 'Diff' not in gdf.columns:
    gdf['Diff'] = np.nan


for idx, row in gdf.iterrows():
    polygon = row.geometry
    for other_idx, other_row in gdf.iterrows():
        if idx < other_idx:
            other_polygon = other_row.geometry
            
            # 判斷是否重疊 生成重疊區域polygon
            if polygon.overlaps(other_polygon):
                print(idx, other_idx)
                overlap_coords1, overlap_coords2 = [], []
                overlap_polygon = polygon.intersection(other_polygon)
                
                # 啊我忘記處理 GeometryCollection 的問題
                if polygon.geom_type != 'GeometryCollection' and overlap_polygon.geom_type != 'GeometryCollection':
                    # 抓重疊範圍的z值
                    for po_xyz in polygon.exterior.coords:
                        if any(po_xyz[0] == ov_xy[0] and po_xyz[1] == ov_xy[1] for ov_xy in overlap_polygon.exterior.coords):
                            overlap_coords1.append(po_xyz[2])
                    
                    for po_xyz in other_polygon.exterior.coords:
                        if any(po_xyz[0] == ov_xy[0] and po_xyz[1] == ov_xy[1] for ov_xy in overlap_polygon.exterior.coords):
                            overlap_coords2.append(po_xyz[2])
                    
                    # 計算兩者之間的高度差
                    if overlap_coords1 != [] and overlap_coords2 != []:
                        difference = abs(avg_z(overlap_coords1) - avg_z(overlap_coords2))
                        print("po1 = ", avg_z(overlap_coords1))
                        print("po2 = ", avg_z(overlap_coords2))
                        print(difference)
                        print()
                    
                    # 寫入gdf新欄位
                    if np.isnan(gdf.loc[idx, 'Diff']) or difference > gdf.loc[idx, 'Diff']:
                        gdf.loc[idx, 'Diff'] = difference
                    if np.isnan(gdf.loc[other_idx, 'Diff']) or difference > gdf.loc[other_idx, 'Diff']:
                        gdf.loc[other_idx, 'Diff'] = difference

# 輸出更新後的shapefile
output_shapefile = output_file_path + '\\output.shp'
gdf.to_file(filename=output_shapefile, driver='ESRI Shapefile', encoding='utf-8')
