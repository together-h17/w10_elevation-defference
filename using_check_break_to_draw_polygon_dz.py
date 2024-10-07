# -*- coding: utf-8 -*-

import geopandas as gpd
import glob
from shapely.strtree import STRtree
import os
import pandas as pd
from tqdm import tqdm
from shapely.geometry import Point, Polygon
from decimal import Decimal, ROUND_HALF_UP

polygon_folder = r'C:\Users\AE133\Desktop\三維道路\test'
frame_fn = r'R:\13018\113年_分配範圍.shp'
ouput_file_path = r'C:\Users\AE133\Desktop\三維道路\test\output'

city_code = ['U', 'G']

files = glob.glob(polygon_folder + '\\*.shp')

frames = gpd.read_file(frame_fn, encoding = 'utf-8')
frames_tree = STRtree(frames['geometry'])
# frames_dict = dict((id(plygn), i) for i, plygn in enumerate(frames['geometry']))

float_acc = 4 # If there are round off problem, using "1"
round_off = lambda num : float(Decimal(str(num)).quantize(Decimal('.{}'.format('0' * float_acc)), ROUND_HALF_UP))

progress = tqdm(total = len(files), position = 0, leave = True)
for file in files:
    progress.update(1)
    polygon = gpd.read_file(file, encoding = 'utf-8')
    
    FRAMEID = polygon['FRAMEID'][0]
    # Here take the first value, need to be careful
    assert FRAMEID[0] == '9'
    
    frame_geo = frames[frames['MapID'] == FRAMEID]['geometry'].values[0]
    frames_i = frames_tree.query(frame_geo.buffer(5))
    
    frames_name = []
    for i in frames_i:
        for city in city_code:
            file_path = polygon_folder + '\\' + \
                        'AREA_RD_{0}_{1}.shp'.format(frames.at[i, 'MapID'],
                                                     city)
            frames_name.append(file_path)
            
    checked_fn = []
    for fn in frames_name:
        if os.path.exists(fn):
            checked_fn.append(fn)
    
    plygns = gpd.GeoDataFrame()
    for fn in checked_fn:
        plygn = gpd.read_file(fn, encoding = 'utf-8')
        # plygn = plygn.drop(columns = ['index'])
        plygns = pd.concat([plygns, plygn], axis = 0)
    plygns = plygns.reset_index(drop=True)
    
    plygns_tree = STRtree(plygns['geometry'])
    
    # 初始化新的欄位 'dz' 和 'status'
    plygns['status'] = None
    plygns['dz'] = None
    
    
    for index, row in polygon.iterrows():
        layer = row['LAYER']
        if row['geometry'].geom_type != 'Polygon':
            exterior = row['geometry'].geoms[0].exterior.coords[:]
            interiors = []
            for i in range(1, len(row['geometry'].geoms)):
                interiors.append(row['geometry'].geoms[i].exterior.coords[:])
            
            row['geometry'] = Polygon(shell=exterior, holes=interiors)
        
        for xyz in row['geometry'].exterior.coords:
            query_plygns_i = plygns_tree.query(Point(xyz).buffer(0.1))
            intersect_plygns_i = []
            for plygn_i in query_plygns_i:
                if plygns['geometry'][plygn_i].touches(Point(xyz)):
                    intersect_plygns_i.append(plygn_i)
            
            if len(intersect_plygns_i) > 1:
                for plygn_i in intersect_plygns_i:
                    plygn = plygns['geometry'][plygn_i]
                    queried_plygn = [(round_off(coord[0]), round_off(coord[1]), round_off(coord[2])) 
                                     for coord in plygn.exterior.coords[:]]
                    xyz = (round_off(xyz[0]), round_off(xyz[1]), round_off(xyz[2]))
                    if xyz not in queried_plygn:
                        dz = None
                        status = 'xy different'
                        coords_2D = [(round_off(coord[0]), round_off(coord[1])) 
                                     for coord in plygn.exterior.coords[:]]
                        closest_point = min(coords_2D, key=lambda i: Point(xyz).distance(Point(i)))
                        if (xyz[0], xyz[1]) in coords_2D:
                            status = 'only z different'
                            dzs=[]
                            for coord in plygn.exterior.coords[:]:
                                coord = (round_off(coord[0]), round_off(coord[1]), round_off(coord[2]))
                                if coord[0] == xyz[0] and coord[1] == xyz[1] and coord[2] != xyz[2]:
                                   dzs.append(abs(coord[2] - xyz[2]))
                            dz = sum(dzs) / len(dzs)
                                   
    
                        # 更新 plygns 中的 dz 和 status 欄位
                        plygns.at[plygn_i, 'status'] = status
                        if pd.isna(plygns.at[plygn_i, 'dz']) or dz > plygns.at[plygn_i, 'dz']:
                            plygns.at[plygn_i, 'dz'] = dz
                        

    output_shapefile = ouput_file_path + '\\Updated_Plygns2.shp'
    plygns.to_file(filename=output_shapefile, driver='ESRI Shapefile', encoding='utf-8')

progress.close()        


