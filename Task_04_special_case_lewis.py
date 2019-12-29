import requests
import tkinter as tk
from tkinter import messagebox
import json
import winsound
import rasterio.mask
from shapely.geometry import Point, Polygon
import shapefile
import numpy as np
import rasterio
from rasterio import plot
import pandas as pd
import networkx as nx


elevation = rasterio.open(r'Material\elevation\SZ.asc')
ele_matrix = elevation.read(1)

cellsize=(elevation.bounds.right-elevation.bounds.left)/9000
coor=elevation.xy(elevation.height // 2, elevation.width // 2)


with open(r'Material\itn\solent_itn.json', "r") as f:
    solent_itn = json.load(f)
road_table = pd.DataFrame(solent_itn['roadlinks'])
print(road_table)

g = nx.DiGraph()
for column_name, row in road_table.iteritems():
    # get coordinates of start point and end point of each road
    point_start_x = row[1][0][0]
    point_start_y = row[1][0][1]
    point_end_x = row[1][-1][0]
    point_end_y = row[1][-1][1]

    # get elevation data of start point and end point of each road
    row_start, col_start = elevation.index(point_start_x, point_start_y)
    row_end, col_end = elevation.index(point_start_x, point_start_y)
    ele_start = ele_matrix[row_start, col_start]
    ele_end = ele_matrix[row_end, col_end]

    # Naismith's rule
    distance = row[0]
    time_walk = distance/5000*60
    time_walk_climb = distance/5000*60+abs(ele_start-ele_end)/10

    # add weighed edges
    if ele_start >= ele_end:
        g.add_edge(row[2], row[3], weight=time_walk)
        g.add_edge(row[3], row[2], weight=time_walk_climb)
    else:
        g.add_edge(row[2], row[3], weight=time_walk_climb)
        g.add_edge(row[3], row[2], weight=time_walk)

# calculate shortest path
path = nx.dijkstra_path(g, source='osgb4000000026219225', target='osgb4000000026240451')
print(path)





#
#

#
# for link in road_links.values():
#
#     point_a_id = link['start']
#     point_b_id = link['end']
#     print(point_a_id)


