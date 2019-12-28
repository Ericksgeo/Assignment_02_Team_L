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




with open(r'Material\elevation\SZ.asc', "r") as elefile:
    ele_lines = elefile.readlines()

elevation = rasterio.open(r'Material\elevation\SZ.asc')
matrix=elevation.read(1)

print(elevation.width)
cellsize=(elevation.bounds.right-elevation.bounds.left)/9000


coor=elevation.xy(elevation.height // 2, elevation.width // 2)
print(coor)

with open(r'Material\itn\solent_itn.json', "r") as f:
    solent_itn = json.load(f)

# road_links = solent_itn['roadlinks']
# road_nodes = solent_itn['roadnodes']
df=pd.DataFrame(solent_itn['roadlinks'])




# g = nx.DiGraph()
#

#
# for link in road_links.values():
#
#     point_a_id = link['start']
#     point_b_id = link['end']
#     print(point_a_id)


