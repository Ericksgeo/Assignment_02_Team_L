from shapely.geometry import Point, Polygon
import shapefile
from rtree import index
import numpy as np


nodes = shapefile.Reader(r'Material\roads\nodes.shp')
feature = nodes.shapeRecords()[0]
shapes = nodes.shapes()
points = []
for i in range(len(shapes)):
    points.append(tuple(shapes[i].points[0]))

idx = index.Index()
for n, point in enumerate(points):
    idx.insert(n, point, str(n))

q = (460000,85000)

for j in idx.nearest(q):
    print(j)




