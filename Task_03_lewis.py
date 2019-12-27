from shapely.geometry import Point, Polygon
import shapefile
from rtree import index
import numpy as np

class ReadRoadsFromShp:
    def __init__(self, roadsfilepath):
        self.__filepath = roadsfilepath

    def get_nodes(self):
        nodes = shapefile.Reader(self.__filepath)
        shapes = nodes.shapes()
        points = []
        for i in range(len(shapes)):
            points.append(tuple(shapes[i].points[0]))
        return points

class NearestITN:
    def __init__(self, user_point, safe_point, road_point):
        self.__user_point = user_point
        self.__safe_point = safe_point
        self.__road_point = road_point

    def Rtree(self):
        idx = index.Index()
        for n, point in enumerate(self.__road_point):
            idx.insert(n, point, str(n))

        for j in idx.nearest(self.__user_point):
            index1 = j

        for k in idx.nearest(self.__safe_point):
            index2 = k

        return index1, index2

def main():
    highest_point = (448800, 93730)
    user_point = (450000, 85000)
    road_point = ReadRoadsFromShp(r'Material\roads\nodes.shp').get_nodes()
    index1, index2 = NearestITN(user_point, highest_point, road_point).Rtree()
    print(index1, index2)

if __name__ == "__main__":
    main()

