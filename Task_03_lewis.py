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

    def get_point_attribute(self):
        nodes = shapefile.Reader(self.__filepath)
        attribute = nodes.records()
        return attribute


class NearestITN:
    def __init__(self, user_point, safe_point, road_point):
        self.__user_point = user_point
        self.__safe_point = safe_point
        self.__road_point = road_point

    def Rtree(self, attribute):
        idx = index.Index()
        for n, point in enumerate(self.__road_point):
            idx.insert(n, point, str(n))

        for j in idx.nearest(self.__user_point):
            index1 = j

        for k in idx.nearest(self.__safe_point):
            index2 = k

        user_point_name = attribute[index1][0]
        safe_point_name = attribute[index2][0]

        return user_point_name, safe_point_name

def main():
    highest_point = (448800, 93730)
    user_point = (450000, 85000)
    road_point = ReadRoadsFromShp(r'Material\roads\nodes.shp').get_nodes()
    attribute = ReadRoadsFromShp(r'Material\roads\nodes.shp').get_point_attribute()

    user_point_name, safe_point_name = NearestITN(user_point, highest_point, road_point).Rtree(attribute)
    print(user_point_name, safe_point_name)


if __name__ == "__main__":
    main()

