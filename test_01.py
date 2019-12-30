import requests
import tkinter as tk
from tkinter import messagebox
import json
import winsound
import rasterio.mask
from shapely.geometry import Point, Polygon, LineString, MultiLineString
import shapefile
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rtree import index
import networkx as nx


# creating UI to input coordinates
# Source: Tkinter tutorial https://www.tutorialsteacher.com/python/create-ui-using-tkinter-in-python

class MyWindow:

    def __init__(self):
        # creating labels and positions
        self.root = tk.Tk()
        self.lbl0 = tk.Label(self.root,
                             text='Please insert coordinates (OSGB36) or postal code from the Isle of Wight:')
        self.lbl1 = tk.Label(self.root, text='EAST or first part of postcode (PO30):')
        self.lbl2 = tk.Label(self.root, text='NORTH or last part of postcode (1QB):')
        self.lbl3 = tk.Label(self.root, text='Status:')
        self.lbl4 = tk.Label(self.root, text='')
        self.t1 = tk.Entry(bd=3)
        self.t2 = tk.Entry(bd=3)
        # creating buttons to calculate and reset input
        self.btn1 = tk.Button(self.root, text='plot')
        self.btn2 = tk.Button(self.root, text='clear')
        self.lbl0.place(x=10, y=10)
        self.lbl1.place(x=10, y=50)
        self.t1.place(x=250, y=50)
        self.lbl2.place(x=10, y=100)
        self.t2.place(x=250, y=100)
        self.b1 = tk.Button(self.root, text='calculate', command=self.calculate)
        self.b2 = tk.Button(self.root, text='reset')
        self.b2.bind('<Button-1>', self.reset)
        self.b1.place(x=100, y=150)
        self.b2.place(x=250, y=150)
        self.lbl3.place(x=100, y=200)
        self.lbl4.place(x=150, y=200)
        self.e = ""
        self.n = ""

    def calculate(self):
        # creating "calculate" method to return E and N
        es = self.t1.get()
        nt = self.t2.get()
        if es != "" and nt != "":
            while True:
                try:
                    x = float(es)
                    y = float(nt)
                except ValueError:
                    # Postcode detector sample: PO30 1QB
                    # SOURCE: https: // api.postcodes.io /
                    self.lbl4.configure(text='Characters in the Input, Checking if Postcode match..')
                    postcode = str(es + nt).replace(" ", "")
                    resp = requests.get('https://api.postcodes.io/postcodes/' +
                                        str(postcode))
                    if resp.status_code != 404:
                        self.lbl4.configure(text='It is a postcode')
                        json_data = json.loads(resp.text)
                        qq = json_data["result"]["eastings"]
                        ww = json_data["result"]["northings"]
                        long = json_data["result"]["longitude"]
                        lat = json_data["result"]["latitude"]
                        self.lbl4.configure(text=str("postcode: \n  E: " + str(qq) +
                                                     "\n  N: " + str(ww) + "\n  Lat: " + str(lat) +
                                                     "\n  Long: " + str(long)))
                        self.e = float(qq)
                        self.n = float(ww)
                        break
                    else:
                        self.lbl4.configure(text="No postcode match, try again.")
                        return
                else:
                    # check if coords are in range
                    if (425000 < x < 470000) and (75000 < y < 100000):
                        self.lbl4.configure(text="Initialising, please wait...")
                        self.e = x
                        self.n = y
                        break
                    else:
                        self.lbl4.configure(
                            text="Error: Please Enter Coordinates in range\n (425000,75000)-(470000,"
                                 "100000)\n or a postcode of the Isle of Wight")
                        return
        elif es == "" or nt == "":
            self.lbl4.configure(text="Incomplete Coordinates,\n try again...(Press Reset to clear)")
            return
        self.root.destroy()
        return

    def reset(self, event):
        # creating "reset" method to clear input
        self.t1.delete(0, 'end')
        self.t2.delete(0, 'end')
        self.lbl4.configure(text="")
        winsound.Beep(440, 100)

    def on_closing(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.root.destroy()
            print("Terminated by user.")
            exit()


class ReadIslandFromShp:
    def __init__(self, shpfilepath):
        self.__filepath = shpfilepath

    def get_island_polygon(self):
        # Read shape file
        sf = shapefile.Reader(self.__filepath)
        feature = sf.shapeRecords()[0]
        first = feature.__geo_interface__
        points = first["geometry"]["coordinates"][3][0]
        # create island vector polygon
        island = Polygon(points)
        return island


class ReadElevation:
    def __init__(self, rasterfilepath):
        self.__filepath = rasterfilepath

    def get_elevation_data(self):
        elevation = rasterio.open(self.__filepath)
        return elevation


class ReadITN:
    def __init__(self, itnfilepath):
        self.__filepath = itnfilepath

    def get_itn(self):
        with open(self.__filepath, "r") as f:
            solent_itn = json.load(f)
        roads_table = pd.DataFrame(solent_itn['roadlinks'])
        nodes_table = pd.DataFrame(solent_itn['roadnodes'])
        return roads_table, nodes_table


class PointCheck:
    def __init__(self, polygon, point):
        self.__polygon = polygon
        self.__point = point

    def test_point(self):
        # test whether input point is on the island
        if self.__point.touches(self.__polygon) or self.__point.within(self.__polygon):
            status = "Available Location. Ready for calculation."
            print(status)
        else:
            status = "Error, point outside island boundaries."
            tk.Tk().withdraw()
            messagebox.showerror("ERROR",
                                 "The inserted point is outside Isle of Wight Boundaries, please restart the program...")
            exit()
        return status


class ClipElevation:

    def __init__(self, elefilepath, user_point, island):
        self.__filepath = elefilepath
        self.__user_point = user_point
        self.__island = island

    def get_area_ele(self):
        safe_area = self.__user_point.buffer(5000).intersection(self.__island)
        user_features = [safe_area.__geo_interface__]
        # use shp to clip the rasterio
        # Source: https://www.cnblogs.com/shoufengwei/p/6437934.html
        with rasterio.open(self.__filepath) as scr:
            out_image, out_transform = rasterio.mask.mask(scr, user_features, crop=False)
            area_ele = out_image.reshape(out_image.shape[1], out_image.shape[2])
        return safe_area, area_ele


class HighestPoint:

    def __init__(self, elevation, area_ele):
        self.__elevation = elevation
        self.__area_ele = area_ele

    def get_highest_point(self, user_point):
        highest_index = np.where(self.__area_ele == np.max(self.__area_ele))
        row_user, col_user = self.__elevation.index(user_point.x, user_point.y)

        # choose the nearest highest point
        nearest_highest_point_index = (highest_index[0][0], highest_index[1][0])
        nearest_highest_point = Point(self.__elevation.xy(highest_index[0][0], highest_index[1][0]))
        distance_max = 10000
        if len(highest_index[0]) > 1:
            for i in range(len(highest_index[0])):
                distance = (abs(highest_index[0][i] - row_user) ** 2 + abs(highest_index[1][i] - col_user) ** 2) ** 0.5
                if distance < distance_max:
                    distance_max = distance
                    nearest_highest_point_index = (highest_index[0][i], highest_index[1][i])
                    nearest_highest_point = Point(self.__elevation.xy(highest_index[0][i], highest_index[1][i]))
        return nearest_highest_point_index, nearest_highest_point


class NearestITN:
    def __init__(self, user_point, nearest_highest_point, nodes_table):
        self.__user_point = (user_point.x, user_point.y)
        self.__nearest_highest_point = (nearest_highest_point.x, nearest_highest_point.y)
        self.__nodes_table = nodes_table

    def get_nearest_itn(self):
        idx = index.Index()
        index1 = 0
        index2 = 0
        point_list = []
        nodes_name_list = self.__nodes_table.columns.tolist()
        for node in self.__nodes_table.loc['coords']:
            point_list.append(tuple(node))

        # use Rtree to find nearest point
        for n, point in enumerate(point_list):
            idx.insert(n, point, str(n))
        for j in idx.nearest(self.__user_point):
            index1 = j
        for k in idx.nearest(self.__nearest_highest_point):
            index2 = k
        start_point_name = nodes_name_list[index1]
        end_point_name = nodes_name_list[index2]
        return start_point_name, end_point_name


class ShortestPath:

    def __init__(self, elevation,roads_table, nodes_table):
        self.__elevation = elevation
        self.__elematrix = self.__elevation.read(1)
        self.__roads_table = roads_table
        self.__nodes_table = nodes_table

    def get_path(self, user_point, start_point_name, end_point_name, area_ele, nearest_highest_point_index, nearest_highest_point):

        g = nx.DiGraph()
        for column_name, row in self.__roads_table.iteritems():
            # get coordinates of start point and end point of each road
            point_start_x = row[1][0][0]
            point_start_y = row[1][0][1]
            point_end_x = row[1][-1][0]
            point_end_y = row[1][-1][1]

            # get elevation data of start point and end point of each road
            row_start, col_start = self.__elevation.index(point_start_x, point_start_y)
            row_end, col_end = self.__elevation.index(point_end_x, point_end_y)
            ele_start = self.__elematrix[row_start, col_start]
            ele_end = self.__elematrix[row_end, col_end]

            # Naismith's rule
            distance = row[0]
            time_walk = distance / 5000 * 60
            time_walk_climb = distance / 5000 * 60 + abs(ele_start - ele_end) / 10

            # add weighed edges
            if ele_start >= ele_end:
                g.add_edge(row[2], row[3], fid=column_name, weight=time_walk)
                g.add_edge(row[3], row[2], fid=column_name, weight=time_walk_climb)
            else:
                g.add_edge(row[2], row[3], fid=column_name, weight=time_walk_climb)
                g.add_edge(row[3], row[2], fid=column_name, weight=time_walk)

        # calculate shortest path
        path = nx.dijkstra_path(g, source=start_point_name, target=end_point_name)
        nearest_highest_point = nearest_highest_point

        # adding the shortest_path_gpd to return
        point_list = []
        nodes_name_list = self.__nodes_table.columns.tolist()
        idx = index.Index()
        for node in self.__nodes_table.loc['coords']:
            point_list.append(tuple(node))
            # use Rtree to find nearest point
        for n, point in enumerate(point_list):
            idx.insert(n, point, str(n))

        while True:
            start_point_name_new = start_point_name
            end_point_name_new = end_point_name
            button = 1
            links = []
            geom = []
            first_node = path[0]
            for node in path[1:]:
                link_fid = g.edges[first_node, node]["fid"]
                links.append(link_fid)
                geom.append(LineString(self.__roads_table[link_fid][1]))
                first_node = node
            # path_gpd = gpd.GeoDataFrame({"fid": links, "geometry": geom})

            path_point_list = []
            for i in range(len(geom)):
                for k in range(len(geom[i].xy[0])):
                    path_point_x = geom[i].xy[0][k]
                    path_point_y = geom[i].xy[1][k]
                    path_point_list.append(Point(path_point_x, path_point_y))

            for points in path_point_list:
                if points.distance(user_point) > 5000:
                    print('check1')
                    while True:
                        print(nearest_highest_point_index)
                        area_ele[nearest_highest_point_index[0], nearest_highest_point_index[1]] = 0
                        print(area_ele[nearest_highest_point_index[0], nearest_highest_point_index[1]])
                        nearest_highest_point_index, nearest_highest_point = HighestPoint(self.__elevation,
                                                                                          area_ele).get_highest_point(
                            user_point)
                        print(nearest_highest_point_index)
                        for k in idx.nearest(nearest_highest_point):
                            index2 = k
                        end_point_name_new = nodes_name_list[index2]
                        print(end_point_name_new, end_point_name)
                        print(area_ele[nearest_highest_point_index[0], nearest_highest_point_index[1]])
                        if end_point_name_new != end_point_name:
                            button = 0
                            break
                if button == 0:
                    break
            start_point_name = start_point_name_new
            end_point_name = end_point_name_new
            path = nx.dijkstra_path(g, source=start_point_name_new, target=end_point_name_new)
            print(path)

            if button == 1:
                break

        return path, nearest_highest_point


def main():
    # creating window using the MyWindow class
    window = MyWindow()
    window.root.title('FLOOD EMERGENCY PLANNING')
    window.root.geometry("450x300+500+200")
    window.root.protocol("WM_DELETE_WINDOW", window.on_closing)
    window.root.mainloop()

    # test whether input point is in the island boundaries
    user_point = Point(window.e, window.n)
    island = ReadIslandFromShp("Material/shape/isle_of_wight.shp").get_island_polygon()
    PointCheck(island, user_point).test_point()

    # use island shp to clip the buffer and get the highest point of available area
    print("Loading elevation dataset...")
    elevation = ReadElevation("Material/elevation/SZ.asc").get_elevation_data()
    safe_area, area_ele = ClipElevation("Material/elevation/SZ.asc", user_point, island).get_area_ele()

    print("Searching for highest point within 5km...")
    nearest_highest_point_index, nearest_highest_point = HighestPoint(elevation, area_ele).get_highest_point(user_point)

    # find nearest point
    print("Searching for available path...")
    roads_table, nodes_table = ReadITN("Material/itn/solent_itn.json").get_itn()
    start_point_name, end_point_name = NearestITN(user_point, nearest_highest_point, nodes_table).get_nearest_itn()

    print("Calculating the shortest path...")
    path, nearest_highest_point = ShortestPath(elevation, roads_table, nodes_table).get_path(user_point, start_point_name, end_point_name, area_ele, nearest_highest_point_index, nearest_highest_point)

    print("The highest point available is:", nearest_highest_point)
    print("Path Start Point ID:", start_point_name, "\nPath End Point ID:", end_point_name)
    print(path)
    print("Path created.")


if __name__ == "__main__":
    main()
