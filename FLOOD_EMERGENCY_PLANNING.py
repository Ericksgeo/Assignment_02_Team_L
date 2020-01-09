import requests
import tkinter as tk
from tkinter import messagebox
import json
import winsound
import rasterio.mask
from matplotlib.font_manager import FontProperties
from shapely.geometry import Point, Polygon, LineString
import shapefile
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rtree import index
import networkx as nx
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


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

        # use set() to record all elevation data and sort them by numerical value
        ele_set = set()
        [row, col] = area_ele.shape
        for i in range(row):
            for j in range(col):
                ele_set.add(area_ele[i, j])
        ele_value_list = list(ele_set)
        ele_value_list.sort()
        return safe_area, area_ele, ele_value_list


class FindPoint:

    def __init__(self, elevation, area_ele):
        self.__elevation = elevation
        self.__area_ele = area_ele

    def get_highest_point(self, user_point):
        highest_index = np.where(self.__area_ele == np.max(self.__area_ele))
        row_user, col_user = self.__elevation.index(user_point.x, user_point.y)

        # choose the nearest highest point
        nearest_highest_point_index = (highest_index[0][0], highest_index[1][0])
        nearest_highest_point = [[(self.__elevation.xy(highest_index[0][0], highest_index[1][0]))]]
        distance_max = 10000
        if len(highest_index[0]) > 1:
            for i in range(len(highest_index[0])):
                distance = (abs(highest_index[0][i] - row_user) ** 2 + abs(highest_index[1][i] - col_user) ** 2) ** 0.5
                if distance < distance_max:
                    distance_max = distance
                    nearest_highest_point_index = (highest_index[0][i], highest_index[1][i])
                    nearest_highest_point = [[(self.__elevation.xy(highest_index[0][i], highest_index[1][i]))]]
        return nearest_highest_point_index, nearest_highest_point

    def get_median_point(self, ele_value_list):
        median_point_list = []
        median_point_list_a = []
        median_point_list_b = []

        # there are two median points when the length of ele_value_list is even
        if len(ele_value_list) % 2 == 0:
            median_num_a = len(ele_value_list) / 2
            median_num_b = len(ele_value_list) / 2 + 1
            ele_value_a = ele_value_list[int(median_num_a - 1)]
            ele_value_b = ele_value_list[int(median_num_b - 1)]

            # find index of all points with a certain height
            median_index_a = np.where(self.__area_ele == ele_value_a)
            median_index_b = np.where(self.__area_ele == ele_value_b)

            # use rasterio spatial index function to query coordinates of each point
            for i in range(len(median_index_a[0])):
                median_point_list_a.append((self.__elevation.xy(median_index_a[0][i], median_index_a[1][i])))
            for j in range(len(median_index_b[0])):
                median_point_list_b.append((self.__elevation.xy(median_index_b[0][j], median_index_b[1][j])))
            median_point_list = [median_point_list_a, median_point_list_b]

            # calculate the length of ele_value_list need to be iterated in next loop
            list_length_next = len(ele_value_list) / 2

        # only one median point when the length of ele_value_list is odd
        else:
            median_num = len(ele_value_list) // 2 + 1
            ele_value = ele_value_list[int(median_num - 1)]
            median_index = np.where(self.__area_ele == ele_value)
            for i in range(len(median_index[0])):
                median_point_list.append((self.__elevation.xy(median_index[0][i], median_index[1][i])))
            median_point_list = [median_point_list]
            list_length_next = len(ele_value_list) // 2

        return median_point_list, list_length_next


class NearestITN:
    def __init__(self, user_point, target_point, nodes_table):
        self.__user_point = (user_point.x, user_point.y)
        self.__target_point = target_point
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
        for k in idx.nearest(self.__target_point[0][0]):
            index2 = k
        start_point_name = nodes_name_list[index1]
        end_point_name = nodes_name_list[index2]
        return start_point_name, end_point_name, idx


class ShortestPath:

    def __init__(self, elevation, area_ele, user_point, start_point_id, end_point_id, roads_table, nodes_table):
        self.__elevation = elevation
        self.__area_ele = area_ele
        self.__user_point = user_point
        self.__start_point_id = start_point_id
        self.__end_point_id = end_point_id
        self.__elematrix = self.__elevation.read(1)
        self.__roads_table = roads_table
        self.__nodes_table = nodes_table

    def get_path_point(self, g):
        # get the segment points of path, put them into geopandas DataFrame
        links = []
        geom = []
        path_length = 0
        first_node = self.__path[0]
        for node in self.__path[1:]:
            link_fid = g.edges[first_node, node]["fid"]
            links.append(link_fid)
            geom.append(LineString(self.__roads_table[link_fid][1]))
            path_length = path_length + self.__roads_table[link_fid][0]
            first_node = node
        path_gpd = gpd.GeoDataFrame({"fid": links, "geometry": geom})

        # put them into a list
        path_point_list = []
        for i in range(len(geom)):
            for k in range(len(geom[i].xy[0])):
                path_point_x = geom[i].xy[0][k]
                path_point_y = geom[i].xy[1][k]
                path_point_list.append(Point(path_point_x, path_point_y))
        return path_point_list, path_gpd, path_length

    def get_path(self):
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

        # calculate shortest path and determine whether it is a special case
        self.__path = nx.dijkstra_path(g, source=self.__start_point_id, target=self.__end_point_id)
        trigger = 0
        path_point_list, path_gpd, path_length = self.get_path_point(g)
        for points in path_point_list:
            if points.distance(self.__user_point) > 5000:
                trigger = 1
        return self.__path, path_gpd, path_length, g, trigger

    def special_case(self, idx, g, ele_value_list):
        # adding the shortest_path_gpd to return
        nodes_name_list = self.__nodes_table.columns.tolist()
        ele_work_point_list = []
        end_point_id_list = []

        # adding the shortest_path_gpd to return
        while True:
            trigger1 = 0
            length = 0
            ele_value_list_new = ele_value_list
            median_point_list, list_length_next = FindPoint(self.__elevation, self.__area_ele).get_median_point(
                ele_value_list)

            # get list length that contains all points need to be tested
            if len(median_point_list) == 2:
                length = len(median_point_list[0]) + len(median_point_list[1])
            elif len(median_point_list) == 1:
                length = len(median_point_list[0])
            count = 0

            # use for loop to test these points
            for m in range(len(median_point_list)):
                for n in range(len(median_point_list[m])):
                    count += 1

                    # get the ITN and path of each point
                    for x in idx.nearest(median_point_list[m][n]):
                        trigger = 0
                        index2 = x
                        end_point_id_new = nodes_name_list[index2]
                        self.__path = nx.dijkstra_path(g, source=self.__start_point_id, target=end_point_id_new)
                        path_point_list, path_gpd, path_length = self.get_path_point(g)

                        # if any segment point of current path is out of the 5km buffer range, break the loop
                        for points in path_point_list:
                            if points.distance(self.__user_point) > 5000:
                                trigger = 1
                                break

                        # if all points of current path is in the 5km buffer range, trigger will not be activated
                        if trigger == 0:
                            trigger1 = 1
                            ele_work_point_list.append(median_point_list[m][n])
                            end_point_id_list.append(end_point_id_new)

                            # clip the list that contains all the elevation values to search from higher part
                            if len(median_point_list) == 2:
                                ele_value_list_new = ele_value_list[int(list_length_next):]
                            elif len(median_point_list) == 1:
                                if int(list_length_next) == 0:
                                    ele_value_list_new = ele_value_list[0]
                                else:
                                    ele_value_list_new = ele_value_list[(int(list_length_next) + 1):]

                        # if no point of current elevation value pass the test, clip the list to search from lower part
                        if count == length and trigger1 != 1:
                            ele_value_list_new = ele_value_list[0:int(list_length_next)]
                            break

            # break the while loop when finish all tests
            ele_value_list = ele_value_list_new
            if int(list_length_next) == 0:
                break

        # choose the highest one from all candidate points and get its path
        height_max = 0
        index = 0
        for i in range(len(ele_work_point_list)):
            candidate_point = Point(ele_work_point_list[i][0], ele_work_point_list[i][1])
            row, col = self.__elevation.index(candidate_point.x, candidate_point.y)
            height = self.__area_ele[row, col]
            if height >= height_max:
                height_max = height
                index = i
        nearest_highest_point = [[ele_work_point_list[index]]]
        end_point_id = end_point_id_list[index]
        self.__path = nx.dijkstra_path(g, source=self.__start_point_id, target=end_point_id)
        path_point_list, path_gpd, path_length = self.get_path_point(g)

        return self.__path, path_gpd, path_length, nearest_highest_point, end_point_id


class MapPlotting:

    def __init__(self):
        self.fig = plt.figure(figsize=(4, 3), dpi=300)
        self.ax = self.fig.add_subplot(1, 1, 1, projection=ccrs.OSGB())

    def show_result(self, bg_file, elevation, ele_box, path, user_point, start_point, end_point):

        # plot background
        x0 = user_point.x
        y0 = user_point.y
        background = rasterio.open(str(bg_file))
        back_array = background.read(1)
        palette = np.array([value for key, value in background.colormap(1).items()])
        background_image = palette[back_array]
        bounds = background.bounds
        extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
        display_extent = [x0 - 5000, x0 + 5000, y0 - 5000, y0 + 5000]
        if x0 - 5000 < bounds.left:
            display_extent[0] = bounds.left
            display_extent[1] = bounds.left + 10000
        if x0 + 5000 > bounds.right:
            display_extent[1] = bounds.right
            display_extent[0] = bounds.right - 10000
        if y0 - 5000 < bounds.bottom:
            display_extent[2] = bounds.bottom
            display_extent[3] = bounds.bottom + 10000
        if y0 + 5000 > bounds.top:
            display_extent[3] = bounds.top
            display_extent[2] = bounds.top - 10000
        self.ax.imshow(background_image, origin="upper", extent=extent, zorder=1)
        self.ax.set_extent(display_extent, crs=ccrs.OSGB())

        # plot buffer
        buffer = plt.Circle((x0, y0), 5000, color="b", alpha=0.1, zorder=3)
        self.ax.add_patch(buffer)

        # plot point
        start, = plt.plot([start_point[0]], [start_point[1]], "o", color="#F54B0C", markersize=4, zorder=5)
        end, = plt.plot([end_point[0]], [end_point[1]], "o", color="green", markersize=4, zorder=5)

        # plot path
        path.plot(ax=self.ax, edgecolor="#4a59ff", linewidth=1, zorder=4)

        # plot elevation and colorbar
        elevation[elevation == 0] = np.NaN
        left, bottom, right, top = ele_box
        ele = self.ax.imshow(elevation, interpolation='nearest', extent=[left, right, bottom, top], origin="upper",
                             cmap='terrain', zorder=2, alpha=0.4)
        elebar = plt.colorbar(ele, fraction=0.046, pad=0.04)
        elebar.ax.tick_params(labelsize=5)

        # plot scale bar
        # Source: https://pypi.org/project/matplotlib-scalebar/
        font0 = FontProperties()
        font0.set_size(5.)
        font0.set_weight("normal")
        scalebar = ScaleBar(5, length_fraction=0.2, frameon=True, location="lower left")
        scalebar.set_font_properties(font0)
        scalebar.set_box_alpha(0.7)
        self.ax.add_artist(scalebar)

        # plot north arrow
        # Source: https://stackoverflow.com/questions/58088841/how-to-add-a-north-arrow-on-a-geopandas-map
        self.ax.annotate(' ', xy=(0.07, 0.9), xytext=(0.07, 0.78),
                         arrowprops=dict(arrowstyle="wedge", facecolor="black"), ha="center", va="center", fontsize=15,
                         xycoords=self.ax.transAxes)
        self.ax.text(0.07, 0.93, "N", horizontalalignment='center', verticalalignment='center', fontsize=6,
                     transform=self.ax.transAxes)

        # plot legend
        # Source: https://matplotlib.org/3.1.1/tutorials/intermediate/legend_guide.html
        blue_patch = mpatches.Patch(color="b", alpha=0.1, label="5km Area")
        blue_line = mlines.Line2D([], [], linewidth=1, color="#4a59ff", markersize=8, label="Shortest Path")
        plt.legend([blue_patch, blue_line, start, end],
                   ["5km Area", "Shortest Path", "Path Start Point", "Path End Point"],
                   loc="upper right", fontsize=5)

        # plot title
        plt.title("Emergency Path Planning", fontsize=8)

        # show
        plt.show()


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
    safe_area, area_ele, ele_value_list = ClipElevation("Material/elevation/SZ.asc", user_point, island).get_area_ele()
    ele_box = elevation.bounds

    print("Searching for highest point within 5km...")
    nearest_highest_point_index, nearest_highest_point = FindPoint(elevation, area_ele).get_highest_point(user_point)

    # find nearest point
    print("Searching for available path...")
    roads_table, nodes_table = ReadITN("Material/itn/solent_itn.json").get_itn()
    start_point_id, end_point_id, idx = NearestITN(user_point, nearest_highest_point, nodes_table).get_nearest_itn()

    print("Calculating the shortest path...")
    path, path_gpd, path_length, g, trigger = ShortestPath(elevation, area_ele, user_point, start_point_id, end_point_id,
                                              roads_table, nodes_table).get_path()
    if trigger == 1:
        path, path_gpd, path_length, nearest_highest_point, end_point_id = ShortestPath(elevation, area_ele, user_point,
                                                                           start_point_id, end_point_id, roads_table,
                                                                           nodes_table).special_case(idx, g,
                                                                                                     ele_value_list)
    start_point = nodes_table[start_point_id][0]
    end_point = nodes_table[end_point_id][0]
    print("-" * 150)
    print("The highest point available is:", nearest_highest_point[0][0])
    print("Path Start Point ID:", start_point_id, "\nPath End Point ID:", end_point_id)
    print("The shortest distance is : %.2f m" % path_length)
    print("Path created.")

    print("-" * 150)
    print("Map plotting...")
    MapPlotting().show_result("Material/background/raster-50k_2724246.tif",
                              area_ele, ele_box, path_gpd, user_point, start_point, end_point)
    print("Finished.")


if __name__ == "__main__":
    main()
