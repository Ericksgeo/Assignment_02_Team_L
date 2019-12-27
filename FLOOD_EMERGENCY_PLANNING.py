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
from rtree import index


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
                        self.lbl4.configure(text="Error: Please Enter Coordinates in range\n (425000,75000)-(470000,"
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


class ReadRoadsFromShp:
    def __init__(self, roadsfilepath):
        self.__filepath = roadsfilepath

    def get_nodes(self):
        nodes = shapefile.Reader(self.__filepath)
        shapes = nodes.shapes()
        pointslist = []
        for i in range(len(shapes)):
            pointslist.append(tuple(shapes[i].points[0]))
        return pointslist

    def get_point_attribute(self):
        nodes = shapefile.Reader(self.__filepath)
        attribute = nodes.records()
        return attribute


class PointCheck:
    def __init__(self, polygon, point):
        self.__polygon = polygon
        self.__point = point

    def test_point(self):
        # test whether input point is on the island
        if self.__point.touches(self.__polygon) or self.__point.within(self.__polygon):
            status = "Available Location. Calculating..."
            print(status)
        else:
            status = "Error, point outside island boundaries."
            tk.Tk().withdraw()
            messagebox.showerror("ERROR",
                                 "The inserted point is outside Isle of Wight Boundaries, please restart the program...")
            exit()
        return status


class Elevation:

    def __init__(self, elefilepath, point):
        self.__point = point
        self.__filepath = elefilepath
        with open(self.__filepath, "r") as elefile:
            self.__lines = elefile.readlines()

    def get_ncols(self):
        line0 = self.__lines[0].strip().split(" ")
        return int(line0[1])

    def get_nrows(self):
        line1 = self.__lines[1].strip().split(" ")
        return int(line1[1])

    def get_xllcorner(self):
        line2 = self.__lines[2].strip().split(" ")
        return int(line2[1])

    def get_yllcorner(self):
        line3 = self.__lines[3].strip().split(" ")
        return int(line3[1])

    def get_cellsize(self):
        line4 = self.__lines[4].strip().split(" ")
        return int(line4[1])

    def get_elevation(self):
        # init elevation data array
        self.__eledata = np.zeros((self.get_nrows(), self.get_ncols()))
        row = 0
        # read elevation data
        for line in self.__lines[5:]:
            eachline = line.strip().split(" ")
            for col in range(0, self.get_ncols()):
                self.__eledata[row, col] = float(eachline[col])
            row = row + 1
        return self.__eledata

    def get_max_point(self, area):
        features = [area.__geo_interface__]
        # use shp to clip the rasterio
        # Source: https://www.cnblogs.com/shoufengwei/p/6437934.html
        with rasterio.open(self.__filepath) as scr:
            out_image, out_transform = rasterio.mask.mask(scr, features, crop=False)
            area_ele = out_image.reshape(out_image.shape[1], out_image.shape[2])
        a = np.where(area_ele == np.max(area_ele))
        shortest = 10000
        # get the nearst point in Euclidean distance as the only highest point
        for m in range(0, len(a[0])):
            xr = a[0][m]
            yr = a[1][m]
            xo = self.get_cellsize() * yr + self.get_cellsize() / 2 + self.get_xllcorner()
            yo = self.get_nrows() * self.get_cellsize() - (
                    self.get_cellsize() * xr + self.get_cellsize() / 2) + self.get_yllcorner()
            xu = self.__point.x
            yu = self.__point.y
            distance = np.sqrt((xu - xo) ** 2 + (yu - yo) ** 2)
            if distance < shortest:
                shortest = distance
                hpoint = Point(xo, yo)
        return hpoint


class NearestITN:
    def __init__(self, user_point, safe_point, road_point):
        self.__user_point = (user_point.x,user_point.y)
        self.__safe_point = (safe_point.x,safe_point.y)
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
    # creating window using the MyWindow class
    window = MyWindow()
    window.root.title('FLOOD EMERGENCY PLANNING')
    window.root.geometry("450x300+100+100")
    window.root.mainloop()

    # test whether input point is in the island boundaries
    user_point = Point(window.e, window.n)
    island = ReadIslandFromShp("Material/shape/isle_of_wight.shp").get_island_polygon()
    PointCheck(island, user_point).test_point()

    # use island shp to clip the buffer and get the highest point of available area
    shparea = user_point.buffer(5000).intersection(island)
    highest_point = Elevation("Material/elevation/SZ.asc", user_point).get_max_point(shparea)
    print("The highest point within 5km is:", highest_point)

    # find nearest point
    road_point = ReadRoadsFromShp(r'Material\roads\nodes.shp').get_nodes()
    attribute = ReadRoadsFromShp(r'Material\roads\nodes.shp').get_point_attribute()

    user_point_name, safe_point_name = NearestITN(user_point, highest_point, road_point).Rtree(attribute)
    print("Path Start Point ID:",user_point_name, "\nPath End Point ID:",safe_point_name)

if __name__ == "__main__":
    main()
