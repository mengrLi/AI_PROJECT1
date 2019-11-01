import cmd
import sys
import time
import math
import geopandas as gpd
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.patches as patches
from shapely.geometry import Polygon
from descartes import PolygonPatch
import methodOne
import methodTwo
import node_class


class creatGraph(cmd.Cmd):
    """
    This class is mainly about create the figure, including load "crime.shp",
    draw grid map with block and non-block grids by different colors.
    After find the path, show the optimal path on the figure as final result.
    """

    def __init__(self):
        super().__init__()
        self.main()

    # This method is used to create the "grid.shp" file, join with the "crime.shp" to get all points in grids
    # use group by function calculate the number of points in each grid, store in grid.shp file as column 'total'
    # sort all the grids by 'total' ,store as grid_sorted and get the value by threshold in grid_sorted
    # draw the final result graph with optimal path
    def showMap(self, grid_size, threshold):
        index = []
        polygons = []
        point_set = []  # store all the point in grid map
        x_list = []
        sorted_point_list = []
        block_list = []  # block_list is used to store all edge in a block grid

        # load crime_dt,.shp as a point map
        crime = gpd.read_file("Shape/crime_dt.shp")
        plt.show()
        # draw a grid map store as grid.shp, each grid is a Polygon object store in polygons[]
        self.draw_grid(grid_size, x_list, point_set, polygons, index)
        grid = gpd.GeoDataFrame({'geometry': polygons, 'index': index})
        grid.to_file("grid.shp")

        crime_map = gpd.sjoin(crime, grid, op='within')
        counts_point = crime_map.groupby('index').size()
        grid = grid.merge(counts_point.reset_index(name="total"), how='left', on='index')
        grid_sorted = grid.sort_values('total', ascending=False)
        grids = len(grid)
        blocks = int(grids * (1 - float(threshold)))
        value = 0
        for i in grid_sorted.iterrows():
            if blocks == 0:
                value = i[1]['total']
                break
            else:
                blocks -= 1
        x_list.append(-73.55)
        self.creat_sorted_point_list(point_set, x_list, sorted_point_list)

        # set the figure value
        fig = plt.figure()
        plt.xlim((-73.59, -73.55))
        plt.ylim((45.49, 45.53))
        my_x_ticks = np.arange(-73.59, -73.55, 0.005)
        my_y_ticks = np.arange(45.49, 45.53, 0.005)
        plt.xticks(my_x_ticks)
        plt.yticks(my_y_ticks)

        map_size = math.ceil(round((-73.55 - (-73.59)) / float(grid_size), 4))

        if len(grid) > map_size * map_size:
            map_size = map_size + 1

        if value > 0:
            # threshole can get a value
            for i in grid.iterrows():
                total = i[1]['total']
                ax = fig.gca()
                ax.add_patch(PolygonPatch(i[1]['geometry'], fc='Blue'))
                if total >= value:
                    ax = fig.gca()
                    ax.add_patch(PolygonPatch(i[1]['geometry'], fc='yellow'))
                    temp_t = Polygon(i[1]['geometry']).bounds
                    block_list = self.creat_blocklist(block_list, temp_t)
        else:
            # threshold is too small, get the value is nan, than only consider the grid with points as a block grid
            for i in grid.iterrows():
                total = i[1]['total']
                ax = fig.gca()
                ax.add_patch(PolygonPatch(i[1]['geometry'], fc='Blue'))
                if total > 0:
                    ax = fig.gca()
                    ax.add_patch(PolygonPatch(i[1]['geometry'], fc='yellow'))
                    temp_t = Polygon(i[1]['geometry']).bounds
                    block_list = self.creat_blocklist(block_list, temp_t)

        # find the optimal path by A star algorithm

        graph = self.buildMap(sorted_point_list, block_list, map_size + 1)

        start_node = (-73.55, 45.49)
        end_node = (-73.59, 45.53)
        start_time = time.time()
        result_one = methodOne.algorithm_one(graph, graph[start_node], graph[end_node]).a_star()
        # result_one = findAnotherPath.Astar_Algorithm(graph, graph[start_node], graph[end_node],block_list).a_star()
        end_time = time.time()
        self.print_result(grid, result_one)
        print("Method1:Time For Searching Path: " + str(end_time - start_time))
        # print("Method1: COST is :" + str(cost))
        if result_one is not None:
            color = 'red'
            self.creat_resultlist(result_one, color)
        plt.show()

    # This method is used to show the result in console
    def print_result(self, grid, result_path):
        total_crimes = np.sum(grid.get('total'))
        average_crimes = total_crimes / len(grid)
        mean_crimes = np.mean(grid.get('total'))
        standard_deviation = np.std(grid.get('total'))
        print('--------------   Total  of the grid is:' + str(total_crimes) + '----------------')
        print('--------------   Average  of the grid  is:' + str(average_crimes) + '----------------')
        print('--------------   Mean of the grid is:' + str(mean_crimes) + '----------------')
        print('--------------   Standard deviation of the grid is:' + str(standard_deviation) + '----------------')

    # This method is create a block list to store all the edges in the block grids
    # temp_t is the list store four poins in a block grid
    def creat_blocklist(self, block_list, temp_t):
        # low bound,two way
        block_list.append(((temp_t[0], temp_t[1]), (temp_t[2], temp_t[1])))
        block_list.append(((temp_t[2], temp_t[1]), (temp_t[0], temp_t[1])))
        # right
        block_list.append(((temp_t[2], temp_t[1]), (temp_t[2], temp_t[3])))
        block_list.append(((temp_t[2], temp_t[3]), (temp_t[2], temp_t[1])))
        # upper
        block_list.append(((temp_t[2], temp_t[3]), (temp_t[0], temp_t[3])))
        block_list.append(((temp_t[0], temp_t[3]), (temp_t[2], temp_t[3])))
        # left
        block_list.append(((temp_t[0], temp_t[3]), (temp_t[0], temp_t[1])))
        block_list.append(((temp_t[0], temp_t[1]), (temp_t[0], temp_t[3])))
        # left to right
        block_list.append(((temp_t[0], temp_t[1]), (temp_t[2], temp_t[3])))
        block_list.append(((temp_t[2], temp_t[3]), (temp_t[0], temp_t[1])))
        # right to left
        block_list.append(((temp_t[2], temp_t[1]), (temp_t[0], temp_t[3])))
        block_list.append(((temp_t[0], temp_t[3]), (temp_t[2], temp_t[1])))

        return block_list

    # This method is used to draw the optimal path as a red line in graph
    def creat_resultlist(self, result_path, color):
        x_list = []
        y_list = []
        for item in result_path:
            x_list.append(item[0])
            y_list.append(item[1])
        ax = plt.gca()
        ax.plot(x_list, y_list, color=color, linewidth=3)

    # This method is used to create a block edges dictionary,
    # which contains all edges in a block grid and calculate the number of times edges contains in a block list
    # key: points(x,y); value:times
    # value is used to define weight
    def creat_block_dict(self, blocklist):
        my_dict = {}
        for item in blocklist:
            if not item in my_dict:
                my_dict[item] = 1
            else:
                temp = my_dict[item] + 1
                my_dict[item] = temp
        return my_dict

    # This method is used to create a map according to the block list.
    # And the data structure of the map is dictionary, storing each point coordination and node.
    def buildMap(self, sort_point_list, block_list, map_size):
        my_dict = self.creat_block_dict(block_list)
        map = {}
        for i in range(len(sort_point_list)):
            for j in range(len(sort_point_list[i])):

                current_point = sort_point_list[i][j]
                point = node_class.Node(node_class.Point(current_point[0], current_point[1]))

                #      build weight list
                weight_list = []
                # left top node (-1, -1) weight = 1.5
                if 0 <= i - 1 < map_size and 0 <= j - 1 < map_size:
                    tmp_point = sort_point_list[i - 1][j - 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, True)

                    # left node (0, -1) weight = 1 or 1.3
                if 0 <= j - 1 < map_size:
                    tmp_point = sort_point_list[i][j - 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, False)

                    # left down node (+1, -1) weight = 1.5
                if 0 <= i + 1 < map_size and 0 <= j - 1 < map_size:
                    tmp_point = sort_point_list[i + 1][j - 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, True)

                    # mediate top node (-1, 0) weight = 1 or 1.3
                if 0 <= i - 1 < map_size:
                    tmp_point = sort_point_list[i - 1][j]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, False)

                    # mediate down node (+1, 0) weight = 1 or 1.3
                if 0 <= i + 1 < map_size:
                    tmp_point = sort_point_list[i + 1][j]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, False)

                    # right top node(-1, +1) weight = 1.5
                if 0 <= i - 1 < map_size and 0 <= j + 1 < map_size:
                    tmp_point = sort_point_list[i - 1][j + 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, True)

                    # right node (0, +1) weight = 1 or 1.3
                if 0 <= j + 1 < map_size:
                    tmp_point = sort_point_list[i][j + 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, False)

                    # right down node(+1, +1) weight = 1.5
                if 0 <= i + 1 < map_size and 0 <= j + 1 < map_size:
                    tmp_point = sort_point_list[i + 1][j + 1]
                    self.getWeight(weight_list, my_dict, current_point, tmp_point, True)

                point.weight_list = weight_list
                map[current_point] = point

        return map

    # Calculating each edge's g value, like 1, 1.3 or 1.5.
    # And put it into weight list.
    def getWeight(self, weight_list, my_dict, start_point, end_point, flag):
        edge = (start_point, end_point)
        time = 0
        if edge in my_dict:
            time = my_dict[edge]

        if self.boundaryline(start_point, end_point) is not True:
            if flag:
                if time == 0:
                    tuple = (end_point, 1.5)
                    weight_list.append(tuple)
            else:
                if time == 0:
                    tuple = (end_point, 1)
                    weight_list.append(tuple)
                elif time == 1:
                    tuple = (end_point, 1.3)
                    weight_list.append(tuple)

    # Judging the edge is boundary line or not. True is boundary line, otherwise is not.
    def boundaryline(self, start_point, end_point):
        if (start_point[0] == -73.59 and end_point[0] == -73.59) or (
                start_point[0] == -73.55 and end_point[0] == -73.55):
            return True
        elif (start_point[1] == 45.49 and end_point[1] == 45.49) or (start_point[1] == 45.53 and end_point[1] == 45.53):
            return True
        else:
            return False

    # This method is used to draw a grid graph
    # consider the grid size is user defined, which may not be exact divided with no reminder,
    # four cases designed to draw a grid with a reminder size
    def draw_grid(self, grid_size, x_list, point_set, polygons, index):
        x = -73.59
        x_size = 1
        grid_num = 0
        if 0.04 % float(grid_size) == 0:
            size = 0.04 / float(grid_size)
        else:
            size = 0.04 // float(grid_size) + 1

        temp_size = round(0.04 - float(grid_size) * (size - 1), 5)
        while x <= -73.55 and x + temp_size <= -73.54999:
            x_list.append(round(x, 5))
            y = 45.49
            y_size = 0
            while y < 45.5299999:
                y_size = y_size + 1
                if y_size == size and x_size != size:

                    rect = patches.Rectangle(
                        (x, y),
                        # (x,y)
                        float(grid_size),
                        # width
                        float(temp_size),
                        # height
                        edgecolor='b',
                        facecolor="w"
                    )
                    polygons.append(
                        Polygon([(round(x, 5), round(y, 5)), (round(x + float(grid_size), 5), round(y, 5)),
                                 (round(x + float(grid_size), 5), round(y + float(temp_size), 5)),
                                 (round(x, 5), round(y + float(temp_size), 5))]))

                    point_set.append((round(x, 5), round(y, 5)))
                    point_set.append((round(x, 5), round(y + float(temp_size), 5)))
                    point_set.append((round(x + float(grid_size), 5), round(y, 5)))
                    point_set.append((round(x + float(grid_size), 5), round(y + float(temp_size), 5)))

                elif y_size != size and x_size != size:
                    rect = patches.Rectangle(
                        (x, y),
                        # (x,y)
                        float(grid_size),
                        # width
                        float(grid_size),
                        # height
                        edgecolor='b',
                        facecolor="w"
                    )
                    polygons.append(
                        Polygon([(round(x, 5), round(y, 5)), (round(x + float(grid_size), 5), round(y, 5)),
                                 (round(x + float(grid_size), 5), round(y + float(grid_size), 5)),
                                 (round(x, 5), round(y + float(grid_size), 5))]))

                    point_set.append((round(x, 5), round(y, 5)))
                    point_set.append((round(x, 5), round(y + float(grid_size), 5)))
                    point_set.append((round(x + float(grid_size), 5), round(y, 5)))
                    point_set.append((round(x + float(grid_size), 5), round(y + float(grid_size), 5)))
                elif x_size == size and y_size != size:
                    rect = patches.Rectangle(
                        (x, y),
                        # (x,y)
                        float(temp_size),
                        # width
                        float(grid_size),
                        # height
                        edgecolor='b',
                        facecolor="w"
                    )
                    polygons.append(
                        Polygon([(round(x, 5), round(y, 5)), (round(x + float(temp_size), 5), round(y, 5)),
                                 (round(x + float(temp_size), 5), round(y + float(grid_size), 5)),
                                 (round(x, 5), round(y + float(grid_size), 5))]))

                    point_set.append((round(x, 5), round(y, 5)))
                    point_set.append((round(x, 5), round(y + float(grid_size), 5)))
                    point_set.append((round(x + float(temp_size), 5), round(y, 5)))
                    point_set.append((round(x + float(temp_size), 5), round(y + float(grid_size), 5)))

                elif x_size == size and y_size == size:
                    rect = patches.Rectangle(
                        (x, y),
                        # (x,y)
                        float(temp_size),
                        # width
                        float(temp_size),
                        # height
                        edgecolor='b',
                        facecolor="w"
                    )
                    polygons.append(
                        Polygon([(round(x, 5), round(y, 5)), (round(x + float(temp_size), 5), round(y, 5)),
                                 (round(x + float(temp_size), 5), round(y + float(temp_size), 5)),
                                 (round(x, 5), round(y + float(temp_size), 5))]))

                    point_set.append((round(x, 5), round(y, 5)))
                    point_set.append((round(x, 5), round(y + float(temp_size), 5)))
                    point_set.append((round(x + float(temp_size), 5), round(y, 5)))
                    point_set.append((round(x + float(temp_size), 5), round(y + float(temp_size), 5)))

                index.append(grid_num)
                grid_num = grid_num + 1
                y = y + float(grid_size)
            x_size = x_size + 1
            x = x + float(grid_size)

    # This method is used to create a points list sorted by order
    # sorted_point_list is used to build a points map
    def creat_sorted_point_list(self, point_set, x_list, sorted_point_list):
        tmp_sort_list = []
        for item in point_set:
            if not item in tmp_sort_list:
                tmp_sort_list.append(item)

        iter = 0
        while iter < len(x_list):
            list = []
            for item in tmp_sort_list:
                if item[0] == x_list[iter]:
                    list.append(item)
            sorted_point_list.append(list)
            iter += 1

    def showMap2(self, grid_size, threshold):
        index = []
        polygons = []
        point_set = []  # store all the point in grid map
        x_list = []
        sorted_point_list = []
        block_list = []  # block_list is used to store all edge in a block grid

        # load crime_dt,.shp as a point map
        crime = gpd.read_file("Shape/crime_dt.shp")
        # draw a grid map store as grid.shp, each grid is a Polygon object store in polygons[]
        self.draw_grid(grid_size, x_list, point_set, polygons, index)
        grid = gpd.GeoDataFrame({'geometry': polygons, 'index': index})
        grid.to_file("grid.shp")
        crime_map = gpd.sjoin(crime, grid, how='left', op='within')
        counts_point = crime_map.groupby('index').size()
        grid = grid.merge(counts_point.reset_index(name="total"), how='left', on='index')
        grid_sorted = grid.sort_values('total', ascending=False)
        grids = len(grid)
        blocks = int(grids * (1 - float(threshold)))
        value = 0
        for i in grid_sorted.iterrows():
            if blocks == 0:
                value = i[1]['total']
                break
            else:
                blocks -= 1
        x_list.append(-73.55)
        self.creat_sorted_point_list(point_set, x_list, sorted_point_list)

        # set the figure value
        fig = plt.figure()
        plt.xlim((-73.59, -73.55))
        plt.ylim((45.49, 45.53))
        my_x_ticks = np.arange(-73.59, -73.55, 0.005)
        my_y_ticks = np.arange(45.49, 45.53, 0.005)
        plt.xticks(my_x_ticks)
        plt.yticks(my_y_ticks)

        map_size = math.ceil(round((-73.55 - (-73.59)) / float(grid_size), 4))

        if len(grid) > map_size * map_size:
            map_size = map_size + 1

        if value > 0:
            # threshole can get a value
            for i in grid.iterrows():
                total = i[1]['total']
                ax = fig.gca()
                ax.add_patch(PolygonPatch(i[1]['geometry'], fc='Blue'))
                if total >= value:
                    ax = fig.gca()
                    ax.add_patch(PolygonPatch(i[1]['geometry'], fc='yellow'))
                    temp_t = Polygon(i[1]['geometry']).bounds
                    block_list = self.creat_blocklist(block_list, temp_t)
        else:
            # threshold is too small, get the value is nan, than only consider the grid with points as a block grid
            for i in grid.iterrows():
                total = i[1]['total']
                ax = fig.gca()
                ax.add_patch(PolygonPatch(i[1]['geometry'], fc='Blue'))
                if total > 0:
                    ax = fig.gca()
                    ax.add_patch(PolygonPatch(i[1]['geometry'], fc='yellow'))
                    temp_t = Polygon(i[1]['geometry']).bounds
                    block_list = self.creat_blocklist(block_list, temp_t)
        # find the optimal path by A star algorithm

        graph = self.buildMap(sorted_point_list, block_list, map_size + 1)

        start_node = (-73.55, 45.49)
        end_node = (-73.59, 45.53)

        start_time = time.time()
        result_two = methodTwo.algorithm_two(graph, graph[start_node], graph[end_node], block_list).a_star()
        end_time = time.time()
        print("Method2:Time For Searching Path: " + str(end_time - start_time))
        if result_two is not None:
            color = 'white'
            self.creat_resultlist(result_two, color)

        plt.show()

    def main(self):
        while True:
            print("----------please give the grid size(0.000): -----------")
            grid_size = sys.stdin.readline().rstrip('\r\n')
            print("----------please give the threshold(0%-100%): -----------")
            threshold = sys.stdin.readline().rstrip('\r\n')
            self.showMap(grid_size, threshold)
            self.showMap2(grid_size, threshold)
            sys.exit(0)


if __name__ == "__main__":
    t = creatGraph()
    t.main()
