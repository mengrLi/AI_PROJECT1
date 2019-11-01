import math


class Point:
    """
    This is a Point Class, storing the coordination(x,y) of the point.
    """
    x = None
    y = None

    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y


class Node:
    """
    This is a Node Class, storing F, G and H value of each point, it's parent node, and adjacent point list.
    point: a coordination of a point
    father: a parent node
    g: the cost of start point to current point
    h: the evaluate cost of current point to end point
    f: the total of g and h
    weight_list: a list sorting all adjacent point and it's g value; Structure: [((x,y),g), ((x`,y`),g`, .....)]
    """
    point = None
    father = None
    g = 0
    h = 0
    f = 0
    s = 0
    weight_list = []

    def __init__(self, point, g=0, h=0, s=0):
        self.point = point
        self.father = None
        self.g = g
        self.h = h
        self.s = s

    # Calculating H value consider as distance to enf point
    def distance(self, endNode):
        t1 = abs(endNode.point.x - self.point.x) + abs(endNode.point.y - self.point.y)
        t2 = math.sqrt((endNode.point.x - self.point.x) ** 2 + (endNode.point.y - self.point.y) ** 2)
        tmp = (t1 + t2) / 2
        self.h = round(tmp, 4)

    # Calculating F value, the formula is F = G + H.
    def updateF(self):
        tmp = self.g + self.h + self.s
        self.f = round(tmp, 4)

    def safety_parameter(self, current_node, block_list):
        end = (self.point.x, self.point.y)
        weight_list = current_node.weight_list
        tmp = 0
        for item in weight_list:
            if item[0] == end:
                if item[1] == 1.3:
                    tmp = 2
                    break

        if self.s == 0:
            for item in block_list:
                if item[0] == end or item[1] == end:
                    tmp = 1
                    break
        else:
            tmp = 0
        self.s = self.s + tmp
