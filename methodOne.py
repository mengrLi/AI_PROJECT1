
class algorithm_one:

    """
    This is a Astar_Algorithm class. It is aim to achieve A star algorithm.
    map: a dictionary storing all Nodes; Structure: {(x,y): Node, (x`,y`): Node`, ...}
    starNode: star point and it's coordination is (-73.55, 45.49)
    endNode: end point and it;s coordination is (-73.59,45.53)
    path: a list storing the result path
    """

    def __init__(self, map, startNode, endNode):
        self.map = map
        self.startNode = startNode
        self.endNode = endNode
        self.path = []

    # Finding a node in open list with smallest F value.
    def sort_openList_byF(self, open_list):
        min_node = open_list[0]
        for node in open_list:
            if node.f < min_node.f:
                min_node = node
        return min_node

    # Updating node information, like father node, g, h and f value.
    def update_node(self, current_node, adj_node, endNode, g):
        adj_node.father = current_node
        adj_node.g = current_node.g + g
        adj_node.distance(endNode)
        adj_node.updateF()

    # Using searching parent node method to get a path from end point to start point.
    def get_path(self, endNode):
        path = []
        tmp = endNode
        while tmp is not None:
            location = (tmp.point.x, tmp.point.y)
            path.append(location)
            tmp = tmp.father
        return path

    # A star algorithm implementation
    def a_star(self):
        if len(self.startNode.weight_list) == 0 or len(self.endNode.weight_list) == 0:
            print(" Due to blocks, no path is found. Please change the map and try again ")
        else:
            open_list = [self.startNode]
            close_list = []
            target = (self.endNode.point.x, self.endNode.point.y)

            while len(open_list) != 0:
                current_node = self.sort_openList_byF(open_list)    # find a node with smallest F value
                open_list.remove(current_node)
                close_list.append(current_node)

                if len(current_node.weight_list) != 0:
                    for item in current_node.weight_list:
                        if item[0] == target:                       # if the end point coordination in weight_list then break
                            self.endNode.father = current_node
                            self.endNode.g = current_node.g + item[1]
                            self.endNode.distance(self.endNode)
                            self.endNode.updateF()
                            break
                        else:
                            adj_node = self.map[item[0]]
                            g = item[1]
                            if adj_node not in close_list:          # next node cannot in close list
                                if adj_node in open_list:           # if next node in open list, then judge the node whether need to be updated info or not
                                    if adj_node.g > current_node.g + g:
                                        adj_node.g = current_node.g + g
                                        adj_node.updateF()
                                        adj_node.father = current_node
                                else:                               # if next node doesn't in open list and close list, then updating it's info and put it into open list
                                    open_list.append(adj_node)
                                    self.update_node(current_node, adj_node, self.endNode, g)

            cost = 0
            path = self.get_path(self.endNode)
            if self.endNode.father is None:
                print(" Due to blocks, no path is found. Please change the map and try again ")
                return path
            else:
                cost = self.endNode.f
                print("Method1: COST is :" + str(cost))
                return path
