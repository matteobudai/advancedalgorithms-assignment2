from cmath import inf
import math
from collections import defaultdict
import time

#Create a node
class Node:
    def __init__(self, tag: int):
        self.tag = tag
        self.key = None
        self.parent = None
        self.isPresent = True
        self.index = tag-1 # Track the index of the node in the heap instead of using list.index() method which is O(n)
        self.adjacencyList = []
        self.xcoord = 0
        self.ycoord = 0
        self.weight_type = None
        
#Transform coordinate data into connected graph (list in format [u, v, weight])
class Graph:
    def __init__(self):
        self.nodes = defaultdict(Node)

    def createNodes(self, nums: int):
        for i in range(1, nums+1): # nums+1 in order to cover the last node
            self.nodes[i] = Node(i)

    def buildGraph(self, input):
        lines = input.readlines()

        #determine what line the coordinate data starts and ends on (some files have extra blank lines at the end)
        #node_data_start = 0
        node_data_end = 0
        while lines[node_data_end].split()[0] != 'EOF':
            if lines[node_data_end].split()[0] == 'NODE_COORD_SECTION':
                node_data_start = node_data_end
            node_data_end += 1

        #determine weight type (always on line 4)
        line_4 = lines[4].split()
        weight_type = line_4[len(line_4)-1]

        #determine dimension (always line 3), calculate number of nodes, call create nodes method 
        line_3 = lines[3].split()
        dimension = int(line_3[len(line_3)-1])
        self.createNodes(dimension)

        connected_graph = []

        #create a connected graph by looping through each node and connecting it to all other nodes
        #calculate the weight between each node depending on the weight_type defined in the graph
        for iterate in range(node_data_start+1, node_data_end):

            node_u = int(lines[iterate].split()[0])
            u_x1_coord = float(lines[iterate].split()[1])
            u_y1_coord = float(lines[iterate].split()[2])

            for interate_v in range(node_data_start+1, node_data_end):
                node_v = int(lines[interate_v].split()[0])
                if node_u == node_v: #do nothing
                    1
                elif node_u != node_v: 
                    v_x2_coord = float(lines[interate_v].split()[1])
                    v_y2_coord = float(lines[interate_v].split()[2])
                    weight_uv = weight(weight_type, u_x1_coord, v_x2_coord, u_y1_coord, v_y2_coord)
                    connected_graph.append([node_u,node_v,weight_uv])
                    self.makeNodes(node_u,u_x1_coord,u_y1_coord,node_v, weight_uv, weight_type)

    #Create all of the nodes
    def makeNodes(self, tag, tag_x, tag_y, adjTag, adjCost, weight_type):
        self.nodes[tag].adjacencyList.append([self.nodes[adjTag], adjCost])
        self.nodes[adjTag].adjacencyList.append([self.nodes[tag], adjCost])
        self.nodes[tag].xcoord = tag_x
        self.nodes[tag].ycoord = tag_y
        self.nodes[tag].weight_type = weight_type

def weight(weight_type, x1, x2, y1, y2):

    if weight_type == 'EUC_2D':
        weight_xy = int(round(math.sqrt((x2 - x1)**2 + (y2 - y1)**2)))
    else:
        PI = 3.141592
        RRR = 6378.388
            # x = latitude
            # y = longitude

        new_x1 = PI * (int(x1) + 5.0 *  (x1 - int(x1)) / 3.0) / 180.0
        new_x2 = PI * (int(x2) + 5.0 *  (x2 - int(x2)) / 3.0) / 180.0
        new_y1 = PI * (int(y1) + 5.0 *  (y1 - int(y1)) / 3.0) / 180.0
        new_y2 = PI * (int(y2) + 5.0 *  (y2 - int(y2)) / 3.0) / 180.0

        q1 = math.cos( new_y1 - new_y2 ); 
        q2 = math.cos( new_x1 - new_x2 ); 
        q3 = math.cos( new_x1 + new_x2 ); 

        weight_xy = int ( RRR * math.acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0)

    return weight_xy

class ArrayHeap(list):
    def __init__(self, array):
        super().__init__(array)
        self.heapSize = len(array)

class MinHeap:
    def __init__(self, array: list, root: Node):
        self.arrayHeap = ArrayHeap(array)
        
        # Check if the root node is not the first
        # If it is not, reset the starting node by
        # removing the root node from it's original position
        # inserting it in the first position
        # and then update all indexes
        if self.arrayHeap[0] != self.arrayHeap[root.tag-1]: 
            rootNode = self.arrayHeap[root.tag-1]
            self.arrayHeap.remove(rootNode)
            self.arrayHeap.insert(0,rootNode)
            for i in range(0,self.arrayHeap.heapSize):
                self.arrayHeap[i].index = i

    def getParentIndex(self,index):
        if index%2 == 0: 
            return index//2 - 1
        else:
            return index//2

    def getLeftChildIndex(self,index):
        return 2*index+1

    def getRightChildIndex(self,index):
        return 2*index+2

    def minHeapify(self,i):
        #minHeapify is always called after we call extractMin(), in order to maintain the min heap structure
        #we always call minHeapify() FIRST with node 0, then it checks if it needs to perform a swap with the 
        #left or right child, and if so it calls minHeapify() on the newly swapped min node until nodes are arranged properly 
        #such that the root node is smaller than it's children
        #if node 0 is already the minimum node then a swap is not performed 
        l = self.getLeftChildIndex(i)
        r = self.getRightChildIndex(i)
        if l <= self.arrayHeap.heapSize-1 and self.arrayHeap[l].key < self.arrayHeap[i].key:
            min = l
        else:
            min = i
        if r <= self.arrayHeap.heapSize-1 and self.arrayHeap[r].key < self.arrayHeap[min].key:
            min = r
        if min != i:
            self.arrayHeap[i].index, self.arrayHeap[min].index = min, i # Update indexes
            self.arrayHeap[i], self.arrayHeap[min] = self.arrayHeap[min], self.arrayHeap[i]
            self.minHeapify(min)

    def shiftUp(self, index):
        parent = self.getParentIndex(index)
        current = index
        while current > 0 and self.arrayHeap[parent].key > self.arrayHeap[current].key:
            self.arrayHeap[current].index, self.arrayHeap[parent].index = parent, current # Update indexes
            self.arrayHeap[current], self.arrayHeap[parent] = self.arrayHeap[parent], self.arrayHeap[current]
            current = parent
            parent = self.getParentIndex(parent)

    def extractMin(self):
        #this function both extracts the min and pops it out of the min heap
        # in a min heap data structure, the root node at index 0 will always be the priority(minimum) 
        #so we extract the node at index zero and define it as the min to pass back
        min=self.arrayHeap[0]
        self.arrayHeap[0].isPresent = False
        #then swap the right most node and first(min) node
        self.arrayHeap[0], self.arrayHeap[self.arrayHeap.heapSize-1] = self.arrayHeap[self.arrayHeap.heapSize-1], self.arrayHeap[0]
        self.arrayHeap[0].index = 0
        #reduce the heapSize by 1 to pop out the min node
        self.arrayHeap.heapSize -=1
        #finally call minheapify() to restructure/maintain the min heap data structure
        self.minHeapify(0)
        return min
            
def Prim (G: Graph, s: Node):
    for u in G.nodes.values():
        u.key= math.inf
    s.key = 0
    MST = 0
    Q = MinHeap(list(G.nodes.values()), s)
    while Q.arrayHeap.heapSize != 0:
        u=Q.extractMin()
        for v in u.adjacencyList:
            if v[0].isPresent and v[1] < v[0].key:
                #if v[0].key is less than infinity it implies it has already been visited and added to the MST
                #we will remove v[0]key, redefine it has v[1], and replace it's value in the MST
                if v[0].key < math.inf and v[0].key != 0:
                    MST -= v[0].key
                v[0].parent = u
                v[0].key = v[1]
                Q.shiftUp(v[0].index)   
                MST += v[0].key

    return G
             
def Preorder(G: Graph, s: Node, preorder_list: list):
    #on the first call to Preorder, list preorder_list is null, and we add the starting node s
    # after the first call we all add the node next visited in the depth first search of G 
    preorder_list.append(s)

    for u in G.nodes.values():
        if u.parent == s:
            Preorder(G, u,preorder_list)

    return preorder_list

def approx_metric_tsp (G: Graph): 
    start = time.time()
    #define the starting node of cycle/root node for MST
    s = G.nodes.get(1)
    #call prims to return the MST of graph G
    MST = Prim(G, s)
    #call preorder to return a list of the vertices in G visited in a depth first search approach
    pre_order = Preorder(MST, s, [])
    tsp_cycle = pre_order
    #add the starting node to the end of the list to complete the cycle
    tsp_cycle.append(s)

    #calculate and sum the weights between each vertex in the cycle
    tsp_cost = 0
    i = 0
    j = 1
    for u in tsp_cycle:
        if tsp_cycle[i] == s and i > 0: #do nothing, at the end of the cycle
            1
        else:
            tsp_cost += weight(tsp_cycle[0].weight_type, tsp_cycle[i].xcoord, tsp_cycle[j].xcoord, tsp_cycle[i].ycoord, tsp_cycle[j].ycoord)
            i += 1
            j += 1
    end = time.time()
    time_cost = (end - start)
    return tsp_cost, time_cost

data = [
    ["tsp_dataset/burma14.tsp", 3323],
    ["tsp_dataset/ulysses16.tsp", 6859],
    ["tsp_dataset/ulysses22.tsp", 7013],
    ["tsp_dataset/eil51.tsp", 426],
    ["tsp_dataset/berlin52.tsp", 7542],
    ["tsp_dataset/kroD100.tsp", 21294],
    ["tsp_dataset/kroA100.tsp", 21282],
    ["tsp_dataset/ch150.tsp", 6528],
    ["tsp_dataset/gr202.tsp", 40160],
    ["tsp_dataset/gr229.tsp", 134602],
    ["tsp_dataset/pcb442.tsp", 50778],
    ["tsp_dataset/d493.tsp", 35002],
    ["tsp_dataset/dsj1000.tsp", 18659688]
]

results = []

for filepath, opt_solution in data:
    
    new = Graph()
    new.buildGraph(open(filepath, "r"))

    solution, time_cost = approx_metric_tsp(new) #, new.nodes.get(startingNode))
    error = ((solution-opt_solution)/opt_solution)*100

    results.append(["2-Approx",filepath, solution, time_cost,error])

    print("File name: ", filepath) 
    print("My solution: ", solution)
    print("Optimal Solution: ", opt_solution)
    print("Execution Time: ","%.5f"%(time_cost))
    print("Error: ", "%.2f"%(error), "%\n")


with open("tsp_dataset/results.txt", "w") as f:
    for line in results:
        f.write(str(line))
        f.write('\n')










