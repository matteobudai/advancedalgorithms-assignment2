import math
import copy
from tracemalloc import start
import time

#Transform coordinate data into connected graph (list in format [u, v, weight])
class Graph:

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
        connected_graph = []
        vertex= []
        
        
        #create a connected graph by looping through each node and connecting it to all other nodes
        #calculate the weight between each node depending on the weight_type defined in the graph
        for iterate in range(node_data_start+1, node_data_end):
            
            node_u = int(lines[iterate].split()[0])
            u_x1_coord = float(lines[iterate].split()[1])
            u_y1_coord = float(lines[iterate].split()[2])
            vertex.append(node_u)

            for interate_v in range(node_data_start+1, node_data_end):
                node_v = int(lines[interate_v].split()[0])
                if node_u == node_v: #do nothing
                    1
                elif node_u != node_v: 
                    #node_v = int(lines[interate_v].split()[0])
                    v_x2_coord = float(lines[interate_v].split()[1])
                    v_y2_coord = float(lines[interate_v].split()[2])
                    weight_uv = self.weight(weight_type, u_x1_coord, v_x2_coord, u_y1_coord, v_y2_coord)
                    connected_graph.append([node_u,node_v,weight_uv])
            
        return(connected_graph, vertex)


    def weight(self, weight_type, x1, x2, y1, y2):

        if weight_type == 'EUC_2D':
           weight = int(round(math.sqrt((x2 - x1) **2 + (y2 - y1) **2))) 
        elif weight_type == 'GEO':
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

            weight = int ( RRR * math.acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0)

        return weight

# Function to get the nearest neighbor
# Input: path
# Ouput: nearest neighbor node, weight of edge 
def minDist(p):
    # Initially setup of the weight 
    dist= math.inf
    # Take the last node of the path
    node= p[-1]
    # Take the length of the vertexes
    l=len(vertex)
    # set the node to None
    vertexF=None
    #set the index where start to find the nearest neighbor
    i=(l-1)*(node-1)
    # While we have visited all neighbors of node
    while(i<(l-1)*node):
        # if the weigh is lower and the node is not in the path
        if graph[i][2]<dist and graph[i][1] not in p:
            # Save lowest weight
            dist=graph[i][2]
            # Save node connected by lowest cost edge
            vertexF=graph[i][1]
        i=i+1
    return dist,vertexF


# Function for Nearest Neighbor algorithm 
# Output: Weight of all path and time of algorithm ran
def NearestNeighbor():
     # Start time point
    start= time.time()
    # Create a list of all nodes, which will be used as visited list
    unvisited = copy.deepcopy(vertex)
    # Add first node to the final path and remove it from the unvisited list
    startingNode=unvisited[0]
    # Remove the starting Node from the unvisited
    unvisited.remove(startingNode)
    # Create the initially path
    path=[startingNode]
    # counter of the weight of the path
    dist=0
    # While we have unvisited nodes do
    while len(unvisited) >0:
        #call the minDist function to get the neirest neighbor node and its weight
        d, vertexF= minDist(path)
        # remove from the unvisited the node found
        unvisited.remove(vertexF)
        # add to the path the new node
        path.append(vertexF)
        # Updatge the total weight
        dist=dist+d
    # update the weight with the distance between the starting node and the last of the path
    dist=dist+graph[path[-1]-2][2]
    # # Add starting node to the end because it should be a hamiltonian circuit
    path.append(startingNode)
    # End time point
    end = time.time()
    # Count time of running 
    time_cost =  end - start
    # Return distance and time required 
    return dist, time_cost

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
    graph, vertex = Graph().buildGraph(open(filepath, "r"))
    solution, time_cost= NearestNeighbor()
    error = ((solution-opt_solution)/opt_solution)*100

    results.append(["Nearest Neighbor",filepath, solution, time_cost,error])

    print("File name: ", filepath) 
    print("My solution: ", solution)
    print("Optimal Solution: ", opt_solution)
    print("Execution Time: ","%.5f"%(time_cost))
    print("Error: ", "%.2f"%(error), "%\n")

with open("tsp_dataset/results.txt", "a") as f:
    for line in results:
        f.write(str(line))
        f.write('\n')
