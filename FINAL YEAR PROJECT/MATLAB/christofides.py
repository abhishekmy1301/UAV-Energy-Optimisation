import networkx as nx
import random
import math


n = 20                     # number of cities
G = nx.complete_graph(n)  # graph with a vertex for each city
# for convenience, pick the city (x,y)-coordinates at random
import random

def getCost(mypos,tour):

    cost=0

    for i in range(len(tour)-1):
        print(i)
        cost=cost+eucl_dist(mypos[i][0],mypos[i][1],mypos[i+1][0],mypos[i+1][1])
        print(cost)
    return cost




my_pos = { i : ( random.random()*1.0, random.random()*1.0 ) for i in G.nodes } # pos[i] = (x_i, y_i)
nx.draw(G, pos=my_pos)

# for convenience, suppose that distances are Euclidean
import math
def eucl_dist(x1,y1,x2,y2):
    return math.sqrt( (x1-x2)**2 + (y1-y2)**2 )

for i,j in G.edges:
    (x1,y1) = my_pos[i]
    (x2,y2) = my_pos[j]
    G.edges[i,j]['length'] = eucl_dist(x1,y1,x2,y2)
# find minimum spanning tree
T = nx.minimum_spanning_tree(G,weight='length')
nx.draw(T, pos=my_pos)

# identify the odd-degree nodes
odd_degree_nodes = [ i for i in T.nodes if T.degree(i) % 2 ]
node_colors = [ T.degree(i) % 2 for i in T.nodes ]
nx.draw(T, pos=my_pos, node_color=node_colors)

# find a minimum-cost perfect matching over the odd-degree nodes
for i,j in G.edges:
    G.edges[i,j]['neg_length'] = - G.edges[i,j]['length']
    
matching = nx.max_weight_matching( G.subgraph(odd_degree_nodes), maxcardinality=True, weight='neg_length')
print(matching)


# draw the matching
nx.draw(G.edge_subgraph(matching),pos=my_pos)



# create a multigraph with edge_set = (spanning tree edges) + (matching)
M = nx.MultiGraph()

M.add_nodes_from(range(n))

M.add_edges_from(T.edges())
M.add_edges_from(matching)

print(M.edges())
print("M has this many edges =",M.number_of_edges())


# find an Eulerian cycle of the multigraph
initial_tour = list ( nx.eulerian_circuit(M,source=0) )
print(initial_tour)


# take shortcuts (avoid repeated nodes)
tour = [ 0 ]
for (i,j) in initial_tour:
    if j not in tour:
        tour.append(j)
print(tour)

# draw the tour
tour_edges = [ (tour[i-1],tour[i]) for i in range(n) ]
nx.draw(G.edge_subgraph(tour_edges), pos=my_pos)

print(getCost(my_pos,tour))