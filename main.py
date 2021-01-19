from read_input import read_eqns, read_spc, read_def, background_spc
from calculations import calculate_weight, calculate_all_weights, calculate_all_r
from chem_graph import ChemGraph
from direct_graph import construct_graph

path = 'isoprene_oxidation_model_v5_190415'
eqn_file = path+'/isoprene_full_v5.eqn'
spc_file = path+'/isoprene_full_v5.spc'
def_file = path+'/isoprene_full_v5.def'
equations = read_eqns(eqn_file)
species = read_spc(spc_file)
inits = read_def(def_file)

# for eqn in equations:
#     weights = calculate_weight(eqn, inits) # a weight_dict for a single equation
#     print(weights)
all_weights_dict = calculate_all_weights(equations, inits, SUN=0.5)
all_r = calculate_all_r(all_weights_dict, equations)
# print(all_r)
nodes, edges, graph = construct_graph(all_r, 0.0)
# print("nodes: ", nodes)
print("size of nodes", len(nodes)) # 382
# print("edges: ", edges)
print("size of edges", len(edges)) # 864
# print("graph: ", graph)


starting_set = []
must_contain = None
chem_graph = ChemGraph(all_r, starting_set, must_contain)


## one epsilon
epsilon = 0.5
reduced_graph = chem_graph.get_dependent_set(epsilon)
## skeleton graph, reduced_graph, dependent set will also be stored within the object
skeleton_graph = chem_graph.skeleton_graph
reduced_graph = chem_graph.reduced_graph
dependent_set = chem_graph.dependent_set


## multiple epsilons
## get a dictionary {epsilon:reduced graph}
epsilons = [0,0.01,0.05,0.1,0.2,0.5] 
skeleton_graphs = chem_graph.get_all_skeleton_graph(epsilons)
reduced_graphs = chem_graph.get_all_reduced_graph(epsilon)






def get_edges_nodes(graph):
    edges, nodes = set(), set()
    for key in graph:
        values = graph[key]
        edge = set((key,val) for val in values)
        edges.update(edge)
        nodes.add(key)
        nodes.update(values)
    return edges, nodes
graph_size = []
for eps in epsilons:
    graph = graphs[eps]
    edges, nodes = get_edges_nodes(graph)
    graph_size.append((len(edges),len(nodes)))
print(graph_size) 
# edges_size = []
# nodes_size = []
# for i in graph_size:
#     edges_size.append(i[0])
#     nodes_size.append(i[1])
# import matplotlib.pyplot as plt
# plt.plot(epsilons,nodes_size,'-b.',color='r',label='nodes_size')


