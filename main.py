from read_input import read_eqns, read_spc, read_def, background_spc
path = 'isoprene_oxidation_model_v5_190415'
eqn_file = path+'/isoprene_full_v5.eqn'
spc_file = path+'/isoprene_full_v5.spc'
def_file = path+'/isoprene_full_v5.def'
# eqn_file = path + '/isoprene_reduced_plus_v5.eqn'
# spc_file = path + '/isoprene_reduced_plus_v5.spc'
# def_file = path + '/isoprene_reduced_plus_v5.def'
equations = read_eqns(eqn_file)
species = read_spc(spc_file)
inits = read_def(def_file)

from calculations import calculate_weight, calculate_all_weights, calculate_all_r
from direct_graph import construct_graph
# for eqn in equations:
#     weights = calculate_weight(eqn, inits) # a weight_dict for a single equation
#     print(weights)
all_weights_dict = calculate_all_weights(equations, inits)
all_r = calculate_all_r(all_weights_dict, equations)
print(all_r)
nodes, edges, graph = construct_graph(all_r, 0.001)
print("nodes: ", nodes)
print("size of nodes", len(nodes))
print("edges: ", edges)
print("size of edges", len(edges))
print("graph: ", graph)
