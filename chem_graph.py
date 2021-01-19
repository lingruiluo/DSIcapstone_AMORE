from collections import defaultdict
import operator
from graphviz import Digraph

def dfs(graph, start, visited=None, reduced_graph=None):
    '''
    Depth-First Search Non-Recursive Function
    :param graph: dictionary 
    :param start: starting node for the search
    :param visited: dictionary for already visited node
    :return: updated dictionary
    '''
    if visited == None:
        visited = defaultdict(lambda: False)
    if visited[start]:
        return visited
    stack = [start]

    if reduced_graph == None:
        reduced_graph = defaultdict(lambda : set())
    while stack:
        vertex = stack.pop()
        if visited[vertex]:
            continue
        visited[vertex] = True
        for s in graph[vertex]:
            neighbor = s[0]
            r = s[1]
            reduced_graph[vertex].add([neighbor, r])
            if not visited[neighbor]:
                stack.append(neighbor)
    return visited, reduced_graph


def visualize_graph(graph, node_to_mark=None, name='Example'):
    g = Digraph(comment=name, format='png')
    if node_to_mark:
        g.attr('node', shape='doublecircle', style='filled', fillcolor='red')
        for n in node_to_mark:
            g.node(n)
    
    g.attr('node', shape='ellipse', fillcolor='white')
    for item in graph.items():
        g.node(item[0])
        for s in item[1]:
            g.node(s[0])
    
    for item in graph.items():
        g.node(item[0])
        for s in item[1]:
            g.edge(item[0], s[0])

    return g

class ChemGraph:

    def __init__(self, dict_all_r_ab, starting_set, must_contain=None):
        self.dict = dict_all_r_ab
        self.starting_set = starting_set
        self.must_contain = must_contain

    def get_skeleton_graph(self, epsilon):
        graph = defaultdict(set)
        for key in self.dict:
            rec_dict = self.dict[key]
            for rec in rec_dict.keys():
                r = rec_dict[rec]
                if rec in self.must_contain:
                    graph[key].add((rec, r))
                    continue
                if r >= epsilon:
                    graph[key].add((rec, r))
        return graph

    def get_all_skeleton_graph(self, epsilons, method='naive'):
        ret_dict = {}
        if method == 'naive':
            for epsilon in epsilons:
                ret_dict[epsilon] = self.get_skeleton_graph(epsilon)
            return ret_dict
        elif method == 'incrementive':
            epsilons = sorted(epsilons)
            tups = [(d[0], rec[0], rec[1]) for d in self.dict.items() for rec in d[1].items()]
            sorted_tups = sorted(tups, key=operator.itemgetter(2))
            start_index = 0
            n = len(sorted_tups)
            for e in epsilons:
                while start_index < n:
                    if sorted_tups[start_index][2] < e:
                        break
                    start_index += 1
                included_tup = sorted_tups[0:start_index]
                temp_graph = defaultdict(set)
                for t in included_tup:
                    temp_graph[t[0]].add(t[1])
                ret_dict[e] = temp_graph
            return ret_dict
        else:
            raise Exception('method not implemented')

    def get_dependent_set(self, epsilon, skeleton_graph=None, update=True):
        if skeleton_graph == None:
            skeleton_graph = self.get_skeleton_graph(epsilon)
        graph = skeleton_graph
        visited = defaultdict(lambda: False)
        for specie in self.starting_set:
            visited, reduced_graph = dfs(graph, specie, visited)
        if update:
            self.skeleton_graph = skeleton_graph
            self.reduced_graph = reduced_graph
            self.dependent_set = visited
        return reduced_graph
    
    def get_all_reduced_graph(self, epsilons):
        skeleton_graph_dict = self.get_all_skeleton_graph(epsilons)
        ret_dict = {}
        for e in epsilons:
            sk_graph = skeleton_graph_dict[e]
            rd_graph = self.get_dependent_set(e, sk_graph, False)
            ret_dict[e] = rd_graph
        return ret_dict
