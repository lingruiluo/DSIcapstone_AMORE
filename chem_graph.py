from collections import defaultdict
import operator
from graphviz import Digraph

def dfs(graph, start, visited=None):
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
    while stack:
        vertex = stack.pop()
        if visited[vertex]:
            continue
        visited[vertex] = True
        for neighbor in graph[vertex]:
            if not visited[neighbor]:
                stack.append(neighbor)
    return visited

def get_dependent_set(graph, starting_set, method='DFS'):
    '''
    Get the dependent set given the starting set
    :param graph: 
    :param starting_set: list like object containing species
    :param method: searching algorithm. Default is depth-first search
    :return: a list of species
    '''
    if method == 'DFS':
        visited = defaultdict(lambda: False)
        for specie in starting_set:
            visited = dfs(graph, specie, visited)
        return list(visited.keys())
    else:
        raise Exception('method not implemented')

def visualize_graph(graph, name='Example'):
    g = Digraph(comment=name, format='png')
    for item in graph.items():
        for rec in item[1]:
            g.edge(item[0], rec)
    return g

class ChemGraph:

    def __init__(self, dict_all_r_ab):
        self.dict = dict_all_r_ab

    def get_skeleton_graph(self, epsilon):
        graph = defaultdict(set)
        for key in self.dict:
            rec_dict = self.dict[key]
            for rec in rec_dict.keys():
                r = rec_dict[rec]
                if r >= epsilon:
                    graph[key].add(rec)
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