from collections import defaultdict

class ChemGraph:

    def __init__(self, dict, all_r_ab):
        self.dict = dict
        self.graph = defaultdict(set)
        self.all_r_ab = all_r_ab

    def construct_graph(self, epsilon):
        for key in self.dict:
            rec_dict = self.all_r_ab[key]
            for rec in rec_dict:
                r = rec_dict[rec]
                if r >= epsilon:
                    self.graph[key].add(rec)

    def get_dependent_set(self, starting_set, method='DFS'):
        '''
        Get the dependent set given the starting set
        :param starting_set: list like object containing species
        :param method: searching algorithm. Default is depth-first search
        :return: a list of species
        '''
        if method == 'DFS':
            visited = defaultdict(lambda: False)
            for specie in starting_set:
                visited = self.dfs(specie, visited)
            return list(visited.keys())
        else:
            raise Exception('method not implemented')

    def dfs(self, start, visited=None):
        '''
        Depth-First Search Non-Recursive Function
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
            for neighbor in self.graph[vertex]:
                stack.append(neighbor)
        return visited
