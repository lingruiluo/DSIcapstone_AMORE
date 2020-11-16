from collections import defaultdict

def construct_graph(all_r_ab, epsilon):
    '''construct graph
    Parameters
    ----------
    all_r_ab: dict
        key: product
        value: dict
            key: reactant
            value: r_ab
    epsilon: float
        threshold value for graph weight
    
    Returns
    ----------
    nodes: list
        nodes of graph
    edges: list
        edges of graph; (product, reactant)
    graph: dict
        key: product
        value: set of reactant
    
    Rules:
    (1) Each node in DRG is uniquely mapped to a species in the detailed mechanism.
    (2) There exists a directed edge from A to B if and only if r_AB is larger than or equal to epsilon.
    '''
    nodes = set()
    edges = set()
    graph = defaultdict(set)
    for key in all_r_ab:
        add_key = False
        rec_dict = all_r_ab[key]
        for rec in rec_dict:
            r = rec_dict[rec]
            if r >= epsilon:
                edges.add((key, rec))
                nodes.add(rec)
                graph[key].add(rec)
                add_key = True
        if add_key:
            nodes.add(key)
    
    return list(nodes), list(edges), graph