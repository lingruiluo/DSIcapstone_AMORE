import numpy as np
import re
from collections import defaultdict
from isoprene_rates import EXP, LOG10, TUN, ALK, NIT, ISO1, ISO2, EPO, KCO, FALL, TROE
from read_input import background_spc

def get_reactants(eqn):
	'''Get reactants from an equation (ignore background species)
    Parameters
    ----------
    eqn: tuple
        An equation.
    
    Returns
    ----------
    reactants_spc: tuple
        A tuple of reactants in the equation
    '''
    import re
    from read_input import background_spc
    find_alpha_index = lambda x:re.search(r'[a-z]', x, re.I).start() # helper function
    reaction, k = eqn
    reactants = reaction.split(' = ')[0].split(' + ') # get reactants
    reactants = [i.strip() for i in reactants]
    # get reactant species name
    reactants_spc = [i[find_alpha_index(i):] for i in reactants]
    reactants_spc = [r for r in reactants_spc if r not in background_spc]
    return(tuple(reactants_spc))

def get_products(eqn):
	'''Get products from an equation (ignore background species)
    Parameters
    ----------
    eqn: tuple
        An equation.
    
    Returns
    ----------
    products_spc: tuple
        A tuple of products in the equation
    '''
    import re
    from read_input import background_spc
    find_alpha_index = lambda x:re.search(r'[a-z]', x, re.I).start() # helper function
    reaction, k = eqn
    products = reaction.split(' = ')[1].split(' + ') # get products
    products = [i.strip() for i in products]
    # get products species name
    products_spc = [i[find_alpha_index(i):] for i in products]
    products_spc = [r for r in products_spc if r not in background_spc]
    return(tuple(products_spc))

def get_properties(eqns):
	'''Get reactants and products from equations
    Parameters
    ----------
    eqns: list
        A List of tuples storing equations
    
    Returns
    ----------
    ret: dict
        A dict of dicts storing reactants and products for each equation
    '''
    from collections import defaultdict
    ret = defaultdict(dict)
    for i in range(0, len(eqns)):
        eqn = eqns[i]
        reactants = get_reactants(eqn)
        products = get_products(eqn)
        ret[i]['reactants'] = reactants
        ret[i]['products'] = products
    return(ret)

'''
if a species is a reactant for the equation, 'r' is marked
if a species is a product for the equation, 'p' is marked
'''
def get_eqns_involve_species(species, eqns):
	'''Get equations involving specified sepcies. 
	If a species is a reactant for the equation, 'r' is marked;
	if a species is a product for the equation, 'p' is marked.

    Parameters
    ----------
    species: str
        Species
    eqns: list
        A List of tuples storing equations
    
    Returns
    ----------
    eqns_idx: list
        A list of tuples storing indicators of reactants or products and index of equations.
    '''
    eqns_idx = []
    properties_dict = get_properties(eqns)
    for i in range(0, len(eqns)):
        
        if species in properties_dict[i]['products']:
            eqns_idx.append(('p', i))
        if species in properties_dict[i]['reactants']:
            eqns_idx.append(('r', i))
    return(eqns_idx)

def calculate_weight(eqn, inits):
    '''Calculate weights for a given equation
    Parameters
    ----------
    eql: tuple
        The first element of the tuple is an equation. The second element is reaction rate.
    inits: dict
        A dictionary storing some initial values
    
    Returns
    ----------
    weight_dict: dict
        A dictionary storing the weights
    '''

    find_alpha_index = lambda x:re.search(r'[a-z]', x, re.I).start() # helper function
    initial_values_dict, TEMP = inits
    CFACTOR = float(initial_values_dict['CFACTOR'])
    reaction, k = eqn
    reactants = reaction.split(' = ')[0].split(' + ') # get reactants
    reactants = [i.strip() for i in reactants] 
    products = reaction.split(' = ')[1].split(' + ') # get products
    products = [i.strip() for i in products]

    # get reactants mole values
    reactants_mole = [float(i[:find_alpha_index(i)]) 
            if find_alpha_index(i)!=0 else 1 for i in reactants] # no idea of how to use this
    # get reactant species name 
    reactants_spc = [i[find_alpha_index(i):] for i in reactants]
    # get products mole values
    products_mole = [float(i[:find_alpha_index(i)]) 
            if find_alpha_index(i)!=0 else 1 for i in products]
    # get products species name 
    products_spc = [i[find_alpha_index(i):] for i in products]

    # v = 1 # assume the stoichiometric coefficient is 1 (might need to fix)
    SUN = 1 # random initial value for sun; need to fix !
    funs_temp_cf = ['ALK', 'NIT','TROE','FALL','EPO'] 
    funs_temp = ['TUN','ISO1','ISO2']
    funs_cf = ['KCO']
    if any([k for i in funs_temp_cf if i in k]):
        k = k[:-1] + ', TEMP, CFACTOR)'
    if any([k for i in funs_temp if i in k]):
        k = k[:-1] + ', TEMP)'
    if any([k for i in funs_cf if i in k]):
        k = k[:-1] + ', CFACTOR)'
    k_val = round(eval(k), 4)
    ls_concentration = []
    for i in reactants:
        if i in initial_values_dict.keys():
            ls_concentration.append(initial_values_dict[i])
        else:
            ls_concentration.append(initial_values_dict['ALL_SPEC'])
    # weight = v * k_val * np.prod(ls_concentration)
    weight_dict = defaultdict(dict) # key: product; value: dict{reactant:weight}
    for product, mole in zip(products_spc,products_mole):
        weight = mole * k_val * np.prod(ls_concentration)
        for reactant in reactants_spc:
            if reactant not in background_spc and product not in background_spc:
                weight_dict[product][reactant] = weight
    # return(weight)
    return weight_dict  


def calculate_all_weights(eqns, inits):
    '''Calculate all weights for the input
    Parameters
    ----------
    eql: tuple
        The first element of the tuple is an equation. The second element is reaction rate.
    inits: dict
        A dictionary storing some initial values
    
    Returns
    ----------
    all_weight_dict: dict
        A dictionary storing the weights
    '''
    from collections import defaultdict
    all_weight_dict = defaultdict(dict)
    for i in range(len(eqns)):
        weight_dict = calculate_weight(eqns[i],inits)
        products = [i for i in weight_dict.keys()]
        inner_dict = defaultdict(dict)
        for product in products:
            reactants = [i for i in weight_dict[product]]
            inner_dict_key = (i, tuple(reactants))
            for reactant in reactants:
                inner_dict[reactant] = weight_dict[product][reactant]
            all_weight_dict[product][inner_dict_key] = inner_dict
    return all_weight_dict

    
def get_weight(eqns, idx, all_weights_dict, reactant_spc):
    eqn = eqns[idx]
    products = get_products(eqn)
    all_weight = 0.0
    for i in products:   ## maybe contains mutiple products, all weights from all products
        for key in list(all_weights_dict[i].keys()):
            if idx in key: 
                all_weight+=all_weights_dict[i][key][reactant_spc]
    return(all_weight)

"""
species_a: product str (main species)
species_b: reactant str
weight_dict: a dictionary from calculate_allweight function
"""
def calculate_r(species_a, species_b, all_weights_dict, eqns):
    a_eqns = get_eqns_involve_species(species_a, eqns)
    b_eqns = get_eqns_involve_species(species_b, eqns)
    a_idx = [eqn_idx for (species_type, eqn_idx) in a_eqns]
    b_idx = [eqn_idx for (species_type, eqn_idx) in b_eqns]
    numerator_eqn_list = a_idx and b_idx # eqns that species a and species b both involve in
    if len(numerator_eqn_list) == 0:
        return('There is no reaction to produce ' + species_a + ' from ' + species_b + '.')
    numerator = 0.0
    for i in numerator_eqn_list:
        for (species_type, eqn_idx) in a_eqns:
            if eqn_idx == i and species_type == 'p':
                for t in all_weights_dict[species_a].keys():
                    idx = t[0] # first position is the eqn_idx
                    if idx == eqn_idx:
                        reactants = list(t[1])
                        value = list(all_weights_dict[species_a][t].values())
                        value = [abs(v) for v in value]
                        numerator += abs(sum(value))
    
    denominator = 0.0
    # species_a as products
    for t in all_weights_dict[species_a].keys():
        reactants = t[1] # second position in the tuple
        reactants = list(reactants)
        for reactant in reactants:
            denominator += abs(all_weights_dict[species_a][t][reactant])
    # species_a as reactants
    for r in a_eqns:
        if species_type == 'r':
            denominator += get_weight(eqns, eqn_idx, all_weights_dict, species_a)
    
    if denominator == 0:
        return(0)
    else:
        return(numerator/denominator)
 

def calculate_all_r(all_weights_dict, eqns):
    from collections import defaultdict
    products_list = list(all_weights_dict.keys())
    rAB_dict = defaultdict(dict)
    for product in products_list:
        reactant_unique_list = list(set([i for key in all_weights_dict[product].keys() for i in list(key[1])]))
        for reactant in reactant_unique_list:
            rAB_dict[product][reactant] = calculate_r(product, reactant, all_weights_dict, eqns)
    return(rAB_dict)
