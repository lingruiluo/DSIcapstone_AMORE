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
    
    # import files
    import numpy as np
    from isoprene_rates import *
    from read_input import background_spc
    import re
    from collections import defaultdict
    
    find_alpha_index = lambda x:re.search(r'[a-z]', x, re.I).start() # helper function
    initial_values_dict, TEMP = inits
    CFACTOR = float(initial_values_dict['CFACTOR'])
    reaction, k = eqn
    reactants = reaction.split(' = ')[0].split(' + ') # get reactants 
    products = reaction.split(' = ')[1].split(' + ') # get products

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
    funs = ['TUN', 'ALK', 'NIT']
    if any([k for i in funs if i in k]):
        k = k[:-1] + ', TEMP, CFACTOR)'
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
