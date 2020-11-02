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
    from isoprene_rates import EXP, LOG10, TUN, ALK, NIT, ISO1, ISO2, EPO, KCO, FALL, TROE
    from read_input import background_spc
    import re
    from collections import defaultdict
    
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
    # import files
    # import numpy as np
    # from isoprene_rates import EXP, LOG10, TUN, ALK, NIT, ISO1, ISO2, EPO, KCO, FALL, TROE
    # from read_input import background_spc
    # import re
    # from collections import defaultdict
    # weight_dict = defaultdict(dict) # key: product; value: dict{reactant:weight}
    # index = 0
    # for eqn in eqns:    
    #     find_alpha_index = lambda x:re.search(r'[a-z]', x, re.I).start() # helper function
    #     initial_values_dict, TEMP = inits
    #     CFACTOR = float(initial_values_dict['CFACTOR'])
    #     reaction, k = eqn
    #     reactants = reaction.split(' = ')[0].split(' + ') # get reactants
    #     reactants = [i.strip() for i in reactants] 
    #     products = reaction.split(' = ')[1].split(' + ') # get products
    #     products = [i.strip() for i in products]

    #     # get reactants mole values
    #     reactants_mole = [float(i[:find_alpha_index(i)]) 
    #             if find_alpha_index(i)!=0 else 1 for i in reactants] # no idea of how to use this
    #     # get reactant species name 
    #     reactants_spc = [i[find_alpha_index(i):] for i in reactants]
    #     # get products mole values
    #     products_mole = [float(i[:find_alpha_index(i)]) 
    #             if find_alpha_index(i)!=0 else 1 for i in products]
    #     # get products species name 
    #     products_spc = [i[find_alpha_index(i):] for i in products]

    #     # v = 1 # assume the stoichiometric coefficient is 1 (might need to fix)
    #     SUN = 1 # random initial value for sun; need to fix !
    #     funs_temp_cf = ['ALK', 'NIT','TROE','FALL','EPO'] 
    #     funs_temp = ['TUN','ISO1','ISO2']
    #     funs_cf = ['KCO']
    #     if any([k for i in funs_temp_cf if i in k]):
    #         k = k[:-1] + ', TEMP, CFACTOR)'
    #     if any([k for i in funs_temp if i in k]):
    #         k = k[:-1] + ', TEMP)'
    #     if any([k for i in funs_cf if i in k]):
    #         k = k[:-1] + ', CFACTOR)'
    #     k_val = round(eval(k), 4)
    #     ls_concentration = []
    #     for i in reactants:
    #         if i in initial_values_dict.keys():
    #             ls_concentration.append(initial_values_dict[i])
    #         else:
    #             ls_concentration.append(initial_values_dict['ALL_SPEC'])
    #     weight_dict = defaultdict(dict) # key: product; value: dict{reactant:weight}
    #     for product, mole in zip(products_spc,products_mole):
    #         weight = mole * k_val * np.prod(ls_concentration)
    #         for reactant in reactants_spc:
    #             if reactant not in background_spc and product not in background_spc:
    #                 tup = tuple([index,reactant])
    #                 index+=1
    #                 kk = defaultdict(dict)
    #                 kk[tup][reactant] = weight
                        
    #                 if product in weight_dict.keys():
    #                     weight_dict[product].update(kk) 
    #                 else:
    #                     weight_dict[product] = kk
    #     # return(weight)
    # return weight_dict 


"""
species_a: product str
species_b: reactant str
weight_dict: a dictionary from calculate_allweight function
"""
def calculate_r(species_a, species_b, weight_dict):
    if species_a not in weight_dict.keys():
        return(species_a + ' is not in the weight dictionary.')
    b_list = [reactant for idx, reactant in weight_dict[species_a].keys()]
    if species_b not in b_list:
        return('There is no reaction to produce ' + species_a + ' from ' + species_b + '.')
    nominator = 0
    denominator = 0
    for idx, reactant in weight_dict[species_a].keys():
        denominator += abs(weight_dict[species_a][(idx, reactant)][reactant])
        if reactant == species_b:
            nominator += abs(weight_dict[species_a][(idx, reactant)][reactant])
    if denominator == 0:
        return(0)
    else:
        return(nominator/denominator)
    
def calculate_all_r(weight_dict):
    products_list = list(weight_dict.keys())
    rAB_dict = defaultdict(dict)
    for product in products_list:
        reactant_unique_list = list(set([reactant for idx, reactant in weight_dict[product].keys()]))
        for reactant in reactant_unique_list:
            rAB_dict[product][reactant] = calculate_r(product, reactant, weight_dict)
    return(rAB_dict)
