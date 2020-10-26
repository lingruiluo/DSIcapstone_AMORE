def calculate_weight(eqn, inits):
    import numpy as np
    from isoprene_rates import EXP, LOG10, TUN, ALK, NIT
    initial_values_dict, TEMP = inits
    CFACTOR = float(initial_values_dict['CFACTOR'])
    reaction, k = eqn
    reactants = reaction.split(' = ')[0].split(' + ')
    v = 1 # assume the stoichiometric coefficient is 1 (might need to fix)
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
    weight = v * k_val * np.prod(ls_concentration)
    return(weight)
