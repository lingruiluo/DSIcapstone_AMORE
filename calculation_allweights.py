#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 13:49:47 2020

@author: zijingwang
"""

def calculate_allweight(eqns, inits):
    import numpy as np
    from isoprene_rates import EXP, LOG10, TUN, ALK, NIT
    from read_input import background_spc
    import re
    from collections import defaultdict
    weight_dict = defaultdict(dict) # key: product; value: dict{reactant:weight}
    index = 0
    for eqn in eqns:
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
        
        for product, mole in zip(products_spc,products_mole):
            weight = mole * k_val * np.prod(ls_concentration)
            for reactant in reactants_spc:
                if reactant not in background_spc and product not in background_spc:
#                         if product not in index_dict.keys():
#                             index = 1
#                             index_dict[product] = index
#                         else:
#                             index = index_dict[product]+1
#                             index_dict[product] = index
#                             print(product,index)
                    tup = tuple([index,reactant])
                    index+=1
                    kk = defaultdict(dict)
                    kk[tup][reactant] = weight
                        
                    if product in weight_dict.keys():
                        weight_dict[product].update(kk) 
                    else:
                        weight_dict[product] = kk
    return weight_dict 