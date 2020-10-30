# Clean up the data based on the formats of equations above

background_spc = ['CO2','CS2','CO','H2CO3','COS','M','O2','OH','HO2',
                  'NO','NO2','NO3','MO2','H2O','H2O2','H2Od','O3',
                  'SO2','H2S','HF','H2S','N2H4','HN3','HI','HBr',
                  'HCl','HCN','H2Se','H2Te','NH2OH','HBrO','HClO','H3PO2','HPO3',
                  'H2O3','OF2','O2F2','NOHSO4','N2O5','N2F4','N2O4','N2O3','HNO3',
                  'HNO2','N2O','NF5','NI3','H2S','H2SO4','H2SO3','SO2Cl2','S4N4','H2SO5',
                  'H2S2O7','S2F10','H3NO3S','Br2S','SF6','SF4'
    ]

def read_eqns(eqn_file):
    equations = None
    with open(eqn_file,'r') as f:
        lines = f.readlines()
        equantions = lines[:]
    equations = [i.strip() for i in equantions[1:]]
    equations = [tuple(i.split(':')) for i in equations if len(i)>0]
    equations = [(i[0].strip(),i[1].strip().strip(';')) for i in equations]
    return(equations)

# Process the raw input into species
def read_spc(spc_file):
    species = None
    with open(spc_file,'r') as f:
        lines = f.readlines()
        species = lines[:]
    species = [s.split('=')[0].strip() for s in species]
    species = [s for s in species if s and s[0]!='#' and s not in background_spc]
    return(species)

# Get initial values for species and temperature
def read_def(def_file):
    import re
    init_values = None
    TEMP = None
    with open(def_file, 'r') as f:
        lines = f.readlines()
        TEMP = [line.strip() for line in lines if line.strip().startswith('TEMP')]
        idx = [i for i in range(len(lines)) if (lines[i].startswith('#INITVALUES')) or (lines[i].startswith('#INLINE F90'))]
        init_values = lines[idx[0]+2:idx[1]-1]

    init_values = [i.strip().split('=') for i in init_values]
    init_values_dict = {}
    for c, val in init_values:
        init_values_dict[c.strip()] = float(re.sub(';', '', val).strip())
    TEMP = float(TEMP[0].split(' = ')[1])
    return(init_values_dict, TEMP)
