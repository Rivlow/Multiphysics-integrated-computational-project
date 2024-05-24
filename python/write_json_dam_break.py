import json
import os
import sys

s = 0.5
L = 1.2
dimension = 2

data = {
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": 1,
        "dt": 0.00001,
        "nsave": 100,
        "kappa": 2,
        "alpha": 0.5,
        "beta": 0,
        "alpha_st": 10,
        "beta_adh": 1.2,
        "dimension": dimension
    },
    
    "domain": {
        "matrix_long": [
            [L-2*s, s/2, L-2*s], # fluid
            [L, s/2, s/2],  
            [L, s/2, s/2],
            [s/2, s/2, L],
            [s/2, s/2, L],
            [s/2, s/2, L],
            [s/2, s/2, L-s/2],

        ],
        "matrix_orig": [
            [s, 0, s], # fluid
            [0, 0, 0], 
            [s/2, 0, s/2], 
            [-s/2, 0, s/2],
            [0, 0, s],
            [L, 0, s],
            [L+s/2, 0, 1.5*s],

        ],
        "vector_type": [1, 0, 0, 0, 0, 0, 0],
        "L_d": [2, 2, 2],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": True,
        "xyz_init": [L/2, 0, s],
        "xyz_end": [L/2, 0, L-2*s]
    },
    "thermo": {
        "rho_0": 1000,
        "rho_moving": 1000,
        "rho_fixed": 1000,
        "T": 298.15,
        "u_init": [0.0, 0.0, 0.0],
        "c_0": 30,
        "gamma": 7,
        "M": 18e-3,
        "R": 8.314,
        "sigma": 52000
    },
    "forces": {
        "gravity": False,
        "surface_tension": False,
        "adhesion": False
    },
    "condition": {
        "print_debug": False,
        "schemeIntegration": {"Euler": True, "RK22": False},
        "stateEquation": {"Ideal gaz law": False, "Quasi incompresible fluid": True},
        "initialCondition": {"Hydrostatic": False, "Constant": True}
    }
}


    
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"dam_break/{dimension}D_dam_break.json"

with open(f'{current_directory}/../tests/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")

