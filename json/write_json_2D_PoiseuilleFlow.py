import json
import os
import sys

s = 0.06
L = 1.2

nb_vtp_output = 250 # the total number of output file desired
dt = 0.0001
nsave = 100
nstepT = 200*nsave


data = {
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": nstepT,
        "dt": dt,
        "nsave": nsave,
        "kappa": 2,
        "alpha": 5,
        "beta": 0,
        "alpha_st": 10,
        "beta_adh": 1.2,
        "dimension": 2
    },
    
    "domain": {
        "matrix_long": [
            [2*L, s/2, L], # fluid
            
            [10*L, s/2, s/2], # floor 1
            [10*L, s/2, s/2], # floor 2
            [10*L, s/2, s/2], # roof wall 1
            [10*L, s/2, s/2], # roof wall 2


        ],
        "matrix_orig": [
            [s, 0, s], # fluid
            
            [-3*L, 0, 0], # floor 1
            [-3*L+s/2, 0, s/2], # floor 2
            [-3*L, 0, L + (3/2)*s], # roof 1 
            [-3*L+s/2, 0, L + 2*s], # roof 2


        ],
        "vector_type": [1, 0, 0, 0, 0],
        "L_d": [3*L, 3*L, 3*L],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": False,
        "xyz_init": [L/2, 0, s],
        "xyz_end": [L/2, 0, L-2*s]
    },
    "thermo": {
        "rho_0": 1000,
        "rho_moving": 1000,
        "rho_fixed": 1000,
        "T": 298.15,
        "u_init": [10.0, 0.0, 0.0],
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
        "print_debug": True,
        "schemeIntegration": {"Euler": True, "RK22": False},
        "stateEquation": {"Ideal gaz law": False, "Quasi incompresible fluid": True},
        "initialCondition": {"Hydrostatic": False, "Constant": True}
    }
}


# Do not modify what is below
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"PoiseuilleFlow/2D_PoiseuilleFlow.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")
