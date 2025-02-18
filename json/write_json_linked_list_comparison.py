import json
import os
import sys
import numpy as np

s = 0.01
L = 1

dt = 0.00005
nsave = 500
nstepT = nsave*300

data = {
    
    "name_file" : "comprarison",

    "omp": {
        "chose_nb_of_threads":False,
        "nb_of_threads":1
    },
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": nstepT,
        "dt": dt,
        "nsave": nsave,
        "kappa": 2,
        "alpha": 0.5,
        "beta": 0,
        "alpha_st": 10,
        "beta_adh": 1.2,
        "dimension": 3,
        "scheme_integration": {"Euler": True, "RK22": False},
        "comparison_algorithm": True,
    },

    "following_part": {
        "part": False,
        "min": False,
        "max": False,
        "particle": 500,
        "pressure": True,
        "rho": True,
        "position": [False, False, False],
        "velocity": [False, False, False],
    },
    
    "domain": {
        "matrix_long": [
            np.round(np.array([L, L, L]), decimals = 4).tolist(), # fluid
            

        ],
        "matrix_orig": [
            np.round(np.array([0, 0, 0]), decimals = 4).tolist(), # fluid
            

        ],
        "sphere": {
            "do": [0],
            "radius": np.round(np.array([0.3]), decimals = 4).tolist()
            },
        "vector_type": [1],
        "L_d": np.round(np.array([L, L, L]), decimals = 4).tolist(),
        "o_d": np.round(np.array([0.0, 0.0, 0.0]), decimals = 4).tolist()
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
        "u_init": np.round(np.array([0.0, 0.0, 0.0]), decimals = 4).tolist(),
        "c_0": 30,
        "gamma": 7,
        "M": 18e-3,
        "R": 8.314,
        "sigma": 52000
    },
    "forces": {
        "gravity": False,
        "surface_tension_1": False,
        "surface_tension_2": False,
        "adhesion": False
    },
    "condition": {
        "print_debug": True,
        "state_equation": {"Ideal gaz law": False, "Quasi incompresible fluid": True},
        "initial_condition": {"Hydrostatic": False, "Constant": True}
    }
}


# Do not modify what is below    
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"Other/linked_list_comparaison.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")
