import json
import os
import sys

# Définir les données complètes du JSON

s = 0.001
L = 0.12

data = {
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": 100000,
        "dt": 0.00001,
        "nsave": 100,
        "kappa": 2,
        "alpha": 0.5,
        "beta": 0,
        "alpha_st": 10,
        "beta_adh": 1.2,
        "dimension": 2
    },
    
    "domain": {
        "matrix_long": [
            [L/2, s/2, L/2], # fluid
            [L+s, s/2, s/2],  
            [L, s/2, s/2],
        ],
        "matrix_orig": [
            [L/4+s/2, 0, L/4], # fluid
            [0.0, 0, 0], 
            [s/2, 0, s/2],  
        ],
        "vector_type": [1, 0, 0],
        "L_d": [(L+s)*1.1, s, (L+s)*1.1],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": False,
        "xyz_init": [s, 0, L/5],
        "xyz_end": [L+s, 0, L/5]
    },
    "thermo": {
        "rho_0": 1000,
        "rho_moving": 1000,
        "rho_fixed": 1000,
        "T": 298.15,
        "u_init": [0.0, 0.0, 0.0],
        "c_0": 340,
        "gamma": 7,
        "M": 18e-3,
        "R": 8.314,
        "sigma": 52000
    },
    "forces": {
        "gravity": True,
        "surface_tension": True,
        "adhesion": True
    },
    "condition": {
        "print_debug": False,
        "schemeIntegration": {"Euler": True, "RK22": False},
        "stateEquation": {"Ideal gaz law": False, "Quasi incompresible fluid": True},
        "initialCondition": {"Hydrostatic": False, "Constant": True}
    }
}

current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
print(current_directory)

print(f"{current_directory}/..")


# Écrire les données dans un fichier JSON
json_src = "surface_tension/2D_surface_tension.json"
with open(f'{current_directory}/../tests/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)

print(f"Data written in '{json_src}'")

