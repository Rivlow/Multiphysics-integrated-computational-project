import json
import os
import sys

# Définir les données complètes du JSON

s = 0.05
data = {
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": 20000,
        "dt": 0.0001,
        "nsave": 100,
        "kappa": 2,
        "alpha": 0.5,
        "beta": 0,
        "alpha_st": 10,
        "beta_adh": 10,
        "dimension": 2
    },
    
    "domain": {
        "matrix_long": [
            [0.6, 0.025, 0.6],
            [1.25, 0.025, 0.02+s],
            [1.2, 0.025, 0.02+s],
            [0.05, 0.025, 1.5],
            [0.05, 0.025, 1.45-s],
            [0.05, 0.025, 1.5],
            [0.05, 0.025, 1.45-s]
        ],
        "matrix_orig": [
            [0.3, 0, 0.6],
            [0.0, 0.0, 0],
            [0.025, 0, 0.025],
            [0.0, 0.0, 0.05],
            [0.025, 0, 0.075],
            [1.25, 0, 0.05],
            [1.225, 0, 0.075]
        ],
        "vector_type": [1, 0, 0, 0, 0, 0, 0],
        "L_d": [3, 3, 3],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": True,
        "xyz_init": [0.3, 0, 0.6],
        "xyz_end": [0.9, 0, 1.2]
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
print(current_directory)

print(f"{current_directory}/..")


# Écrire les données dans un fichier JSON
json_src = "splash/2D_splash.json"
with open(f'{current_directory}/../tests/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)

print(f"Data written in '{json_src}'")

