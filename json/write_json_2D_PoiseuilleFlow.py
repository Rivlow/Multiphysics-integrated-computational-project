import json
import os
import sys
import numpy as np

s = 0.02
L = 1.2

dt = 0.00001
nsave = 100
nstepT = 300*nsave


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
        "dimension": 2,
        "schemeIntegration": {"Euler": True, "RK22": False},
    },
    
    "domain": {
        "matrix_long": [

            np.round(np.array([6*L, s/2, L]),  decimals = 4).tolist(), # fluid
            np.round(np.array([10*L, s/2, s/2]),  decimals = 4).tolist(), # floor 1
            np.round(np.array([10*L, s/2, s/2]),  decimals = 4).tolist(), # floor 2
            np.round(np.array([10*L, s/2, s/2]),  decimals = 4).tolist(), # roof wall 1
            np.round(np.array([10*L, s/2, s/2]),  decimals = 4).tolist(), # roof wall 2


        ],
        "matrix_orig": [
            np.round(np.array([L, 0, s]),  decimals = 4).tolist(), # fluid
            np.round(np.array([0, 0, 0]),  decimals = 4).tolist(), # floor 1
            np.round(np.array([s/2, 0, s/2]),  decimals = 4).tolist(), # floor 2
            np.round(np.array([s/2, 0, L + (3/2)*s]),  decimals = 4).tolist(), # roof 1 
            np.round(np.array([0, 0, L + 2*s]),  decimals= 4 ).tolist(), # roof 2


        ],
        "vector_type": [1, 0, 0, 0, 0],
        "L_d": np.round(np.array([9*L, 3*L, 3*L]),  decimals = 4).tolist(),
        "o_d": np.round(np.array([0.0, 0.0, 0.0]),  decimals = 4).tolist()
    },
    "post_process": {
        "do": True,
        "xyz_init": np.round(np.array([4*L, 0, s]),  decimals = 4).tolist(),
        "xyz_end": np.round(np.array([4*L, 0, L + 2*s]),  decimals = 4).tolist()
    },
    "thermo": {
        "rho_0": 1000,
        "rho_moving": 1000,
        "rho_fixed": 1000,
        "T": 298.15,
        "u_init": [-0.5, 0.0, 0.0],
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
