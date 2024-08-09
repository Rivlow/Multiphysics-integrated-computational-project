import json
import os
import sys
import numpy as np

s = 0.01
L = 1.2

nb_vtp_output = 250 # the total number of output file desired
dt = 0.0001
nstepT = 100000
nsave = 500 

data = {
    
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
        "dimension": 2,
        "schemeIntegration": {"Euler": True, "RK22": False},
        "comparison_algorithm": False,
    },

    "following_part": {
        "part": True,
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
            np.round(np.array([L/2, s/2, L/2]), decimals = 4).tolist(), # fluid
            np.round(np.array([L+3*s, s/2, s/2]), decimals = 4).tolist(), # floor 1
            np.round(np.array([L+2*s, s/2, s/2]), decimals = 4).tolist(), # floor 2
            np.round(np.array([s/2, s/2, 2*L]), decimals = 4).tolist(), # left wall 1
            np.round(np.array([s/2, s/2, 2*L+s]), decimals = 4).tolist(), # left wall 2
            np.round(np.array([s/2, s/2, 2*L+s]), decimals = 4).tolist(), # right wall 1
            np.round(np.array([s/2, s/2, 2*L]), decimals = 4).tolist(), # right wall 2
            np.round(np.array([L+s,s/2,s/2]), decimals = 4).tolist(), # upper roof
            np.round(np.array([L,s/2,s/2]), decimals = 4).tolist(), # lower roof


        ],
        "matrix_orig": [
            np.round(np.array([s*3/2+L/4, 0, 1.5*s+L/2]), decimals = 4).tolist(), # fluid
            np.round(np.array([0, 0, 0]), decimals = 4).tolist(), # floor 1
            np.round(np.array([s/2, 0, s/2]), decimals = 4).tolist(), # floor 2
            np.round(np.array([s/2, 0, 3*s/2]), decimals = 4).tolist(),# left wall 1
            np.round(np.array([0, 0, s]), decimals = 4).tolist(), # left wall 2
            np.round(np.array([L+3*s, 0, s]), decimals = 4).tolist(), # right wall 1
            np.round(np.array( [L+5*s/2, 0, 3*s/2]), decimals = 4).tolist(), # right wall 2
            np.round(np.array( [s,0,2*L+2*s]), decimals = 4).tolist(), # right wall 2
            np.round(np.array( [3*s/2,0,2*L+3*s/2]), decimals = 4).tolist(), # right wall 2

        ],
        "sphere": {
                "do": [0, 0, 0, 0, 0, 0, 0, 0, 0],
                "radius": [0.3]
            },
        "vector_type": [1, 0, 0, 0, 0, 0, 0, 0, 0],
        "L_d": np.round(np.array([3*L, 4*s, 2*L]), decimals = 4).tolist(),
        "o_d": np.round(np.array([0.0, 0.0, 0.0]), decimals = 4).tolist()
    },
    "post_process": {
        "do": False,
        "xyz_init": np.round(np.array([L/2, 0, s]), decimals = 4).tolist(),
        "xyz_end": np.round(np.array([L/2, 0, L-2*s]), decimals = 4).tolist()
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
        "gravity": True,
        "surface_tension_1": False,
        "surface_tension_2": False,
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
json_src = f"splash/2D_splash.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")
