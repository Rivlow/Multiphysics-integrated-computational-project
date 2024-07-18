import json
import os
import sys

s = 0.02
L = 0.7
Lz = 0.62
dimension = 2

nb_vtp_output = 100 # the total number of output file desired
dt = 0.00001
nstepT = 2*25000
nsave = nb_vtp_output/(dt*nstepT)  


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
        "dimension": dimension
    },
    
    "domain": {
        "matrix_long": [
            [L-3*s, s/2, Lz-s], # fluid
            [L+s, s/2, s/2], # floor 1
            [L, s/2, s/2], # floor 2
            [s/2, s/2, Lz], # left wall 1
            [s/2, s/2, Lz], # left wall 2
            [s/2, s/2, Lz], # right wall 1
            [s/2, s/2, Lz], # right wall 2

        ],
        "matrix_orig": [
            [2*s, 0, 2*s], # fluid
            [0, 0, 0], # floor 1
            [s/2, 0, s/2], # floor 2
            [s/2, 0, 3*s/2],# left wall 1
            [0, 0, s], # left wall 2
            [L+s, 0, s], # right wall 1
            [L+s/2, 0, 3*s/2], # right wall 2

        ],
        "vector_type": [1, 0, 0, 0, 0, 0, 0],
        "L_d": [2*L, 4.1*s, 2*Lz],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": True,
        "xyz_init": [L/2, 0, 3*s],
        "xyz_end": [L/2, 0, Lz+3*s]
    },

    "following_part": {
        "part": False,
        "min": False,
        "max": True,
        "particle" : 50,
        "pressure" : 0,
        "rho" : 0,
        "position" :[0, 0, 1],
        "velocity" :[0, 0, 0]
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


# Do not modify what is below    
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"hydrostatic/{dimension}D_hydrostatic.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")

