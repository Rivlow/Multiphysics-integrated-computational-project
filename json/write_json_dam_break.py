import json
import os
import sys

s = 0.01
L = 1
dimension = 2

nb_vtp_output = 250 # the total number of output file desired
dt = 0.0001
nstepT = 100000
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
            [0.4*L, s/2, 0.5*L], # fluid
            [L+s, s/2, s/2], # floor 1
            [L, s/2, s/2], # floor 2
            [s/2, s/2, L], # left wall 1
            [s/2, s/2, L], # left wall 2
            [s/2, s/2, L], # right wall 1
            [s/2, s/2, L], # right wall 2

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
        "L_d": [2*L, 2*L, 2*L],
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
        "initialCondition": {"Hydrostatic": True, "Constant": False}
    }
}


# Do not modify what is below    
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"dam_break/{dimension}D_dam_break.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")

