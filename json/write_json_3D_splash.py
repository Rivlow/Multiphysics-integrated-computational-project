import json
import os
import sys

s = 0.025
L = 1.2

nb_vtp_output = 250 # the total number of output file desired
dt = 1e-5
nstepT = 200000
nsave = 1000


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
        "dimension": 3
    },
    
    "domain": {
        "matrix_long": [
            [L/2, L/2, L/2], # fluid
            [L+3*s, L+3*s, s/2], # floor 1
            [L+2*s, L+2*s, s/2], # floor 2
            
            [s/2, L+3*s, L+2*s], # left wall 1
            [s/2, L+2*s, L+s], # left wall 2
            [s/2, L+3*s, L+2*s], # right wall 1
            [s/2, L+2*s, L+s], # right wall 2
            
            [L+s, s/2, L+2*s], # back wall 1
            [L, s/2, L+s], # back wall 2
            [L+s, s/2, L+2*s], # front wall 1
            [L, s/2, L+s], # front wall 2

        ],
        "matrix_orig": [
            [s*3/2+L/4, s*3/2+L/4, 1.5*s+L/2], # fluid
            [0, 0, 0], # floor 1
            [s/2, s/2, s/2], # floor 2
            
            [0, 0, s], # left wall 1
            [s/2, s/2, 1.5*s],# left wall 2
            [L+3*s, 0, s], # right wall 1
            [L+5*s/2, s/2, 1.5*s], # right wall 2
            
            [s, 0, s], # back wall 1
            [s*3/2, s/2, 1.5*s], # back wall 2
            [s, L+3*s, s], # front wall 1
            [3*s/2,L+5*s/2, 1.5*s], # front wall 2

        ],
        "vector_type": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        "L_d": [2*L, 2*L, 2*L],
        "o_d": [0.0, 0.0, 0.0]
    },
    "post_process": {
        "do": False,
        "xyz_init": [L/2 +3*s/2, L/2 +3*s/2, s],
        "xyz_end": [L/2 +3*s/2, L/2 +3*s/2, L]
    },

    "following_part": {
        "part": True,
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
json_src = f"splash/3D_splash.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")

