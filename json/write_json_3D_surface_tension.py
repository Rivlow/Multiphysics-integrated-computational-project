import json
import os
import sys
import numpy as np

s = 0.015
L = 0.1
dt = 1e-5
nsave = 1000 
nstepT = nsave*200


data = {
    
    "name_file" : "3D_surface_tension",
    
    "simulation": {
        "theta": 0.5,
        "s": s,
        "nstepT": nstepT,
        "dt": dt,
        "nsave": nsave,
        "kappa": 2,
        "alpha": 0.5,
        "beta": 0,
        "alpha_st": 1,
        "beta_adh": 1.2,
        "dimension": 3,
        "schemeIntegration": {"Euler": True, "RK22": False},
        "comparison_algorithm": False,
    },

    "following_part": {
        "part": False,
        "min": True,
        "max": True,
        "particle": 500,
        "pressure": False,
        "rho": False,
        "position": [True, True, True],
        "velocity": [False, False, False],
    },
    

  "domain":{
    "matrix_long" : np.round(np.array([[L, L, L]]), decimals = 7).tolist(),
    "matrix_orig" : np.round(np.array([[L, L, L]]), decimals = 7).tolist(),
    "sphere": {
                "do": [0],
                "radius": [L/2]
              },
    "vector_type" : [1],
    "L_d": np.round(np.array([3*L, 3*L, 3*L]), decimals = 7).tolist(),
    "o_d": np.round(np.array([0.0, 0.0, 0.0]), decimals = 7).tolist()
  },

  "post_process":{
    "do": True,
    "xyz_init": np.round(np.array([0, 1.5*L, 1.5*L]), decimals = 7).tolist(),
    "xyz_end": np.round(np.array([3*L, 1.5*L, 1.5*L]), decimals = 7).tolist()
  },


  "thermo":{
    "rho_0": 1000,
    "rho_moving": 1000,
    "rho_fixed" : 1000,
    "T": 298.15,
    "u_init": np.round(np.array([0.0, 0.0, 0.0]), decimals = 7).tolist(),
    "c_0": 10, 
    "gamma": 7, 
    "M": 18e-3, 
    "R":8.314,
    "sigma":0.1
  },


  "forces":{
    "gravity":False,
    "surface_tension_1": True,
    "surface_tension_2": False,
    "adhesion":False
  },

  "condition":{
    "print_debug":False,
    "stateEquation" : {"Ideal gaz law":False, "Quasi incompresible fluid":True},
    "initialCondition" : {"Hydrostatic":False, "Constant":True}
  }
}

# Do not modify what is below
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"surface_tension/3D_cube_to_sphere.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")