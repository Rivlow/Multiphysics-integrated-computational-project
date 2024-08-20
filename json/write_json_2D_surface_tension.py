import json
import os
import sys
import numpy as np

s = 5e-4
L = 0.01

dt = 1e-5
nsave = 1500
nstepT = nsave*250


data = {
    
    "name_file" : "2D_surface_tension",

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
        "alpha": 1.5,
        "beta": 0,
        "alpha_st": 1,
        "beta_adh": 1.2,
        "dimension": 2,
        "scheme_integration": {"Euler": True, "RK22": False},
        "comparison_algorithm": False,
    },

    "following_part": {
        "part": False,
        "min": True,
        "max": True,
        "particle": 500,
        "pressure": True,
        "rho": True,
        "position": [True, False, True],
        "velocity": [False, False, False],
    },
    

  "domain":{
    "matrix_long" : np.round(np.array([[L, s/2, L]]), decimals = 7).tolist(),
    "matrix_orig" : np.round(np.array([[L, 0, L]]), decimals = 7).tolist(),
    "sphere": {
                "do": [1],
                "radius": [L/2]
              },
    "vector_type" : [1],
    "L_d": np.round(np.array([3*L, 5*s, 3*L]), decimals = 7).tolist(),
    "o_d": np.round(np.array([0.0, 0.0, 0.0]), decimals = 7).tolist()
  },

  "post_process":{
    "do": True,
    "xyz_init": np.round(np.array([0, 0, 1.5*L]), decimals = 7).tolist(),
    "xyz_end": np.round(np.array([3*L, 0, 1.5*L]), decimals = 7).tolist()
  },


  "thermo":{
    "rho_0": 1000,
    "rho_moving": 1000,
    "rho_fixed" : 1000,
    "T": 298.15,
    "u_init": np.round(np.array([0.0, 0.0, 0.0]), decimals = 7).tolist(),
    "c_0": 15, 
    "gamma": 7, 
    "M": 18e-3, 
    "R":8.314,
    "sigma":0.072
  },


  "forces":{
    "gravity":False,
    "surface_tension_1": False,
    "surface_tension_2": True,
    "adhesion":False
  },

  "condition":{
    "print_debug":False,
    "state_equation" : {"Ideal gaz law":False, "Quasi incompresible fluid":True},
    "initial_condition" : {"Hydrostatic":False, "Constant":True}
  }
}

# Do not modify what is below
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"surface_tension/2D_cube_to_sphere.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")