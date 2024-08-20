import json
import os
import sys
import numpy as np

s = 1e-2
L = 1e-1

dt = 1e-5
nsave = 500
nstepT = nsave*100

data = {
    
    "name_file" : "2D_no_adhesion",

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
        "alpha": 0.05,
        "beta": 0,
        "alpha_st": 1000,
        "beta_adh": 0.5,
        "dimension": 2,
        "scheme_integration": {"Euler": True, "RK22": False},
        "comparison_algorithm": False,
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
    
  "domain":{
    "matrix_long" : [np.round(np.array([L, s/2, L]), decimals = 4).tolist(), 
                     np.round(np.array([3*L, s/2, s/4]), decimals = 4).tolist(),
                     np.round(np.array([3*L-s,s/2, s/4]), decimals = 4).tolist()
                    ],
    "matrix_orig" : [np.round(np.array([L, 0, L + s/4]), decimals = 4).tolist(), 
                     np.round(np.array([0, 0, 2*L+s+s/4]), decimals = 4).tolist(),
                     np.round(np.array([s/2,0, 2*L+0.75*s]), decimals = 4).tolist()
                    ],

    "sphere": {
          "do": [0, 0, 0],
          "radius": [0.3]
      },

    "vector_type" : [1, 0, 0],
    "L_d": np.round(np.array([3*L, 4*s, 3*L]), decimals = 4).tolist(),
    "o_d": [0.0, 0.0, 0.0],

  },

  "post_process":{
    "do": False,
    "xyz_init": np.round(np.array([0, 0, 1.5*L]), decimals = 4).tolist(),
    "xyz_end": np.round(np.array([3*L, 0, 1.5*L]), decimals = 4).tolist()
  },


  "thermo":{
    "rho_0": 1000,
    "rho_moving": 1000,
    "rho_fixed" : 1000,
    "T": 298.15,
    "u_init": np.round(np.array([0.0, 0.0, 0.0]), decimals = 4).tolist(),
    "c_0": 15, 
    "gamma": 7, 
    "M": 18e-3, 
    "R":8.314,
    "sigma":0.0728
  },


  "forces":{
    "gravity":True,
    "surface_tension_1": False,
    "surface_tension_2": False,
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
json_src = f"adhesion/2D_no_adhesion.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")
