import json
import os
import sys

s = 0.00125
L = 1e-2
dimension = 3

nb_vtp_output = 250 # the total number of output file desired
dt = 0.0001
nstepT = 2500
nsave = nb_vtp_output/(dt*nstepT)  

data = {

  "domain":{
    "matrix_long" : [[L, L, L], 
                     [3*L, 3*L, s/2],
                     [3*L, 3*L, s/2]
                    ],
    "matrix_orig" : [[L, L, s], 
                     [0, 0, 0],
                     [s/2, s/2, s/2]
                    ],
    "vector_type" : [1, 0, 0],
    "L_d": [4*L, 4*L, 4*L],
    "o_d": [0.0, 0.0, 0.0]
  },

  "post_process":{
    "do": True,
    "xyz_init": [0, 1.5*L, 1.5*L],
    "xyz_end": [3*L, 1.5*L, 1.5*L]
  },

  "simulation":{
    "theta" :0.5,
    "s": s,
    "nstepT": nstepT,
    "dt": dt,
    "nsave": nsave,
    "kappa": 2,
    "alpha": 1500,
    "beta": 0,
    "alpha_st":0.15,
    "beta_adh": 10,
    "dimension": dimension
  },

  "thermo":{
    "rho_0": 1000,
    "rho_moving": 1000,
    "rho_fixed" : 1000,
    "T": 298.15,
    "u_init": [0.0, 0.0, 0.0],
    "c_0": 30, 
    "gamma": 7, 
    "M": 18e-3, 
    "R":8.314,
    "sigma":0.0728
  },


  "forces":{
    "gravity":True,
    "surface_tension":True,
    "adhesion":True
  },

  "condition":{
    "print_debug":False,
    "schemeIntegration": {"Euler":True,"RK22":False},
    "stateEquation" : {"Ideal gaz law":False, "Quasi incompresible fluid":True},
    "initialCondition" : {"Hydrostatic":False, "Constant":True}
  }
}

# Do not modify what is below
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
json_src = f"adhesion/{dimension}D_adhesion.json"

with open(f'{current_directory}/{json_src}', 'w') as json_file:
    json.dump(data, json_file, indent=4)


print(f"Data written in '{json_src}'")