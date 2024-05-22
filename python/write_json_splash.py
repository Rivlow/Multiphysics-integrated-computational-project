{

  "domain":{
    "matrix_long" : [[0.01, 0.01, 0.01]],
    "matrix_orig" : [[0.01, 0.01, 0.01]],
    "vector_type" : [1],
    "L_d": [0.03, 0.03, 0.03],
    "o_d": [0.0, 0.0, 0.0]
  },

  "post_process":{
    "do": true,
    "xyz_init": [0, 0.015, 0.015],
    "xyz_end": [0.03,0.015, 0.015]
  },

  "simulation":{
    "theta" :0.5,
    "s": 0.00125,
    "nstepT": 500000,
    "dt": 0.0000001,
    "nsave": 100,
    "kappa": 2,
    "alpha": 0.01,
    "beta": 0,
    "alpha_st":1,
    "beta_adh": 10,
    "dimension": 3
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
    "gravity":false,
    "surface_tension":true,
    "adhesion":false
  },

  "condition":{
    "print_debug":false,
    "schemeIntegration": {"Euler":true,"RK22":false},
    "stateEquation" : {"Ideal gaz law":false, "Quasi incompresible fluid":true},
    "initialCondition" : {"Hydrostatic":false, "Constant":true}
  }
}