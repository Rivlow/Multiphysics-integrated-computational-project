# How to use python scripts ?
Depending of the kind of simulation desired, chose the associated python script. You can directly modify the particle spacing **s** or the caracteristic lenght **L** while the other dimension will automatically be adapted.

Once selected, simply run the python script and the associated json file will be modified.

As mentionned in the main readme, you can next run the simulation in the console by the following coding line (here as an example, one simulate splash in 2D):
```solver.exe ../../json/splash/2D_splash.json```


# Simulation Models

This project utilizes various simulation models to study different fluid phenomena. Below is a list of the available models:

## 1) Splash

### a) Description
The Splash model simulates the fall of a single cube of particles into a box containing a floor, walls, and a ceiling. 

### b) Forces used 
**Only gravity and artificial viscosity** are considered.


## 2) Dam break

### a) Description
The dam break model simulates a rectangle of fluid particles (next to a wall) subjected to gravity with initial hydrostatic density assigned.

### b) Forces used 
**Only gravity and artificial viscosity** are considered.


## 3) Hydrostatic

### a) Description
The hydrostatic model simulates a rectangle of fluid particles subjected to gravity with initial constant densidy assigned. A hydrostatic pressure is expected at the end of the simulation.

### b) Forces used 
**Only gravity and artificial viscosity** are considered.


## 4) Surface Tension

### a) Description
The Surface Tension model simulates a single cube of particles in a **vacuum**. One expects the cube to turn into a sphere.

### b) Forces used
**Surface tension force and artificial viscosity** are considered.


## 5) Adhesion

### a) Description
A bunch of fluid particles is attached to a single ceiling with no floor. A fraction of the particles will experience adhesion force due to the fixed particles. One expects some particles to stick to the ceiling.

### b) Forces used
**Gravity, adhesion and artificial viscosity** are considered. 

## 5) Other

### a) Description
DO NOT USE THIS FOLDER !!! A bunch of fluid particles is present in a box of fixed particles. This file is only use to compar the efficiency of the linked list algorithm with respect to the naive neighbouring search. In order to make this file run, you have to disable all functions in the main.cpp file and enable the functions associated to naive_algo.

### b) Forces used
**No Gravity neither adhesion and artificial viscosity** considered. 