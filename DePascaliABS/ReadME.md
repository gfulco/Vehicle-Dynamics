# ReadMe file for the simulation ABS with discrete dynamics
In this folder you can find the Matlab and Simulink files for the simulation of a braking maneuver with ABS control, considering actuators with discrete dynamics (Hydarulic Actuated Brake HAB). The project is written with Matlab and Simulink, version r2016b.

The folder contains:

1. `HAB_control.slx`: Simulink file containing the vehicle model (non-linear single corner model) and the ABS controller (state machine implementing the switching logic).
2. `burckhardt.m`: Matlab function implementing _Burckhardt_ friction model. It takes in input longitudinal wheel slip and the road condition and returns the value of the longitudinal friction coefficient. This model is an approximate form of the _Pacejka Magic Formula_. 
3. `mainABS.m`: Matlab script that defines general parameters for the Simulation, runs the Simulink model `ABS_contDyn_model.slx` and processes simulation results.

To simulate the braking maneuver, run `mainABS.m` script.



