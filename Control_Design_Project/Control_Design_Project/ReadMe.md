# ReadMe file for the project Vehicle Longitudinal and Lateral Dynamics Control design 
This file contains information about the files contained in the folder `Project\`. The vehicle model has been written using Simulink and __Matlab r2016b__.
Open `mainVehicleModel.m` and follow the instructions therein described. Execute the file to initialize parameters and run the Simulink model `FullVehicleModel_ABS_ESP.slx`.
The Matlab function `FullVehicleModel_V3_sfun.m` contains the Matlab S-function that models the dynamics of a full vehicle model (formula SAE size). The vehicle model implements also suitable adjustments to describe the slip dynamics at low speed `low_speed_slip.m`.

Tire-road contact forces are modeled with Pacejka magic formula contained in the Matlab function `pacejka_model.m`.