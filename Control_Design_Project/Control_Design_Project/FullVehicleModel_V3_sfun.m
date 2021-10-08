function FullVehicleModel_V3_sfun(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 5;
block.NumOutputPorts = 35;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts
  block.InputPort(i).Dimensions        = 1;
  block.InputPort(i).DatatypeID        = 0;  % double
  block.InputPort(i).Complexity        = 'Real';
  block.InputPort(i).DirectFeedthrough = false;
end

% states_name = {'u','omega_r','omega_f','Fx_t','Fx_r','lambda_r','lambda_r'};

% Override output port properties
for i = 1:block.NumOutputPorts
  block.OutputPort(i).Dimensions      = 1;
  block.OutputPort(i).DatatypeID      = 0; % double
  block.OutputPort(i).Complexity      = 'Real';
  block.OutputPort(i).SamplingMode    = 'Sample';
end

% SetInputPortSamplingMode:
%   Functionality    : Check and set input and output port
%                      attributes and specify whether the port is operating
%                      in sample-based or frame-based mode
%   C-Mex counterpart: mdlSetInputPortFrameData.
%   (The DSP System Toolbox is required to set a port as frame-based)
%
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);

% Register parameters
block.NumDialogPrms = 3;
%   Parameter 1 --> vehicle_data
%   Parameter 2 --> initial_conditions_data
%   Parameter 3 --> auxiliary_data

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Setup Number of Continuous states
block.NumContStates = 27;

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
% block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
% block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
% function DoPostPropSetup(block)
% block.NumDworks = 1;
%
%   block.Dwork(1).Name            = 'x1';
%   block.Dwork(1).Dimensions      = 1;
%   block.Dwork(1).DatatypeID      = 0;      % double
%   block.Dwork(1).Complexity      = 'Real'; % real
%   block.Dwork(1).UsedAsDiscState = true;


%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is
%%                      present in an enabled subsystem configured to reset
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)


%   ___ _        _
%  / __| |_ __ _| |_ ___ ___
%  \__ \  _/ _` |  _/ -_|_-<
%  |___/\__\__,_|\__\___/__/

% x(1)  = u(t)        -> Longitudinal vehicle velocity
% x(2)  = v(t)        -> Lateral vehicle velocity
% x(3)  = Omega(t)    -> Yaw rate
% x(4)  = phi(t)      -> Roll angle
% x(5)  = delta(t)    -> Steering angle
% x(6)  = alpha_rr(t) -> Side slip angle, rear right wheel
% x(7)  = alpha_rl(t) -> Side slip angle, rear left wheel
% x(8)  = alpha_fr(t) -> Side slip angle, front right wheel
% x(9)  = alpha_fl(t) -> Side slip angle, front left wheel
% x(10) = kappa_rr(t) -> Longitudinal slip, rear right wheel
% x(11) = kappa_rl(t) -> Longitudinal slip, rear left wheel
% x(12) = kappa_fr(t) -> Longitudinal slip, front right wheel
% x(13) = kappa_fl(t) -> Longitudinal slip, front left wheel
% x(14) = Fz_rr(t)    -> Vertical load, rear right wheel
% x(15) = Fz_rl(t)    -> Vertical load, rear rear left wheel
% x(16) = Fz_fr(t)    -> Vertical load, front right wheel
% x(17) = Fz_fl(t)    -> Vertical load, front left wheel
% x(18) = phi_dot(t)  -> Roll rate
% x(19) = omega_rr(t) -> angular velocity rear right wheel
% x(20) = omega_rl(t) -> angular velocity rear left wheel
% x(21) = omega_fr(t) -> angular velocity front right wheel
% x(22) = omega_fl(t) -> angular velocity front left wheel

% Copy Parameters from Dialog Box
vehicle_data            = block.DialogPrm(1).Data;
initial_conditions_data = block.DialogPrm(2).Data;
auxiliary_data          = block.DialogPrm(3).Data;

g     = auxiliary_data.environment.g;
Lf    = vehicle_data.vehicle.Lf;
Lr    = vehicle_data.vehicle.Lr;
L     = Lf+Lr;
tau_D = vehicle_data.steering_system.tau_D;
Rf    = vehicle_data.front_wheel.Rf; % Front Wheel Radius
Rr    = vehicle_data.rear_wheel.Rr; % Rear Wheel Radius
m     = vehicle_data.vehicle.m; % Vehicle Mass
Wf    = vehicle_data.vehicle.Wf;    % [m] Width of front wheels axle 
Wr    = vehicle_data.vehicle.Wr;    % [m] Width of rear wheels axle 

% Aerodynamics
CAx   = vehicle_data.aerodynamics.CAx;
CAzf  = vehicle_data.aerodynamics.CAzf;
CAzr  = vehicle_data.aerodynamics.CAzr;

% Initial conditions
u0       = initial_conditions_data.longSpeed;
v0       = initial_conditions_data.latSpeed;
Omega0   = initial_conditions_data.Omega;
phi0     = initial_conditions_data.Phi;
delta0   = initial_conditions_data.delta / tau_D;
phi_dot0 = initial_conditions_data.Phi_dot;

% Aerodynamics Forces
FAz_r0 = +CAzr*u0^2;
FAz_f0 = +CAzf*u0^2;

% Longitudinal Wheel Speed
vx_rr0 = u0 - Omega0*Wr/2;
vx_rl0 = u0 + Omega0*Wr/2;
vx_fr0 = u0 - Omega0*Wf/2;
vx_fl0 = u0 + Omega0*Wf/2;

% Lateral Wheel Speed
vy_rr0 = v0 - Omega0*Lr;
vy_rl0 = v0 - Omega0*Lr;
vy_fr0 = v0 + Omega0*Lf;
vy_fl0 = v0 + Omega0*Lf;

% Wheel Side Slip Angle
alpha_rr0 = atan(vy_rr0/vx_rr0);
alpha_rl0 = atan(vy_rl0/vx_rl0);
alpha_fr0 = atan(vy_fr0/vx_fr0)-delta0;
alpha_fl0 = atan(vy_fl0/vx_fl0)-delta0;

% Wheel angular speed
omega_rr0 = vx_rr0/Rr;
omega_rl0 = vx_rl0/Rr;
omega_fr0 = vx_fr0/Rf;
omega_fl0 = vx_fl0/Rf;

% Wheel longitudinal speed (Condition of no longitudinal speed)
kappa_rr0 = 0;
kappa_rl0 = 0;
kappa_fr0 = 0;
kappa_fl0 = 0;

% Vertical load
Fz_rr0 = (m * g * Lf + FAz_r0 * L) / L / 0.2e1;
Fz_rl0 = (m * g * Lf + FAz_r0 * L) / L / 0.2e1;
Fz_fr0 = (Lr * g * m + FAz_f0 * L) / L / 0.2e1;
Fz_fl0 = (Lr * g * m + FAz_f0 * L) / L / 0.2e1;

u = zeros(block.NumInputPorts,1);
for i = 1:block.NumInputPorts
  u(i) = block.InputPort(i).Data;
end

block.ContStates.Data = [...
  u0; v0; Omega0; phi0 ; delta0 ; ...
  alpha_rr0 ; alpha_rl0 ; alpha_fr0 ; alpha_fl0 ;...
  kappa_rr0 ; kappa_rl0 ; kappa_fr0 ; kappa_fl0 ; ...
  Fz_rr0 ; Fz_rl0 ; Fz_fr0 ; Fz_fl0 ; ...
  phi_dot0; ...
  omega_rr0 ; omega_rl0 ; omega_fr0 ; omega_fl0; u];
%end InitializeConditions


%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this
%%                      is the place to do it.
%%   Required         : No
%%   C-MEX counterpart: mdlStart
%%
% function Start(block)
%
% block.Dwork(1).Data = 0;

%end Start

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)

x = block.ContStates.Data;

for i = 1:23
  block.OutputPort(i).Data = x(i);
end

for i=24:27
  block.OutputPort(i).Data = x(i)*regSign(x(i-5));
end

% Copy Parameters from Dialog Box
vehicle_data            = block.DialogPrm(1).Data;
auxiliary_data          = block.DialogPrm(3).Data;


g     = auxiliary_data.environment.g;
Lf    = vehicle_data.vehicle.Lf;
Lr    = vehicle_data.vehicle.Lr;
L     = Lf+Lr;
tau_D = vehicle_data.steering_system.tau_D;
Rf    = vehicle_data.front_wheel.Rf; % Front Wheel Radius
Rr    = vehicle_data.rear_wheel.Rr; % Rear Wheel Radius
m     = vehicle_data.vehicle.m; % Vehicle Mass
Wf    = vehicle_data.vehicle.Wf;    % [m] Width of front wheels axle 
Wr    = vehicle_data.vehicle.Wr;    % [m] Width of rear wheels axle 

%  _____               ___       _        
% |_   _|  _ _ _ ___  |   \ __ _| |_ __ _ 
%   | || || | '_/ -_) | |) / _` |  _/ _` |
%   |_| \_, |_| \___| |___/\__,_|\__\__,_|
%       |__/                              

lx_rear            = vehicle_data.rear_tyre.long_relax_length;
ly_rear            = vehicle_data.rear_tyre.lat_relax_length;
tau_N_rear         = vehicle_data.rear_tyre.tau_N;
Vth_rear           = vehicle_data.rear_tyre.Vth; % [m/s] Threshod for low speed pacejka model
pacejkaParam_rear  = vehicle_data.rear_tyre.pacejkaParam;

lx_front           = vehicle_data.front_tyre.long_relax_length;
ly_front           = vehicle_data.front_tyre.lat_relax_length;
tau_N_front        = vehicle_data.front_tyre.tau_N;
Vth_front          = vehicle_data.front_tyre.Vth; % [m/s] Threshod for low speed pacejka model
pacejkaParam_front = vehicle_data.front_tyre.pacejkaParam;

% Camber angle
grr_cc = vehicle_data.rear_suspension.camber_right_wheel; % [-0.53447e-1, +0.839386,    -0.201305,    -0.87213e-1];
grl_cc = vehicle_data.rear_suspension.camber_left_wheel ;%  [ 0.53447e-1, +0.839386,     0.201305,    -0.87213e-1];
gfr_cc = vehicle_data.front_suspension.camber_right_wheel; % [-0.53447e-1, +0.839386,    -0.201305,    -0.87213e-1];
gfl_cc = vehicle_data.front_suspension.camber_left_wheel ;%  [ 0.53447e-1, +0.839386,     0.201305,    -0.87213e-1];

phi  = 0*x(4); %  actual roll angle value
gamma_rr = grr_cc(1)+grr_cc(2)*phi+grr_cc(3)*phi^2+grr_cc(4)*phi^3;
gamma_rl = grl_cc(1)+grl_cc(2)*phi+grl_cc(3)*phi^2+grl_cc(4)*phi^3;
gamma_fr = gfr_cc(1)+gfr_cc(2)*phi+gfr_cc(3)*phi^2+gfr_cc(4)*phi^3;
gamma_fl = gfl_cc(1)+gfl_cc(2)*phi+gfl_cc(3)*phi^2+gfl_cc(4)*phi^3;


% Steering system model
delta_fl = x(5); % + missing roll steer correction
delta_fr = x(5); % + missing roll steer correction

% Fy : Lateral Force (N)
% Fx : Longitudinal Force (N)
% Mz : Aligning Torque (N-m)


% The followig part aims at damping the oscillations arising at low speed
% See Chapter 8.6 Pacejka
Vx_rr = ( (x(3) * Wr) / 0.2e1 + x(1));
Vx_rl = (-(x(3) * Wr) / 0.2e1 + x(1));
Vx_fr = ( (x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf));
Vx_fl = (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf));

Vsx_rr = ((x(19) * Rr) - ( (x(3) * Wr) / 0.2e1 + x(1)));
Vsx_rl = ((x(20) * Rr) - (-(x(3) * Wr) / 0.2e1 + x(1)));
Vsx_fr = ((x(21) * Rf) - ( (x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf)));
Vsx_fl = ((x(22) * Rf) - (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf)));

Vsy_rr = 2 * x(3) * Lr - 2 * x(2);
Vsy_rl = 2 * x(3) * Lr - 2 * x(2);
Vsy_fr =  delta_fr * x(3) * Wf - 2 * x(3) * Lf + 2 * delta_fr * x(1) - 2 * x(2);
Vsy_fl = -delta_fl * x(3) * Wf - 2 * x(3) * Lf + 2 * delta_fl * x(1) - 2 * x(2);

[kappa_eff_rr, alpha_eff_rr] = low_speed_slip(Vx_rr,Vsx_rr,Vsy_rr,x(10),x(6),x(14),gamma_rr,Vth_rear,pacejkaParam_rear);
[kappa_eff_rl, alpha_eff_rl] = low_speed_slip(Vx_rl,Vsx_rl,Vsy_rl,x(11),x(7),x(15),gamma_rl,Vth_rear,pacejkaParam_rear);
[kappa_eff_fr, alpha_eff_fr] = low_speed_slip(Vx_fr,Vsx_fr,Vsy_fr,x(12),x(8),x(16),gamma_fr,Vth_front,pacejkaParam_front);
[kappa_eff_fl, alpha_eff_fl] = low_speed_slip(Vx_fl,Vsx_fl,Vsy_fl,x(13),x(9),x(17),gamma_fl,Vth_front,pacejkaParam_front);

% Rear Right Wheel
[Fx_rr,Fy_rr,~] = pacejka_model(kappa_eff_rr,x(14),alpha_eff_rr,gamma_rr,pacejkaParam_rear);
% Fx_rr=0;
% Fy_rr=0;
% AM_rr=0;

% Rear Left Wheel
[Fx_rl,Fy_rl,~] = pacejka_model(kappa_eff_rl,x(15),alpha_eff_rl,gamma_rl,pacejkaParam_rear);
% Fx_rl=0;
% Fy_rl=0;
% AM_rl=0;

% Front Right Wheel
[Fx_fr,Fy_fr,~] = pacejka_model(kappa_eff_fr,x(16),alpha_eff_fr,gamma_fr,pacejkaParam_front);
% Fx_fr=0;
% Fy_fr=0;
% AM_fr=0;

% Front Left Wheel
[Fx_fl,Fy_fl,~] = pacejka_model(kappa_eff_fl,x(17),alpha_eff_fl,gamma_fl,pacejkaParam_front);

block.OutputPort(28).Data = Fx_rr;
block.OutputPort(29).Data = Fy_rr;
block.OutputPort(30).Data = Fx_rl;
block.OutputPort(31).Data = Fy_rl;
block.OutputPort(32).Data = Fx_fr;
block.OutputPort(33).Data = Fy_fr;
block.OutputPort(34).Data = Fx_fl;
block.OutputPort(35).Data = Fy_fl;

%end Outputs

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlUpdate
%%
% function Update(block)
%
% block.Dwork(1).Data = block.InputPort(1).Data;

%end Update

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C-MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

% System states
% x(1)  = u(t)        -> Longitudinal vehicle velocity
% x(2)  = v(t)        -> Lateral vehicle velocity
% x(3)  = Omega(t)    -> Yaw rate
% x(4)  = phi(t)      -> Roll angle
% x(5)  = delta(t)    -> Steering angle
% x(6)  = alpha_rr(t) -> Side slip angle, rear right wheel
% x(7)  = alpha_rl(t) -> Side slip angle, rear left wheel
% x(8)  = alpha_fr(t) -> Side slip angle, front right wheel
% x(9)  = alpha_fl(t) -> Side slip angle, front left wheel
% x(10) = kappa_rr(t) -> Longitudinal slip, rear right wheel
% x(11) = kappa_rl(t) -> Longitudinal slip, rear left wheel
% x(12) = kappa_fr(t) -> Longitudinal slip, front right wheel
% x(13) = kappa_fl(t) -> Longitudinal slip, front left wheel
% x(14) = Fz_rr(t)    -> Vertical load, rear right wheel
% x(15) = Fz_rl(t)    -> Vertical load, rear rear left wheel
% x(16) = Fz_fr(t)    -> Vertical load, front right wheel
% x(17) = Fz_fl(t)    -> Vertical load, front left wheel
% x(18) = phi_dot(t)  -> Roll rate
% x(19) = omega_rr(t) -> angular velocity rear right wheel
% x(20) = omega_rl(t) -> angular velocity rear left wheel
% x(21) = omega_fr(t) -> angular velocity front right wheel
% x(22) = omega_fl(t) -> angular velocity front left wheel

% Inputs
% u(1) = delta_D -> Desired steering angle
% u(2) = Tw_rr   -> Torque at rear right wheel
% u(3) = Tw_rl   -> Torque at rear left wheel
% u(4) = Tw_fr   -> Torque at front right wheel
% u(5) = Tw_fl   -> Torque at front left wheel

% Initialize input and output
x = block.ContStates.Data;
u = zeros(block.NumInputPorts,1);
for i = 1:block.NumInputPorts
  u(i) = block.InputPort(i).Data;
end

% Copy Parameters from Dialog Box
vehicle_data            = block.DialogPrm(1).Data;
auxiliary_data          = block.DialogPrm(3).Data;

%   ___       _
%  |   \ __ _| |_ __ _
%  | |) / _` |  _/ _` |
%  |___/\__,_|\__\__,_|
%

g     = auxiliary_data.environment.g;
Lf    = vehicle_data.vehicle.Lf;
Lr    = vehicle_data.vehicle.Lr;
L     = Lf+Lr;
tau_D = vehicle_data.steering_system.tau_D;
Rf    = vehicle_data.front_wheel.Rf; % Front Wheel Radius
Rr    = vehicle_data.rear_wheel.Rr; % Rear Wheel Radius
m     = vehicle_data.vehicle.m; % Vehicle Mass
Wf    = vehicle_data.vehicle.Wf;    % [m] Width of front wheels axle 
Wr    = vehicle_data.vehicle.Wr;    % [m] Width of rear wheels axle 

hs    = vehicle_data.vehicle.hs;
hr    = vehicle_data.vehicle.hr;
hrr   = vehicle_data.rear_suspension.hrr;
Ks_r  = vehicle_data.rear_suspension.Ks_r; % Rear suspension stiffness
Cs_r  = vehicle_data.rear_suspension.Cs_r; % Rear suspension damping
Ks_f  = vehicle_data.front_suspension.Ks_f; % Front suspension stiffness
Cs_f  = vehicle_data.front_suspension.Cs_f; % Front suspension dumping
hrf   = vehicle_data.front_suspension.hrf;

tau_H = vehicle_data.steering_system.tau_H;
% Steering system model
delta_fl = x(5);
delta_fr = x(5);

%  _____               ___       _        
% |_   _|  _ _ _ ___  |   \ __ _| |_ __ _ 
%   | || || | '_/ -_) | |) / _` |  _/ _` |
%   |_| \_, |_| \___| |___/\__,_|\__\__,_|
%       |__/                              

lx_rear            = vehicle_data.rear_tyre.long_relax_length;
ly_rear            = vehicle_data.rear_tyre.lat_relax_length;
tau_N_rear         = vehicle_data.rear_tyre.tau_N;
Vth_rear           = vehicle_data.rear_tyre.Vth; % [m/s] Threshod for low speed pacejka model
pacejkaParam_rear  = vehicle_data.rear_tyre.pacejkaParam;

lx_front           = vehicle_data.front_tyre.long_relax_length;
ly_front           = vehicle_data.front_tyre.lat_relax_length;
tau_N_front        = vehicle_data.front_tyre.tau_N;
Vth_front          = vehicle_data.front_tyre.Vth; % [m/s] Threshod for low speed pacejka model
pacejkaParam_front = vehicle_data.front_tyre.pacejkaParam;

% Camber angle
% Camber angle
grr_cc = vehicle_data.rear_suspension.camber_right_wheel; % [-0.53447e-1, +0.839386,    -0.201305,    -0.87213e-1];
grl_cc = vehicle_data.rear_suspension.camber_left_wheel ;%  [ 0.53447e-1, +0.839386,     0.201305,    -0.87213e-1];
gfr_cc = vehicle_data.front_suspension.camber_right_wheel; % [-0.53447e-1, +0.839386,    -0.201305,    -0.87213e-1];
gfl_cc = vehicle_data.front_suspension.camber_left_wheel ;%  [ 0.53447e-1, +0.839386,     0.201305,    -0.87213e-1];

phi  = 0*x(4); %  actual roll angle value
gamma_rr = grr_cc(1)+grr_cc(2)*phi+grr_cc(3)*phi^2+grr_cc(4)*phi^3;
gamma_rl = grl_cc(1)+grl_cc(2)*phi+grl_cc(3)*phi^2+grl_cc(4)*phi^3;
gamma_fr = gfr_cc(1)+gfr_cc(2)*phi+gfr_cc(3)*phi^2+gfr_cc(4)*phi^3;
gamma_fl = gfl_cc(1)+gfl_cc(2)*phi+gfl_cc(3)*phi^2+gfl_cc(4)*phi^3;

%   __  __                           _    ___      ___
%  |  \/  |__ _ ______  __ _ _ _  __| |  / __|___ / __|
%  | |\/| / _` (_-<_-< / _` | ' \/ _` | | (__/ _ \ (_ |
%  |_|  |_\__,_/__/__/ \__,_|_||_\__,_|  \___\___/\___|

% Wheel mass
m_wr = vehicle_data.rear_wheel.mass;
m_wf = vehicle_data.front_wheel.mass;
% Front unsprung mass
m_uf = vehicle_data.front_unsprung.mass;
% Rear unsprung mass
m_ur = vehicle_data.rear_unsprung.mass;
% Sprung Mass
ms = m - m_uf - m_ur;

% Position of the center of mass
h_G = (hs + hr) / m * ms + (m_uf * Rf + m_ur * Rr) / m;


%   ___              _   _
%  |_ _|_ _  ___ _ _| |_(_)__ _
%   | || ' \/ -_) '_|  _| / _` |
%  |___|_||_\___|_|  \__|_\__,_|

% VEHICLE
i_xx = vehicle_data.vehicle.i_xx;
i_yy = vehicle_data.vehicle.i_yy;
i_zz = vehicle_data.vehicle.i_zz;
I_xz = vehicle_data.vehicle.I_xz;

% CHASSIS
% is =  |  is_xx   0   -is_xz |
%       |    0   is_yy    0   |
%       | -is_xz   0    is_zz |
is_xx = vehicle_data.chassis.is_xx;
is_yy = vehicle_data.chassis.is_yy;
is_zz = vehicle_data.chassis.is_zz;
is_xz = vehicle_data.chassis.is_xz;

% Moment of inertia about X axis of RF0
IXXs0 = hs^2*ms+is_xx;             % chassis moment of inertia w.r.t. roll axis
IYZs0 = is_zz -hs^2*ms-is_yy;      % chassis moment of inertia w.r.t. roll axis
IXX0  = (hr+hs)*hs*ms+is_xx;       % chassis moment of inertia w.r.t. RF0
IYZ0  = ms*(hr+hs)*hs+is_yy-is_zz; % chassis moment of inertia w.r.t. RF0


% WHEEL
% iwd = | iwd   0  0  |
%       |  0  iwa  0  |
%       |  0   0  iwd |
w_wf = vehicle_data.front_wheel.width; % Wheel width
w_wr = vehicle_data.rear_wheel.width; % Wheel width
% Front Wheels
iwd_f = m_wf/12 * (3*Rf^2 + w_wf^2);
iwa_f = m_wf/2 * Rf^2;
% Rear Wheels
iwd_r = m_wr/12 * (3*Rr^2 + w_wr^2);
iwa_r = m_wr/2 * Rr^2;


%     _                   _                      _
%    /_\  ___ _ _ ___  __| |_  _ _ _  __ _ _ __ (_)__ ___
%   / _ \/ -_) '_/ _ \/ _` | || | ' \/ _` | '  \| / _(_-<
%  /_/ \_\___|_| \___/\__,_|\_, |_||_\__,_|_|_|_|_\__/__/
%                           |__/

CAx   = vehicle_data.aerodynamics.CAx;
CAzf  = vehicle_data.aerodynamics.CAzf;
CAzr  = vehicle_data.aerodynamics.CAzr;
FAxc  = +1/2*1.1*CAx*x(1)^2;
FAz_r = +CAzr*x(1)^2;
FAz_f = +CAzf*x(1)^2;

%% Lateral acceleration
ay = x(3)*x(1);

%% Suspension forces
Ms_f = -Ks_f * x(4) - Cs_r*x(18);
Ms_r = -Ks_r * x(4) - Cs_f*x(18);


%   ___              _ _          __  __         _     _
%  | _ \__ _ __ ___ (_) |____ _  |  \/  |___  __| |___| |
%  |  _/ _` / _/ -_)| | / / _` | | |\/| / _ \/ _` / -_) |
%  |_| \__,_\__\___|/ |_\_\__,_| |_|  |_\___/\__,_\___|_|
%                 |__/

% Fy : Lateral Force (N)
% Fx : Longitudinal Force (N)
% Mz : Aligning Torque (N-m)


% The followig part aims at damping the oscillations arising at low speed
% See Chapter 8.6 Pacejka
Vx_rr = ( (x(3) * Wr) / 0.2e1 + x(1));
Vx_rl = (-(x(3) * Wr) / 0.2e1 + x(1));
Vx_fr = ( (x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf));
Vx_fl = (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf));

Vsx_rr = ((x(19) * Rr) - ( (x(3) * Wr) / 0.2e1 + x(1)));
Vsx_rl = ((x(20) * Rr) - (-(x(3) * Wr) / 0.2e1 + x(1)));
Vsx_fr = ((x(21) * Rf) - ( (x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf)));
Vsx_fl = ((x(22) * Rf) - (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf)));

Vsy_rr = 2 * x(3) * Lr - 2 * x(2);
Vsy_rl = 2 * x(3) * Lr - 2 * x(2);
Vsy_fr =  delta_fr * x(3) * Wf - 2 * x(3) * Lf + 2 * delta_fr * x(1) - 2 * x(2);
Vsy_fl = -delta_fl * x(3) * Wf - 2 * x(3) * Lf + 2 * delta_fl * x(1) - 2 * x(2);

[kappa_eff_rr, alpha_eff_rr] = low_speed_slip(Vx_rr,Vsx_rr,Vsy_rr,x(10),x(6),x(14),gamma_rr,Vth_rear,pacejkaParam_rear);
[kappa_eff_rl, alpha_eff_rl] = low_speed_slip(Vx_rl,Vsx_rl,Vsy_rl,x(11),x(7),x(15),gamma_rl,Vth_rear,pacejkaParam_rear);
[kappa_eff_fr, alpha_eff_fr] = low_speed_slip(Vx_fr,Vsx_fr,Vsy_fr,x(12),x(8),x(16),gamma_fr,Vth_front,pacejkaParam_front);
[kappa_eff_fl, alpha_eff_fl] = low_speed_slip(Vx_fl,Vsx_fl,Vsy_fl,x(13),x(9),x(17),gamma_fl,Vth_front,pacejkaParam_front);


% Rear Right Wheel
[Fx_rr,Fy_rr,AM_rr] = pacejka_model(kappa_eff_rr,x(14),alpha_eff_rr,gamma_rr,pacejkaParam_rear);

% Rear Left Wheel
[Fx_rl,Fy_rl,AM_rl] = pacejka_model(kappa_eff_rl,x(15),alpha_eff_rl,gamma_rl,pacejkaParam_rear);

% Front Right Wheel
[Fx_fr,Fy_fr,AM_fr] = pacejka_model(kappa_eff_fr,x(16),alpha_eff_fr,gamma_fr,pacejkaParam_front);

% Front Left Wheel
[Fx_fl,Fy_fl,AM_fl] = pacejka_model(kappa_eff_fl,x(17),alpha_eff_fl,gamma_fl,pacejkaParam_front);

% Tota lself-aligning moment
AM = AM_rl + AM_rr + AM_fl + AM_fr;


%% Equations
% State x(1) --> Longitudinal speed: u(t)
dotx(1) = ((x(4) * hs * ms * Wr * Fx_rr * IXXs0 - x(4) * hs * ms * Wf * Fx_fl * IXXs0 + 4 * IXXs0 * x(3) * x(18) * hs * i_zz * ms - 2 * x(4) ^ 2 * hs * ms * IYZs0 * x(3) ^ 2 * is_xz + 2 * x(4) * hs * ms * Lf * Fy_fr * IXXs0 + 2 * x(4) * hs * ms * Lf * Fy_fl * IXXs0 - 2 * x(4) * hs * ms * Lr * Fy_rl * IXXs0 - 2 * x(4) * hs * ms * Lr * Fy_rr * IXXs0 - x(4) * hs * ms * Wr * Fx_rl * IXXs0 + x(4) * hs * ms * Wf * Fx_fr * IXXs0 + 2 * x(4) * hs ^ 2 * ms ^ 2 * ay * is_xz + 2 * x(4) * hs * ms * AM_rl * IXXs0 + 2 * x(4) * hs * ms * AM_rr * IXXs0 + 2 * x(4) * hs * ms * AM_fl * IXXs0 + 2 * x(4) * hs * ms * AM_fr * IXXs0 + 2 * x(4) * hs * ms * Ms_f * is_xz + 2 * x(4) * hs * ms * Ms_r * is_xz + 2 * x(4) ^ 2 * hs ^ 2 * ms ^ 2 * g * is_xz - 4 * x(3) * x(18) * hs * is_xz ^ 2 * ms - 2 * IXXs0 * x(3) * x(2) * i_zz * m + 2 * x(3) * x(2) * is_xz ^ 2 * m + 2 * Fx_rl * is_xz ^ 2 - 2 * FAxc * is_xz ^ 2 + 2 * Fx_rr * is_xz ^ 2 + 2 * Fx_fl * is_xz ^ 2 + 2 * Fx_fr * is_xz ^ 2 - 2 * delta_fl * Fy_fl * is_xz ^ 2 - 2 * delta_fr * Fy_fr * is_xz ^ 2 - 2 * IXXs0 * Fx_rl * i_zz + 2 * IXXs0 * FAxc * i_zz - 2 * IXXs0 * Fx_rr * i_zz - 2 * IXXs0 * Fx_fl * i_zz - 2 * IXXs0 * Fx_fr * i_zz + 2 * IXXs0 * delta_fl * Fy_fl * i_zz + 2 * IXXs0 * delta_fr * Fy_fr * i_zz + 2 * x(4) * hs * ms * delta_fl * Lf * Fx_fl * IXXs0 + 2 * x(4) * hs * ms * delta_fr * Lf * Fx_fr * IXXs0 + x(4) * hs * ms * delta_fl * Wf * Fy_fl * IXXs0 - x(4) * hs * ms * delta_fr * Wf * Fy_fr * IXXs0) / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m)) / 0.2e1;

% State x(2) --> Lateral speed: v(t)
dotx(2) = -((2 * IXXs0 * delta_fl * Fx_fl * i_zz * m + 2 * IXXs0 * delta_fr * Fx_fr * i_zz * m + 2 * hs * ms * is_xz * AM_fr * m - 2 * x(4) ^ 3 * g * hs ^ 4 * ms ^ 4 - 2 * ay * x(4) ^ 2 * hs ^ 4 * ms ^ 4 - 2 * x(4) ^ 2 * Ms_f * hs ^ 3 * ms ^ 3 - 2 * x(4) ^ 2 * Ms_r * hs ^ 3 * ms ^ 3 + 2 * x(1) * x(3) * is_xz ^ 2 * m ^ 2 + 2 * IXXs0 * Fy_fr * i_zz * m + 2 * IXXs0 * Fy_rr * i_zz * m + 2 * IXXs0 * Fy_fl * i_zz * m - 2 * delta_fl * Fx_fl * is_xz ^ 2 * m - 2 * delta_fr * Fx_fr * is_xz ^ 2 * m + 2 * IXXs0 * Fy_rl * i_zz * m - 2 * hs ^ 2 * ms ^ 2 * is_xz * Fx_rl * x(4) - 2 * hs ^ 2 * ms ^ 2 * is_xz * Fx_fr * x(4) + 2 * hs ^ 2 * ms ^ 2 * is_xz * FAxc * x(4) - 2 * hs ^ 2 * ms ^ 2 * is_xz * x(4) * Fx_rr - 2 * hs ^ 2 * ms ^ 2 * is_xz * x(4) * Fx_fl + 2 * IYZs0 * x(3) ^ 2 * x(4) ^ 3 * hs ^ 3 * ms ^ 3 + 2 * IXXs0 * x(3) ^ 2 * x(4) ^ 3 * hs ^ 3 * ms ^ 3 - 2 * IXXs0 * Fy_rl * x(4) ^ 2 * hs ^ 2 * ms ^ 2 - 2 * IXXs0 * Fy_fr * x(4) ^ 2 * hs ^ 2 * ms ^ 2 - 2 * IXXs0 * Fy_rr * x(4) ^ 2 * hs ^ 2 * ms ^ 2 - 2 * IXXs0 * x(4) ^ 2 * Fy_fl * hs ^ 2 * ms ^ 2 - 2 * IXXs0 * x(1) * x(3) * i_zz * m ^ 2 + 2 * hs * ms * is_xz * AM_rl * m + 2 * hs * ms * is_xz * AM_rr * m + 2 * hs * ms * is_xz * AM_fl * m + 2 * ay * hs ^ 2 * i_zz * m * ms ^ 2 + 2 * Ms_f * hs * i_zz * m * ms + 2 * Ms_r * hs * i_zz * m * ms - 2 * IXXs0 * delta_fr * Fx_fr * x(4) ^ 2 * hs ^ 2 * ms ^ 2 + 4 * hs ^ 3 * ms ^ 3 * is_xz * x(3) * x(4) * x(18) + 2 * hs ^ 2 * ms ^ 2 * is_xz * delta_fl * x(4) * Fy_fl + 2 * hs ^ 2 * ms ^ 2 * is_xz * delta_fr * Fy_fr * x(4) + 2 * x(3) ^ 2 * x(4) * hs * is_xz ^ 2 * m * ms + 2 * x(4) * g * hs ^ 2 * i_zz * m * ms ^ 2 - 2 * IXXs0 * delta_fl * x(4) ^ 2 * Fx_fl * hs ^ 2 * ms ^ 2 - 2 * hs ^ 2 * ms ^ 2 * is_xz * x(3) * x(2) * x(4) * m - 2 * IYZs0 * x(3) ^ 2 * x(4) * hs * i_zz * m * ms + 2 * IXXs0 * x(1) * x(3) * x(4) ^ 2 * hs ^ 2 * m * ms ^ 2 - 2 * IXXs0 * x(3) ^ 2 * x(4) * hs * i_zz * m * ms + 2 * hs * ms * is_xz * delta_fl * Lf * Fx_fl * m + 2 * hs * ms * is_xz * delta_fr * Lf * Fx_fr * m + hs * ms * is_xz * delta_fl * Wf * Fy_fl * m - hs * ms * is_xz * delta_fr * Wf * Fy_fr * m - 2 * Fy_rl * is_xz ^ 2 * m - 2 * Fy_fr * is_xz ^ 2 * m - 2 * Fy_rr * is_xz ^ 2 * m - 2 * Fy_fl * is_xz ^ 2 * m + 2 * hs * ms * is_xz * Lf * Fy_fl * m - 2 * hs * ms * is_xz * Lr * Fy_rl * m - 2 * hs * ms * is_xz * Lr * Fy_rr * m - hs * ms * is_xz * Wr * Fx_rl * m + hs * ms * is_xz * Wf * Fx_fr * m + hs * ms * is_xz * Wr * Fx_rr * m - hs * ms * is_xz * Wf * Fx_fl * m + 2 * hs * ms * is_xz * Lf * Fy_fr * m) / m / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m)) / 0.2e1;

% State x(3) --> Yaw rate: Omega(t)
dotx(3) = -(1 / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m) * (2 * x(4) * g * hs * is_xz * ms * m + 4 * IXXs0 * x(3) * x(18) * x(4) * hs ^ 2 * ms ^ 2 + 2 * IXXs0 * delta_fl * x(4) * Fy_fl * hs * ms + 2 * IXXs0 * delta_fr * Fy_fr * x(4) * hs * ms - 2 * IYZs0 * x(3) ^ 2 * x(4) * is_xz * m - 2 * IXXs0 * Fx_rl * x(4) * hs * ms - 2 * IXXs0 * Fx_fr * x(4) * hs * ms + 2 * IXXs0 * FAxc * x(4) * hs * ms - 2 * IXXs0 * Fx_rr * x(4) * hs * ms - 2 * IXXs0 * x(4) * Fx_fl * hs * ms + 2 * delta_fl * Lf * Fx_fl * IXXs0 * m + 2 * delta_fr * Lf * Fx_fr * IXXs0 * m + delta_fl * Wf * Fy_fl * IXXs0 * m - delta_fr * Wf * Fy_fr * IXXs0 * m + 2 * ay * hs * is_xz * ms * m + 2 * Lf * Fy_fr * IXXs0 * m + 2 * Lf * Fy_fl * IXXs0 * m - 2 * Lr * Fy_rl * IXXs0 * m - 2 * Lr * Fy_rr * IXXs0 * m - Wr * Fx_rl * IXXs0 * m + Wf * Fx_fr * IXXs0 * m + Wr * Fx_rr * IXXs0 * m - Wf * Fx_fl * IXXs0 * m + 2 * AM_rr * IXXs0 * m + 2 * AM_fl * IXXs0 * m + 2 * AM_fr * IXXs0 * m + 2 * Ms_f * is_xz * m + 2 * Ms_r * is_xz * m + 2 * AM_rl * IXXs0 * m - 2 * IXXs0 * x(3) * x(2) * x(4) * hs * m * ms)) / 0.2e1;

% State x(4) --> Roll angle: phi(t)
dotx(4) = x(18);

% State x(5) --> Steering angle: delta(t)
dotx(5) = -0.1e1 / tau_D / tau_H * (x(5) * tau_D - x(23));

% State x(6) --> Side slip angle, rear right wheel: alpha_rr(t)
if (x(3) * Wr / 0.2e1 + x(1)) <= Vth_rear
  x6_ss = 2*( -(-x(3) * Lr + x(2))) / (Vth_rear + (x(3) * Wr / 0.2e1 + x(1))^2 / Vth_rear);
else
  x6_ss = -(-x(3) * Lr + x(2)) / (x(3) * Wr / 0.2e1 + x(1));
end
dotx(6) = -abs(x(1))*(x(6) - x6_ss)/ly_rear;

% State x(7) --> Side slip angle, rear left wheel: alpha_rl(t)
if (x(3) * Wr / 0.2e1 + x(1)) <= Vth_rear
  x7_ss = 2*(-(-x(3) * Lr + x(2))) / (Vth_rear + (x(3) * Wr / 0.2e1 + x(1))^2 / Vth_rear);
else
  x7_ss = -(-x(3) * Lr + x(2)) / (x(3) * Wr / 0.2e1 + x(1));
end
dotx(7) = -abs(x(1))*(x(7) - x7_ss)/ly_rear;

% State x(8) --> Side slip angle, front right wheel: alpha_fr(t)
if (x(3) * Wf / 0.2e1 + x(1) + delta_fr * x(3) * Lf) <= Vth_front
  x8_ss = 2*(-(x(3) * Lf + x(2) + (-x(3) * Wf / 0.2e1 - x(1)) * delta_fr)) / (Vth_front + (x(3) * Wf / 0.2e1 + x(1) + delta_fr * x(3) * Lf)^2 / Vth_front);
else
  x8_ss = -(x(3) * Lf + x(2) + (-x(3) * Wf / 0.2e1 - x(1)) * delta_fr) / (x(3) * Wf / 0.2e1 + x(1) + delta_fr * x(3) * Lf);
end
dotx(8) = -abs(x(1))*(x(8) - x8_ss)/ly_front;

% State x(9) --> Side slip angle, front left wheel: alpha_fl(t)
if (x(3) * Wf / 0.2e1 + x(1) + delta_fl * x(3) * Lf) <= Vth_front
  x9_ss = 2*(-(x(3) * Lf + x(2) + (x(3) * Wf / 0.2e1 - x(1)) * delta_fl)) / (Vth_front + (x(3) * Wf / 0.2e1 + x(1) + delta_fl * x(3) * Lf)^2 / Vth_front);
else
  x9_ss = -(x(3) * Lf + x(2) + (x(3) * Wf / 0.2e1 - x(1)) * delta_fl) / (x(3) * Wf / 0.2e1 + x(1) + delta_fl * x(3) * Lf);
end
dotx(9) = -abs(x(1))*(x(9) - x9_ss)/ly_front;

% State x(10) --> Longitudinal slip, rear right wheel: kappa_rr(t)
if max(abs((x(19) * Rr)),  abs((x(3) * Wr) / 0.2e1 + x(1))) <= Vth_rear
  x10_ss = 2*((x(19) * Rr) - ((x(3) * Wr) / 0.2e1 + x(1))) / (Vth_rear + (max(abs((x(19) * Rr)),  abs((x(3) * Wr) / 0.2e1 + x(1))))^2 / Vth_rear);
else
  x10_ss = ((x(19) * Rr) - ((x(3) * Wr) / 0.2e1 + x(1))) / max(abs((x(19) * Rr)),  abs((x(3) * Wr) / 0.2e1 + x(1)));
end
dotx(10) = -abs(x(1))*(x(10) - x10_ss)/lx_rear;

% if (x(10)>=1 && dotx(10)>=0) || (x(10)<=-1 && dotx(10)<=0)
%   dotx(10) = 0;
% end

% State x(11) --> Longitudinal slip, rear left wheel: kappa_rl(t)
if max(abs((x(20) * Rr)), abs(-(x(3) * Wr) / 0.2e1 + x(1))) <= Vth_rear
  x11_ss = 2*((x(20) * Rr) - (-(x(3) * Wr) / 0.2e1 + x(1))) / (Vth_rear + (max(abs((x(20) * Rr)), abs(-(x(3) * Wr) / 0.2e1 + x(1))))^2 / Vth_rear);
else
  x11_ss = ((x(20) * Rr) - (-(x(3) * Wr) / 0.2e1 + x(1))) / max(abs((x(20) * Rr)), abs(-(x(3) * Wr) / 0.2e1 + x(1)));
end
dotx(11) = -abs(x(1))*(x(11) - x11_ss)/lx_rear;

% if (x(11)>=1 && dotx(11)>=0) || (x(11)<=-1 && dotx(11)<=0)
%   dotx(11) = 0;
% end


% State x(12) --> Longitudinal slip, front right wheel: kappa_fr(t)
if max(abs((x(21) * Rf)),  abs((x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf))) <= Vth_front
  x12_ss = 2*((x(21) * Rf) - ((x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf))) / (Vth_front + (max(abs((x(21) * Rf)),  abs((x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf))))^2 / Vth_front);
else
  x12_ss = ((x(21) * Rf) - ((x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf))) / max(abs((x(21) * Rf)),  abs((x(3) * Wf) / 0.2e1 + x(1) + (delta_fr * x(3) * Lf)));
end
dotx(12) = -abs(x(1))*(x(12) - x12_ss)/lx_front;

if (x(12)>=1 && dotx(12)>=0) || (x(12)<=-1 && dotx(12)<=0)
  dotx(12) = 0;
end

% State x(13) --> Longitudinal slip, front left wheel: kappa_fl(t)
if max(abs((x(22) * Rf)), abs(-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf))) <= Vth_front
  x13_ss = 2*((x(22) * Rf) - (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf))) / (Vth_front + (max(abs((x(22) * Rf)), abs(-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf))))^2/Vth_front );
else
  x13_ss = ((x(22) * Rf) - (-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf))) / max(abs((x(22) * Rf)), abs(-(x(3) * Wf) / 0.2e1 + x(1) + (delta_fl * x(3) * Lf)));
end
dotx(13) = -abs(x(1))*(x(13) - x13_ss)/lx_front;

if (x(13)>=1 && dotx(13)>=0) || (x(13)<=-1 && dotx(13)<=0)
  dotx(13) = 0;
end

% State x(14) --> Vertical load, rear right wheel: Fz_rr(t)
dotx(14) = ((4 * IXXs0 * L * x(14) * Wr * i_zz * m - 4 * IXXs0 * L * x(4) ^ 2 * Ms_r * hs ^ 2 * ms ^ 2 + 4 * IXXs0 * AM * x(4) ^ 2 * hr * hs ^ 2 * ms ^ 2 + 4 * Lf * x(3) * x(1) * hr * is_xz ^ 2 * m ^ 2 + 2 * x(3) * x(2) * Wr * h_G * is_xz ^ 2 * m ^ 2 - 2 * x(4) ^ 2 * hs * ms * h_G * Wr * IYZs0 * x(3) ^ 2 * is_xz * m + 2 * x(4) * hs * ms * h_G * Wr * Lf * Fy_fr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * Lf * Fy_fl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * Lr * Fy_rl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * Lr * Fy_rr * IXXs0 * m + x(4) * hs * ms * h_G * Wr * Wf * Fx_fr * IXXs0 * m - x(4) * hs * ms * h_G * Wr * Wf * Fx_fl * IXXs0 * m + 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * Rr * hs ^ 2 * m_ur * ms ^ 2 - 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m_ur * ms ^ 2 + 4 * IXXs0 * Lf * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m * ms ^ 2 + 4 * IXXs0 * x(3) * x(18) * Wr * h_G * hs * i_zz * m * ms + 2 * Lf * Wr * g * is_xz ^ 2 * m ^ 2 + 2 * L * FAz_r * Wr * is_xz ^ 2 * m + 2 * Fx_rl * Wr * h_G * is_xz ^ 2 * m + 2 * Fx_fr * Wr * h_G * is_xz ^ 2 * m - 2 * FAxc * Wr * h_G * is_xz ^ 2 * m + 2 * Fx_rr * Wr * h_G * is_xz ^ 2 * m + 2 * Fx_fl * Wr * h_G * is_xz ^ 2 * m + 4 * IXXs0 * L * Ms_r * i_zz * m - 4 * IXXs0 * AM * hr * i_zz * m - 2 * delta_fl * Fy_fl * Wr * h_G * is_xz ^ 2 * m - 2 * delta_fr * Fy_fr * Wr * h_G * is_xz ^ 2 * m - 2 * IXXs0 * L * FAz_r * Wr * i_zz * m - 2 * IXXs0 * Fx_rl * Wr * h_G * i_zz * m - 2 * IXXs0 * Fx_fr * Wr * h_G * i_zz * m + 2 * IXXs0 * FAxc * Wr * h_G * i_zz * m - 2 * IXXs0 * Fx_rr * Wr * h_G * i_zz * m - 2 * IXXs0 * Fx_fl * Wr * h_G * i_zz * m - 2 * IXXs0 * Lf * Wr * g * i_zz * m ^ 2 + 2 * x(4) * hs * ms * h_G * Wr * AM_rr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * AM_rl * IXXs0 * m + 2 * IXXs0 * Lf * x(4) ^ 2 * Wr * g * hs ^ 2 * m * ms ^ 2 + 2 * x(4) * hs * ms * h_G * Wr * Ms_r * is_xz * m + 2 * x(4) * hs * ms * h_G * Wr * Ms_f * is_xz * m + 2 * x(4) * hs * ms * h_G * Wr * AM_fl * IXXs0 * m - 4 * x(3) * x(18) * Wr * h_G * hs * is_xz ^ 2 * m * ms + 2 * x(4) ^ 2 * hs ^ 2 * ms ^ 2 * h_G * Wr * g * is_xz * m - x(4) * hs * ms * h_G * Wr ^ 2 * Fx_rl * IXXs0 * m + 2 * x(4) * hs ^ 2 * ms ^ 2 * h_G * Wr * ay * is_xz * m - 4 * IXXs0 * L * x(3) * x(1) * Rr * i_zz * m * m_ur + 4 * IXXs0 * L * x(3) * x(1) * hr * i_zz * m * m_ur + x(4) * hs * ms * h_G * Wr ^ 2 * Fx_rr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * AM_fr * IXXs0 * m - 4 * L * Ms_r * is_xz ^ 2 * m + 4 * AM * hr * is_xz ^ 2 * m - 4 * L * x(3) * x(1) * hr * is_xz ^ 2 * m * m_ur - 4 * IXXs0 * L * x(14) * x(4) ^ 2 * Wr * hs ^ 2 * ms ^ 2 + 2 * IXXs0 * L * FAz_r * x(4) ^ 2 * Wr * hs ^ 2 * ms ^ 2 - 4 * IXXs0 * Lf * x(3) * x(1) * hr * i_zz * m ^ 2 - 2 * IXXs0 * x(3) * x(2) * Wr * h_G * i_zz * m ^ 2 + 4 * L * x(3) * x(1) * Rr * is_xz ^ 2 * m * m_ur - 4 * L * x(14) * Wr * is_xz ^ 2 * m + 2 * IXXs0 * delta_fl * Fy_fl * Wr * h_G * i_zz * m + 2 * IXXs0 * delta_fr * Fy_fr * Wr * h_G * i_zz * m + 2 * x(4) * hs * ms * h_G * Wr * delta_fl * Lf * Fx_fl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * delta_fr * Lf * Fx_fr * IXXs0 * m + x(4) * hs * ms * h_G * Wr * delta_fl * Wf * Fy_fl * IXXs0 * m - x(4) * hs * ms * h_G * Wr * delta_fr * Wf * Fy_fr * IXXs0 * m) / L / tau_N_rear / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m) / Wr) / 0.4e1;

% State x(15) --> Vertical load, rear left wheel: Fz_rl(t)
dotx(15) = -((-4 * IXXs0 * L * x(15) * Wr * i_zz * m - 4 * IXXs0 * L * x(4) ^ 2 * Ms_r * hs ^ 2 * ms ^ 2 + 4 * IXXs0 * AM * x(4) ^ 2 * hr * hs ^ 2 * ms ^ 2 + 4 * Lf * x(3) * x(1) * hr * is_xz ^ 2 * m ^ 2 - 2 * x(3) * x(2) * Wr * h_G * is_xz ^ 2 * m ^ 2 + 2 * x(4) ^ 2 * hs * ms * h_G * Wr * IYZs0 * x(3) ^ 2 * is_xz * m - 2 * x(4) * hs * ms * h_G * Wr * Lf * Fy_fr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * Lf * Fy_fl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * Lr * Fy_rl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wr * Lr * Fy_rr * IXXs0 * m - x(4) * hs * ms * h_G * Wr * Wf * Fx_fr * IXXs0 * m + x(4) * hs * ms * h_G * Wr * Wf * Fx_fl * IXXs0 * m + 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * Rr * hs ^ 2 * m_ur * ms ^ 2 - 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m_ur * ms ^ 2 + 4 * IXXs0 * Lf * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m * ms ^ 2 - 4 * IXXs0 * x(3) * x(18) * Wr * h_G * hs * i_zz * m * ms - 2 * Lf * Wr * g * is_xz ^ 2 * m ^ 2 - 2 * L * FAz_r * Wr * is_xz ^ 2 * m - 2 * Fx_rl * Wr * h_G * is_xz ^ 2 * m - 2 * Fx_fr * Wr * h_G * is_xz ^ 2 * m + 2 * FAxc * Wr * h_G * is_xz ^ 2 * m - 2 * Fx_rr * Wr * h_G * is_xz ^ 2 * m - 2 * Fx_fl * Wr * h_G * is_xz ^ 2 * m + 4 * IXXs0 * L * Ms_r * i_zz * m - 4 * IXXs0 * AM * hr * i_zz * m + 2 * delta_fl * Fy_fl * Wr * h_G * is_xz ^ 2 * m + 2 * delta_fr * Fy_fr * Wr * h_G * is_xz ^ 2 * m + 2 * IXXs0 * L * FAz_r * Wr * i_zz * m + 2 * IXXs0 * Fx_rl * Wr * h_G * i_zz * m + 2 * IXXs0 * Fx_fr * Wr * h_G * i_zz * m - 2 * IXXs0 * FAxc * Wr * h_G * i_zz * m + 2 * IXXs0 * Fx_rr * Wr * h_G * i_zz * m + 2 * IXXs0 * Fx_fl * Wr * h_G * i_zz * m + 2 * IXXs0 * Lf * Wr * g * i_zz * m ^ 2 - 2 * x(4) * hs * ms * h_G * Wr * AM_rr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * AM_rl * IXXs0 * m - 2 * IXXs0 * Lf * x(4) ^ 2 * Wr * g * hs ^ 2 * m * ms ^ 2 - 2 * x(4) * hs * ms * h_G * Wr * Ms_r * is_xz * m - 2 * x(4) * hs * ms * h_G * Wr * Ms_f * is_xz * m - 2 * x(4) * hs * ms * h_G * Wr * AM_fl * IXXs0 * m + 4 * x(3) * x(18) * Wr * h_G * hs * is_xz ^ 2 * m * ms - 2 * x(4) ^ 2 * hs ^ 2 * ms ^ 2 * h_G * Wr * g * is_xz * m + x(4) * hs * ms * h_G * Wr ^ 2 * Fx_rl * IXXs0 * m - 2 * x(4) * hs ^ 2 * ms ^ 2 * h_G * Wr * ay * is_xz * m - 4 * IXXs0 * L * x(3) * x(1) * Rr * i_zz * m * m_ur + 4 * IXXs0 * L * x(3) * x(1) * hr * i_zz * m * m_ur - x(4) * hs * ms * h_G * Wr ^ 2 * Fx_rr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * AM_fr * IXXs0 * m - 4 * L * Ms_r * is_xz ^ 2 * m + 4 * AM * hr * is_xz ^ 2 * m + 4 * IXXs0 * L * x(15) * x(4) ^ 2 * Wr * hs ^ 2 * ms ^ 2 - 4 * L * x(3) * x(1) * hr * is_xz ^ 2 * m * m_ur - 2 * IXXs0 * L * FAz_r * x(4) ^ 2 * Wr * hs ^ 2 * ms ^ 2 - 4 * IXXs0 * Lf * x(3) * x(1) * hr * i_zz * m ^ 2 + 2 * IXXs0 * x(3) * x(2) * Wr * h_G * i_zz * m ^ 2 + 4 * L * x(3) * x(1) * Rr * is_xz ^ 2 * m * m_ur - 2 * IXXs0 * delta_fl * Fy_fl * Wr * h_G * i_zz * m - 2 * IXXs0 * delta_fr * Fy_fr * Wr * h_G * i_zz * m - 2 * x(4) * hs * ms * h_G * Wr * delta_fl * Lf * Fx_fl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wr * delta_fr * Lf * Fx_fr * IXXs0 * m - x(4) * hs * ms * h_G * Wr * delta_fl * Wf * Fy_fl * IXXs0 * m + x(4) * hs * ms * h_G * Wr * delta_fr * Wf * Fy_fr * IXXs0 * m + 4 * L * x(15) * Wr * is_xz ^ 2 * m) / L / tau_N_rear / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m) / Wr) / 0.4e1;

% State x(16) --> Vertical load, front right wheel: Fz_fr(t)
dotx(16) = ((-4 * IXXs0 * AM * x(4) ^ 2 * hr * hs ^ 2 * ms ^ 2 - x(4) * hs * ms * h_G * Wf * Wr * Fx_rr * IXXs0 * m + 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * Rf * hs ^ 2 * m_uf * ms ^ 2 - 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m_uf * ms ^ 2 + 4 * IXXs0 * Lr * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m * ms ^ 2 - 4 * IXXs0 * x(3) * x(18) * Wf * h_G * hs * i_zz * m * ms + 2 * x(4) ^ 2 * hs * ms * h_G * Wf * IYZs0 * x(3) ^ 2 * is_xz * m - x(4) * hs * ms * h_G * Wf ^ 2 * delta_fl * Fy_fl * IXXs0 * m + x(4) * hs * ms * h_G * Wf ^ 2 * delta_fr * Fy_fr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * Lf * Fy_fr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * Lf * Fy_fl * IXXs0 * m + 2 * Lr * Wf * g * is_xz ^ 2 * m ^ 2 + 4 * IXXs0 * AM * hr * i_zz * m + 2 * delta_fl * Fy_fl * Wf * h_G * is_xz ^ 2 * m + 2 * delta_fr * Fy_fr * Wf * h_G * is_xz ^ 2 * m - 2 * IXXs0 * L * FAz_f * Wf * i_zz * m + 2 * IXXs0 * Fx_rl * Wf * h_G * i_zz * m + 2 * IXXs0 * Fx_fr * Wf * h_G * i_zz * m - 2 * IXXs0 * FAxc * Wf * h_G * i_zz * m + 2 * IXXs0 * Fx_rr * Wf * h_G * i_zz * m + 2 * IXXs0 * Fx_fl * Wf * h_G * i_zz * m - 2 * IXXs0 * Lr * Wf * g * i_zz * m ^ 2 - x(4) * hs * ms * h_G * Wf ^ 2 * Fx_fr * IXXs0 * m - 2 * x(4) * hs ^ 2 * ms ^ 2 * h_G * Wf * ay * is_xz * m - 2 * x(4) * hs * ms * h_G * Wf * Ms_f * is_xz * m - 2 * x(4) ^ 2 * hs ^ 2 * ms ^ 2 * h_G * Wf * g * is_xz * m + 2 * IXXs0 * Lr * x(4) ^ 2 * Wf * g * hs ^ 2 * m * ms ^ 2 + 4 * IXXs0 * L * x(3) * x(1) * hr * i_zz * m * m_uf - 2 * x(4) * hs * ms * h_G * Wf * AM_rl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * AM_fr * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * AM_fl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * AM_rr * IXXs0 * m + x(4) * hs * ms * h_G * Wf ^ 2 * Fx_fl * IXXs0 * m - 4 * IXXs0 * L * x(3) * x(1) * Rf * i_zz * m * m_uf + 4 * x(3) * x(18) * Wf * h_G * hs * is_xz ^ 2 * m * ms - 2 * x(4) * hs * ms * h_G * Wf * Ms_r * is_xz * m + 4 * IXXs0 * L * x(16) * Wf * i_zz * m - 4 * IXXs0 * L * x(4) ^ 2 * Ms_f * hs ^ 2 * ms ^ 2 + 4 * Lr * x(3) * x(1) * hr * is_xz ^ 2 * m ^ 2 - 2 * x(3) * x(2) * Wf * h_G * is_xz ^ 2 * m ^ 2 - 4 * L * Ms_f * is_xz ^ 2 * m - 4 * AM * hr * is_xz ^ 2 * m - 4 * IXXs0 * L * x(16) * x(4) ^ 2 * Wf * hs ^ 2 * ms ^ 2 + 2 * IXXs0 * L * FAz_f * x(4) ^ 2 * Wf * hs ^ 2 * ms ^ 2 - 4 * IXXs0 * Lr * x(3) * x(1) * hr * i_zz * m ^ 2 + 2 * IXXs0 * x(3) * x(2) * Wf * h_G * i_zz * m ^ 2 + 4 * L * x(3) * x(1) * Rf * is_xz ^ 2 * m * m_uf - 4 * L * x(3) * x(1) * hr * is_xz ^ 2 * m * m_uf - 4 * L * x(16) * Wf * is_xz ^ 2 * m + 2 * x(4) * hs * ms * h_G * Wf * Lr * Fy_rl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * Lr * Fy_rr * IXXs0 * m + x(4) * hs * ms * h_G * Wf * Wr * Fx_rl * IXXs0 * m - 2 * IXXs0 * delta_fl * Fy_fl * Wf * h_G * i_zz * m - 2 * IXXs0 * delta_fr * Fy_fr * Wf * h_G * i_zz * m + 2 * L * FAz_f * Wf * is_xz ^ 2 * m - 2 * Fx_rl * Wf * h_G * is_xz ^ 2 * m - 2 * Fx_fr * Wf * h_G * is_xz ^ 2 * m + 2 * FAxc * Wf * h_G * is_xz ^ 2 * m - 2 * Fx_rr * Wf * h_G * is_xz ^ 2 * m - 2 * Fx_fl * Wf * h_G * is_xz ^ 2 * m + 4 * IXXs0 * L * Ms_f * i_zz * m - 2 * x(4) * hs * ms * h_G * Wf * delta_fl * Lf * Fx_fl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * delta_fr * Lf * Fx_fr * IXXs0 * m) / L / tau_N_front / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m) / Wf) / 0.4e1;

% State x(17) --> Vertical load, fromt left wheel: Fz_fl(t)
dotx(17) = -((-4 * IXXs0 * AM * x(4) ^ 2 * hr * hs ^ 2 * ms ^ 2 + x(4) * hs * ms * h_G * Wf * Wr * Fx_rr * IXXs0 * m + 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * Rf * hs ^ 2 * m_uf * ms ^ 2 - 4 * IXXs0 * L * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m_uf * ms ^ 2 + 4 * IXXs0 * Lr * x(3) * x(1) * x(4) ^ 2 * hr * hs ^ 2 * m * ms ^ 2 + 4 * IXXs0 * x(3) * x(18) * Wf * h_G * hs * i_zz * m * ms - 2 * x(4) ^ 2 * hs * ms * h_G * Wf * IYZs0 * x(3) ^ 2 * is_xz * m + x(4) * hs * ms * h_G * Wf ^ 2 * delta_fl * Fy_fl * IXXs0 * m - x(4) * hs * ms * h_G * Wf ^ 2 * delta_fr * Fy_fr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * Lf * Fy_fr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * Lf * Fy_fl * IXXs0 * m - 2 * Lr * Wf * g * is_xz ^ 2 * m ^ 2 + 4 * IXXs0 * AM * hr * i_zz * m - 2 * delta_fl * Fy_fl * Wf * h_G * is_xz ^ 2 * m - 2 * delta_fr * Fy_fr * Wf * h_G * is_xz ^ 2 * m + 2 * IXXs0 * L * FAz_f * Wf * i_zz * m - 2 * IXXs0 * Fx_rl * Wf * h_G * i_zz * m - 2 * IXXs0 * Fx_fr * Wf * h_G * i_zz * m + 2 * IXXs0 * FAxc * Wf * h_G * i_zz * m - 2 * IXXs0 * Fx_rr * Wf * h_G * i_zz * m - 2 * IXXs0 * Fx_fl * Wf * h_G * i_zz * m + 2 * IXXs0 * Lr * Wf * g * i_zz * m ^ 2 + x(4) * hs * ms * h_G * Wf ^ 2 * Fx_fr * IXXs0 * m + 2 * x(4) * hs ^ 2 * ms ^ 2 * h_G * Wf * ay * is_xz * m + 2 * x(4) * hs * ms * h_G * Wf * Ms_f * is_xz * m + 2 * x(4) ^ 2 * hs ^ 2 * ms ^ 2 * h_G * Wf * g * is_xz * m - 2 * IXXs0 * Lr * x(4) ^ 2 * Wf * g * hs ^ 2 * m * ms ^ 2 + 4 * IXXs0 * L * x(3) * x(1) * hr * i_zz * m * m_uf + 2 * x(4) * hs * ms * h_G * Wf * AM_rl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * AM_fr * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * AM_fl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * AM_rr * IXXs0 * m - x(4) * hs * ms * h_G * Wf ^ 2 * Fx_fl * IXXs0 * m - 4 * IXXs0 * L * x(3) * x(1) * Rf * i_zz * m * m_uf - 4 * x(3) * x(18) * Wf * h_G * hs * is_xz ^ 2 * m * ms + 2 * x(4) * hs * ms * h_G * Wf * Ms_r * is_xz * m - 4 * IXXs0 * L * x(17) * Wf * i_zz * m - 4 * IXXs0 * L * x(4) ^ 2 * Ms_f * hs ^ 2 * ms ^ 2 + 4 * Lr * x(3) * x(1) * hr * is_xz ^ 2 * m ^ 2 + 2 * x(3) * x(2) * Wf * h_G * is_xz ^ 2 * m ^ 2 + 4 * L * x(17) * Wf * is_xz ^ 2 * m - 4 * L * Ms_f * is_xz ^ 2 * m - 4 * AM * hr * is_xz ^ 2 * m + 4 * IXXs0 * L * x(17) * x(4) ^ 2 * Wf * hs ^ 2 * ms ^ 2 - 2 * IXXs0 * L * FAz_f * x(4) ^ 2 * Wf * hs ^ 2 * ms ^ 2 - 4 * IXXs0 * Lr * x(3) * x(1) * hr * i_zz * m ^ 2 - 2 * IXXs0 * x(3) * x(2) * Wf * h_G * i_zz * m ^ 2 + 4 * L * x(3) * x(1) * Rf * is_xz ^ 2 * m * m_uf - 4 * L * x(3) * x(1) * hr * is_xz ^ 2 * m * m_uf - 2 * x(4) * hs * ms * h_G * Wf * Lr * Fy_rl * IXXs0 * m - 2 * x(4) * hs * ms * h_G * Wf * Lr * Fy_rr * IXXs0 * m - x(4) * hs * ms * h_G * Wf * Wr * Fx_rl * IXXs0 * m + 2 * IXXs0 * delta_fl * Fy_fl * Wf * h_G * i_zz * m + 2 * IXXs0 * delta_fr * Fy_fr * Wf * h_G * i_zz * m - 2 * L * FAz_f * Wf * is_xz ^ 2 * m + 2 * Fx_rl * Wf * h_G * is_xz ^ 2 * m + 2 * Fx_fr * Wf * h_G * is_xz ^ 2 * m - 2 * FAxc * Wf * h_G * is_xz ^ 2 * m + 2 * Fx_rr * Wf * h_G * is_xz ^ 2 * m + 2 * Fx_fl * Wf * h_G * is_xz ^ 2 * m + 4 * IXXs0 * L * Ms_f * i_zz * m + 2 * x(4) * hs * ms * h_G * Wf * delta_fl * Lf * Fx_fl * IXXs0 * m + 2 * x(4) * hs * ms * h_G * Wf * delta_fr * Lf * Fx_fr * IXXs0 * m) / L / tau_N_front / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m) / Wf) / 0.4e1;

% State x(18) --> Roll rate: phi_dot(t)
dotx(18) = -((-2 * IYZs0 * x(3) ^ 2 * x(4) * i_zz * m - 2 * x(4) ^ 3 * g * hs ^ 3 * ms ^ 3 - 2 * x(4) ^ 2 * Ms_f * hs ^ 2 * ms ^ 2 - 2 * x(4) ^ 2 * Ms_r * hs ^ 2 * ms ^ 2 - 2 * ay * x(4) ^ 2 * hs ^ 3 * ms ^ 3 + 2 * is_xz * delta_fl * Lf * Fx_fl * m + 2 * is_xz * delta_fr * Lf * Fx_fr * m + is_xz * delta_fl * Wf * Fy_fl * m - is_xz * delta_fr * Wf * Fy_fr * m + 2 * ay * hs * i_zz * m * ms + 2 * is_xz * AM_rl * m + 2 * is_xz * AM_rr * m + 2 * is_xz * AM_fl * m + 2 * is_xz * AM_fr * m + 2 * Ms_f * i_zz * m + 2 * Ms_r * i_zz * m + is_xz * Wf * Fx_fr * m + is_xz * Wr * Fx_rr * m - is_xz * Wf * Fx_fl * m + 2 * is_xz * Lf * Fy_fr * m + 2 * is_xz * Lf * Fy_fl * m - 2 * is_xz * Lr * Fy_rl * m - 2 * is_xz * Lr * Fy_rr * m - is_xz * Wr * Fx_rl * m + 4 * is_xz * x(3) * x(18) * x(4) * hs ^ 2 * ms ^ 2 + 2 * is_xz * delta_fl * x(4) * Fy_fl * hs * ms + 2 * is_xz * delta_fr * Fy_fr * x(4) * hs * ms + 2 * x(4) * g * hs * i_zz * m * ms - 2 * is_xz * Fx_rl * x(4) * hs * ms - 2 * is_xz * Fx_fr * x(4) * hs * ms + 2 * is_xz * FAxc * x(4) * hs * ms - 2 * is_xz * Fx_rr * x(4) * hs * ms - 2 * is_xz * x(4) * Fx_fl * hs * ms + 2 * IYZs0 * x(3) ^ 2 * x(4) ^ 3 * hs ^ 2 * ms ^ 2 - 2 * is_xz * x(3) * x(2) * x(4) * hs * m * ms) / (hs ^ 2 * ms ^ 2 * x(4) ^ 2 * IXXs0 - i_zz * IXXs0 * m + is_xz ^ 2 * m)) / 0.2e1;

% State x(19) -> angular velocity rear right wheel: omega_rr(t)
dotx(19) = -0.1e1 / iwa_r * (Fx_rr * Rr - x(24)*regSign(x(19)));
% dotx(19) = -0.1e1 / iwa_r * (Fx_rr * Rr - x(24));


% State x(20) -> angular velocity rear left wheel: omega_rl(t)
dotx(20) = -0.1e1 / iwa_r * (Fx_rl * Rr - x(25)*regSign(x(20)));
% dotx(20) = -0.1e1 / iwa_r * (Fx_rl * Rr - x(25));


% State x(21) -> angular velocity front right wheel: omega_fr(t)
dotx(21) = -0.1e1 / iwa_f * (Fx_fr * Rf - x(26)*regSign(x(21)));
% dotx(21) = -0.1e1 / iwa_f * (Fx_fr * Rf - x(26));


% State x(22) -> angular velocity front left wheel: omega_fl(t)
dotx(22) = -0.1e1 / iwa_f * (Fx_fl * Rf - x(27)*regSign(x(22)));
% dotx(22) = -0.1e1 / iwa_f * (Fx_fl * Rf - x(27));

tau = 1e-2;
% State x(23) -> Steering angle from the actuator (low-pass filtered)
dotx(23) = (u(1) - x(23))/1e-3; % delta
dotx(24) = (u(2) - x(24))/tau; % Tw_rr
dotx(25) = (u(3) - x(25))/tau; % Tw_rl
dotx(26) = (u(4) - x(26))/tau; % Tw_fr
dotx(27) = (u(5) - x(27))/tau; % Tw_fl

% if(x(19)) <= 1
%   dotx(24)=0;
% end
% if(x(20)) <= 1
%   dotx(25)=0;
% end
% if(x(21)) <= 1
%   dotx(26)=0;
% end
% if(x(22)) <= 1
%   dotx(27)=0;
% end

% fprintf('Tw_rr = %3.3f  Fx_rr*Rr = %3.3f  ',u(2),Fx_rr * Rr)
% fprintf('Tw_rr_real = %3.3f\n----------------\n',x(24))

block.Derivatives.Data = dotx;
%end Derivatives

function SetInpPortFrameData(block, idx, fd)

block.InputPort(idx).SamplingMode = fd;
block.OutputPort(1).SamplingMode  = fd;

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

