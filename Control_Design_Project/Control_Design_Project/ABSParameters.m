function parametri_ABS = ABSParameters(vehicle_data)
%ABSPARAMETERS Summary of this function goes here
%   Detailed explanation goes here

g     = 9.81;
Lf    = vehicle_data.vehicle.Lf;
Lr    = vehicle_data.vehicle.Lr;
L     = Lf+Lr;
Rf    = vehicle_data.front_wheel.Rf;        % Front Wheel Radius
Rr    = vehicle_data.rear_wheel.Rr;         % Rear Wheel Radius
m     = vehicle_data.vehicle.m;             % Vehicle Mass
Wf    = vehicle_data.vehicle.Wf;            % [m] Width of front wheels axle 
Wr    = vehicle_data.vehicle.Wr;            % [m] Width of rear wheels axle 

% Aerodynamics
CAx   = vehicle_data.aerodynamics.CAx;
CAzf  = vehicle_data.aerodynamics.CAzf;
CAzr  = vehicle_data.aerodynamics.CAzr;

%%  Continuous 
para.v_bar = 40;             % [m/s]
para.tau_delay = 10e-3;      % [ms] Delay of the actuation system (caliper)
para.omega_act = 70;         % [rad/s] Bandwidth of the actuation system (caliper)
para.omega_LPF = 100;        % [rad/s] Bandwidth of the LowPass filter
para.r_w = 0.3;              % [m] Wheel radius
para.m = 225;                % [kg] Single corner mass
para.length = 4.387;         % [m] Vehicle length
para.width  = 1.768;         % [m] Vehicle width
para.g = 9.81;               % [m/s^2] Gravitational aceeleration
para.m_wheel = 5;            % [kg] Wheel mass 
para.J = 1;                  % [kg m^2] wheel inertia
para.init_speed = 100;       % [km/h] Initial vehicle speed

%% Discrete


% Parameters State Machine (as defined in class)
para.rho_off = 0.9;
para.rho_on = 1.02;
para.rho = 0.97;
para.lambda1_th = 0.18;
para.lambda2_th = 0.07;
para.Tb_dot_max = 5e3;
para.v_on = 2.5;
para.hv = 0.5;
para.v_stop = 0.1;

para.Fzr0 = (m * g * Lf) / L 
para.Fzf0 = (Lr * g * m) / L

parametri_ABS=para;

end

