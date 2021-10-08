function vehicle_data = getVehicleDataStruct()
%  ___                           _            ___       _        
% / __|_  _ ____ __  ___ _ _  __(_)___ _ _   |   \ __ _| |_ __ _ 
% \__ \ || (_-< '_ \/ -_) ' \(_-< / _ \ ' \  | |) / _` |  _/ _` |
% |___/\_,_/__/ .__/\___|_||_/__/_\___/_||_| |___/\__,_|\__\__,_|
%             |_|                                                

% REAR
rear_suspension.hrr           = 0.1;             % [m] Heigth of the rear roll axle
rear_suspension.Ks_r  = 20000/1.45;      % [N/m] Rear suspension stiffness
rear_suspension.Cs_r  = 3000.0000/1.45;  % [N*s/m] Rear suspension damping
% suspension kinematics
% gamma_r = a0+a1*phi+a2*phi^2+a3*phi^3
rear_suspension.camber_right_wheel     = [-0.53447e-1, +0.839386,    -0.201305,    -0.87213e-1];
rear_suspension.camber_left_wheel      = [ 0.53447e-1, +0.839386,     0.201305,    -0.87213e-1];
rear_suspension.roll_steer_right_wheel = [-0.26463e-2, -0.428094e-2, -0.119244e-1, +0.65563e-1];
rear_suspension.roll_steer_left_wheel  = [ 0.26463e-2, -0.428094e-2,  0.119244e-1, +0.65563e-1];

% FRONT
front_suspension.Ks_f         = 20000/1.45;      % [N*s/m] Front suspension stiffness
front_suspension.Cs_f         = 3000.0000/1.45;  % [N*s/m] Front suspension dumping
front_suspension.hrf          = 0.1;             % [m] Heigth of the front wheels axle
front_suspension.camber_right_wheel     = [-0.53447e-1 ,+0.839386, -0.201305, -0.87213e-1];
front_suspension.camber_left_wheel      = [ 0.53447e-1 ,+0.839386,   0.201305, -0.87213e-1];
front_suspension.roll_steer_right_wheel = [-0.26463e-2, -0.428094e-2, -0.119244e-1, +0.65563e-1];
front_suspension.roll_steer_left_wheel  = [ 0.26463e-2, -0.428094e-2,  0.119244e-1, +0.65563e-1];


%   ___ _               _      ___       _        
%  / __| |_  __ _ _____(_)___ |   \ __ _| |_ __ _ 
% | (__| ' \/ _` (_-<_-< (_-< | |) / _` |  _/ _` |
%  \___|_||_\__,_/__/__/_/__/ |___/\__,_|\__\__,_|
%                                                 
% CHASSIS IS THE SPRUNG BODY

% CHASSIS
% is =  |  is_xx   0   -is_xz |
%       |    0   is_yy    0   |
%       | -is_xz   0    is_zz |
chassis.is_xx = 85;
chassis.is_yy = 850;
chassis.is_zz = 800;
chassis.is_xz = 60;

%  _   _                                   ___       _        
% | | | |_ _  ____ __ _ _ _  _ _ _  __ _  |   \ __ _| |_ __ _ 
% | |_| | ' \(_-< '_ \ '_| || | ' \/ _` | | |) / _` |  _/ _` |
%  \___/|_||_/__/ .__/_|  \_,_|_||_\__, | |___/\__,_|\__\__,_|
%               |_|                |___/                      
% UNSRPUNG BODY IS MADE OF THE 4 WHEELS AND THE SUSPENSION, TRANSMISSION AND BRAKE MASSES ATTACHED TO WHEELS

% WHEEL
% iwd = | iwd   0  0  |
%       |  0  iwa  0  |
%       |  0   0  iwd |

% REAR
m_wr = 5;      % [kg] Wheel mass
w_wr = 0.2;    % [m] Wheel width
rr   = 0.328;  % [m] Rear Wheel Radius
rear_wheel.Rr      = rr;           % [m] Rear Wheel Radius
rear_wheel.width  = w_wr;                     % [m] Wheel width
rear_wheel.mass    = m_wr;         % [kg] Wheel mass
rear_wheel.iwd_r   = m_wr/12 * (3*rr^2 + w_wr^2);
rear_wheel.iwa_r   = m_wr/2 * rr^2;

m_ur = 15;
rear_unsprung.mass = m_ur;              % [kg] Rear unsprung mass

% FRONT
m_wf = 5;      % [kg] Wheel mass
w_wf = 0.2;    % [m] Wheel width
rf   = 0.328;  % [m] Front Wheel Radius

front_wheel.Rf      = rf;                       % [m] Front Wheel Radius
front_wheel.width   = w_wf;                     % [m] Wheel width
front_wheel.mass    = m_wf;                     % [kg] Wheel mass
front_wheel.iwd_f   = m_wf/12 * (3*rf^2 + w_wf^2); 
front_wheel.iwa_f   = m_wf/2 * rf^2;

m_uf = 15;
front_unsprung.mass = m_uf;                     % [kg] Front unsprung mass


%    ___                   _ _  __   __   _    _    _       ___       _        
%   / _ \__ _____ _ _ __ _| | | \ \ / /__| |_ (_)__| |___  |   \ __ _| |_ __ _ 
%  | (_) \ V / -_) '_/ _` | | |  \ V / -_) ' \| / _| / -_) | |) / _` |  _/ _` |
%   \___/ \_/\___|_| \__,_|_|_|   \_/\___|_||_|_\__|_\___| |___/\__,_|\__\__,_|
%                                                                              

% VEHICLE
L  = 3.1200;
Lr = 1.2 ;
Lf = L-Lr;
hs = 0.2950-0.1 ;
hr = (front_suspension.hrf-rear_suspension.hrr)*Lr/L + rear_suspension.hrr;
  
vehicle.Lf    = Lf;        % [m] Distance between vehicle CoG and front wheels axle
vehicle.Lr    = Lr;        % [m] Distance between vehicle CoG and front wheels axle
vehicle.L     = Lf+Lr;     % [m] Vehicle length
vehicle.hs    = hs;        % [m] Position of chassis CoG from roll axis
vehicle.hr    = hr;        % [m] Heigth of vehicle roll axis
vehicle.Wf    = 1.4780;    % [m] Width of front wheels axle 
vehicle.Wr    = 1.4150;    % [m] Width of rear wheels axle 
% Inertia and mass
vehicle.m     = 357.5;      % [kg] Vehicle Mass
vehicle.g     = 9.8100;     % [m/s^2] Gravitational acceleration
vehicle.i_xx  = 85;
vehicle.i_yy  = 850;
vehicle.i_zz  = 800;
vehicle.I_xz  = 60;

%   __  __                           _    ___      ___
%  |  \/  |__ _ ______  __ _ _ _  __| |  / __|___ / __|
%  | |\/| / _` (_-<_-< / _` | ' \/ _` | | (__/ _ \ (_ |
%  |_|  |_\__,_/__/__/ \__,_|_||_\__,_|  \___\___/\___|

vehicle.ms   = vehicle.m - m_uf - m_ur; % [kg] Sprung Mass
vehicle.h_G  = (vehicle.hs + vehicle.hr) / vehicle.m * vehicle.ms + (m_uf * rf + m_ur * rr) / vehicle.m; % [m] Position of the center of mass

%   ___              _   _
%  |_ _|_ _  ___ _ _| |_(_)__ _
%   | || ' \/ -_) '_|  _| / _` |
%  |___|_||_\___|_|  \__|_\__,_|

% Moment of inertia about X axis of RF0
vehicle.IXXs0 = vehicle.hs^2*vehicle.ms+chassis.is_xx;             % [kg*m^2] chassis moment of inertia w.r.t. roll axis
vehicle.IYZs0 = chassis.is_zz -vehicle.hs^2*vehicle.ms-chassis.is_yy;      % [kg*m^2] chassis moment of inertia w.r.t. roll axis
vehicle.IXX0  = (vehicle.hr+vehicle.hs)*vehicle.hs*vehicle.ms+chassis.is_xx;       % [kg*m^2] chassis moment of inertia w.r.t. RF0
vehicle.IYZ0  = vehicle.ms*(vehicle.hr+vehicle.hs)*vehicle.hs+chassis.is_yy-chassis.is_zz; % [kg*m^2] chassis moment of inertia w.r.t. RF0


%     _                   _                      _
%    /_\  ___ _ _ ___  __| |_  _ _ _  __ _ _ __ (_)__ ___
%   / _ \/ -_) '_/ _ \/ _` | || | ' \/ _` | '  \| / _(_-<
%  /_/ \_\___|_| \___/\__,_|\_, |_||_\__,_|_|_|_|_\__/__/
%                           |__/

aerodynamics.CAx   = 1.5563680;   % [N*s^2/m^2]
aerodynamics.CAzf  = .6412236160; % [N*s^2/m^2]
aerodynamics.CAzr  = .9151443840; % [N*s^2/m^2]

%   _____                       _       _            ___       _        
%  |_   _| _ __ _ _ _  ____ __ (_)_____(_)___ _ _   |   \ __ _| |_ __ _ 
%    | || '_/ _` | ' \(_-< '  \| (_-<_-< / _ \ ' \  | |) / _` |  _/ _` |
%    |_||_| \__,_|_||_/__/_|_|_|_/__/__/_\___/_||_| |___/\__,_|\__\__,_|
%                                                                       

%   ___ _               _             ___             ___       _        
%  / __| |_ ___ ___ _ _(_)_ _  __ _  / __|_  _ ___   |   \ __ _| |_ __ _ 
%  \__ \  _/ -_) -_) '_| | ' \/ _` | \__ \ || (_-<_  | |) / _` |  _/ _` |
%  |___/\__\___\___|_| |_|_||_\__, | |___/\_, /__(_) |___/\__,_|\__\__,_|
%                             |___/       |__/                           
steering_system.tau_D = 12;             % [-] Steering angle ratio
steering_system.tau_H = 0.1;            % [s] Time constant steering wheel

%  _____               ___       _        
% |_   _|  _ _ _ ___  |   \ __ _| |_ __ _ 
%   | || || | '_/ -_) | |) / _` |  _/ _` |
%   |_| \_, |_| \___| |___/\__,_|\__\__,_|
%       |__/                              


tau_N = 0.1;            % [s] Time constant vertical load
Vth   = 5;                % [m/s] Threshold velocity for low speed model switch
ly    = 0.12;              % [m] wheel lateral relaxation length
lx    = 0.05;              % [m] wheel longitudinal relaxation length

% Load Pacejka parameter model 1996 from external file.
% Cam be different or equal for front/rear tyres
pacejkaParam_r = loadPacejkaParam();
pacejkaParam_f = loadPacejkaParam();

rear_tyre.long_relax_length = lx ;
rear_tyre.lat_relax_length  = ly ;
rear_tyre.tau_N             = tau_N ;
rear_tyre.Vth               = Vth ;
rear_tyre.pacejkaParam      = pacejkaParam_r;

front_tyre.long_relax_length = lx ;
front_tyre.lat_relax_length  = ly ;
front_tyre.tau_N             = tau_N ;
front_tyre.Vth               = Vth ;
front_tyre.pacejkaParam       = pacejkaParam_f;




% __   __   _    _    _       ___ _               _     ___       _        
% \ \ / /__| |_ (_)__| |___  / __| |_ _ _ _  _ __| |_  |   \ __ _| |_ __ _ 
%  \ V / -_) ' \| / _| / -_) \__ \  _| '_| || / _|  _| | |) / _` |  _/ _` |
%   \_/\___|_||_|_\__|_\___| |___/\__|_|  \_,_\__|\__| |___/\__,_|\__\__,_|
% 
% Store all sub-structures of sub system data in vehicle structure
%
vehicle_data.chassis          = chassis;
vehicle_data.aerodynamics     = aerodynamics;
vehicle_data.steering_system  = steering_system;

%
vehicle_data.rear_unsprung    = rear_unsprung;
vehicle_data.front_unsprung   = front_unsprung;
%
vehicle_data.rear_wheel       = rear_wheel;
vehicle_data.front_wheel      = front_wheel;
%
vehicle_data.rear_suspension  = rear_suspension;
vehicle_data.front_suspension = front_suspension;
%
vehicle_data.rear_tyre        = rear_tyre;
vehicle_data.front_tyre       = front_tyre;

vehicle_data.vehicle          = vehicle;

end


