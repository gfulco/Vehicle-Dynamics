%% Example closed loop trajectories

clear
close all
clc

% Define LaTeX as interpreter for titlr, labels and legend in plots
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Parameter definition
v_bar = 40;             % [m/s]
tau_delay = 10e-3;      % [ms] Delay of the actuation system (caliper)
omega_act = 70;         % [rad/s] Bandwidth of the actuation system (caliper)
omega_LPF = 100;        % [rad/s] Bandwidth of the LowPass filter
auxdata.r_w = 0.3;      % [m] Wheel radius
auxdata.m = 225;        % [kg] Single corner mass
length = 4.387;         % [m] Vehicle length
width  = 1.768;         % [m] Vehicle width
auxdata.g = 9.81;       % [m/s^2] Gravitational aceeleration
m_wheel = 5;            % [kg] Wheel mass 
auxdata.J = 1;          % [kg m^2] wheel inertia
init_speed = 100;       % [km/h] Initial vehicle speed

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Select road  condition for the simulation:
%   road_condition = 1; --> Dry Asphalt
%   road_condition = 2; --> Wet Asphalt
%   road_condition = 3; --> Cobblestone
%   road_condition = 4; --> Snow
road_condition_names = {' Dry Asphalt' , ' Wet Asphalt' , ' Cobblestone' , ' Snow'};
auxdata.road_condition = 1; % Road Condition

fprintf(strcat('Road condition: ',road_condition_names{auxdata.road_condition},'\n'));
lambda = linspace(0,1,100);

%% Friction coefficients at different road conditions
figure('Name','Friction Coefficient','NumberTitle','off'),
hold on
grid on
road_condition = {'Dry Asphalt' , 'Wet Asphalt' , 'Cobblestone' , 'Snow'};
for i = 1:numel(road_condition)
    plot(lambda,burckhardt(lambda,i))
end
xlabel('$\lambda$')
ylabel('$\mu(\lambda)$')
title('Friction coefficient for different road conditions')
legend(road_condition)


%% Trajectories with positive braking torque

k = 5e3; % [Nm/s]
psi_lambda = @(lambda,Fz) (auxdata.r_w + auxdata.J/(auxdata.r_w*auxdata.m).*(1-lambda)).*Fz.*burckhardt(lambda,1);

% ~~~~~~~ Uncomment to simulate Limit Cycle ~~~~~~~ %
lambda_min = 0.08;
lambda_max = 0.35;
Tb_min = 600;
Tb_max = 900;


% % ~~~~~~~ Uncomment to simulate violation of condition 1 for the existence of the limit cycle ~~~~~~~ %
% lambda_min = 0.08;
% lambda_max = 0.35;
% Tb_min = 600;
% Tb_max = 750;

% ~~~~~~~ Uncomment to simulate violation of condition 2 for the existence of the limit cycle ~~~~~~~ %
% lambda_min = 0.08;
% lambda_max = 0.35;
% Tb_min = 900;
% Tb_max = 1100;

% ~~~~~~~ Uncomment to simulate violation of condition 3 for the existence of the limit cycle ~~~~~~~ %
% lambda_min = 0.08;
% lambda_max = 0.35;
% Tb_min = 700;
% Tb_max = 1200;

%% Simulate and Resample signals

Tb_init = 0;
lambda_init = 0;

t0 = 0;    % [s]
tf = 3;    % [s]
dt = 1e-3; % [s]
time = (t0:dt:tf)';

sim('HAB_control')


Tb_traj = interp1(Tb_sim.Time,Tb_sim.Data,time);
lambda_traj = interp1(lambda_sim.Time,lambda_sim.Data,time);

%% Plot State Trajectory
figure('Name','State Trajectory','NumberTitle','off'),
hold on
plot(lambda,psi_lambda(lambda,auxdata.m*auxdata.g),'--k','LineWidth',2)
title('Sate Trajectory with Switching Surfaces')
xL = get(gca,'XLim');
plot(xL,[Tb_min Tb_min],'--k')
plot(xL,[Tb_max Tb_max],'--k')
yL = get(gca,'YLim');
plot([lambda_min lambda_min],yL,'--k')
plot([lambda_max lambda_max],yL,'--k')
xlabel('$\lambda$ [-]')
ylabel('$T_b$ [Nm]')
grid on

fprintf('Press Enter to continue.\n')
pause

for i = 1:numel(Tb_traj)
    plot(lambda_traj(1:i),Tb_traj(1:i),'linewidth',1,'Color','magenta')
    circ = plot(lambda_traj(i),Tb_traj(i),'X','LineWidth',2);
    pause(dt)
    set(circ,'Visible','off')
end




