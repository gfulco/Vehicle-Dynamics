
%% Parameters definition for PID Simulink

clear
close all
clc

% Define LaTeX as interpreter for titlr, labels and legend in plots
set(0,'defaulttextinterpreter','latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath('Figures')

%% Parameter definition

lambda = linspace(0,1,100); % [-] lambda vector

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

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Parameters State Machine (as defined in class)
rho_off = 0.9;
rho_on = 1.02;
rho = 0.97;
lambda1_th = 0.2;
lambda2_th = 0.1;
Tb_dot_max = 5e3;
v_on = 2.5;
hv = 0.5;
v_stop = 0.1;

% Simulation parameters
t0 = 0;             % [s] Initial simulation time
tf = 5;             % [s] Final simulation time
dt = 1e-3;          % [s] Time interval
time = (t0:dt:tf)'; % [s] Time vector
brake_time = 1;     % [s] Brake time instant 

%% Simulate and resample data
sim('ABS_contDyn_model')

simTime = ABS.v.Time;

% Vehicle with ABS
v_ABS      = interp1(simTime,ABS.v.Data,time);
x_ABS      = interp1(simTime,ABS.x.Data,time);
lambda_ABS = interp1(simTime,ABS.lambda.Data,time);
Tb_ABS     = interp1(simTime,ABS.Tb.Data,time);

% Vehicle without ABS
v_noABS      = interp1(simTime,noABS.v.Data,time);
x_noABS      = interp1(simTime,noABS.x.Data,time);
lambda_noABS = interp1(simTime,noABS.lambda.Data,time);
Tb_noABS     = interp1(simTime,noABS.Tb.Data,time);

fprintf('Press enter to continue.\n')
pause
%% Plot Lambda values
psi_lambda = @(lambda,Fz) (auxdata.r_w + auxdata.J/(auxdata.r_w*auxdata.m).*(1-lambda)).*Fz.*burckhardt(lambda,auxdata.road_condition);

figure('Name','Long. Wheel Slip during braking','NumberTitle','off'),
hold on
plot(lambda,psi_lambda(lambda,1.5*auxdata.m*auxdata.g),'--k','LineWidth',2)
grid on
xlabel('$\lambda [-]$')
ylabel('$T_b [Nm]$')
yL = get(gca,'YLim');

hABS = plot([lambda_ABS(1),lambda_ABS(1)],yL,'color','red');
hnoABS = plot([lambda_noABS(1),lambda_ABS(1)],yL,'color','blue');
% legend([hABS,hnoABS],{'ABS','noABS'})
for i = 1:4:3/dt
    hABS = plot([lambda_ABS(i),lambda_ABS(i)],yL,'color','red','LineWidth',2);
    hnoABS = plot([lambda_noABS(i),lambda_noABS(i)],yL,'color','blue','LineWidth',2);
    pause(dt)
    set(hABS,'Visible','off')
    set(hnoABS,'Visible','off')
end
hABS = plot([lambda_ABS(i),lambda_ABS(i)],yL,'color','red','LineWidth',2);
hnoABS = plot([lambda_noABS(i),lambda_noABS(i)],yL,'color','blue','LineWidth',2);

fprintf('Press enter to continue.\n')
pause
%%
figure('Name','Long. Wheel Slip and Torque during braking','NumberTitle','off'),
hold on
plot(lambda,psi_lambda(lambda,1.5*auxdata.m*auxdata.g),'--k','LineWidth',2)
grid on
ylim([0 1500])
xlabel('$\lambda [-]$')
ylabel('$T_b [Nm]$')

for i = 1:4:3/dt
    hABS = plot(lambda_ABS(1:i),Tb_ABS(1:i),'color','red','Linewidth',2);
    hnoABS = plot(lambda_noABS(1:i),Tb_noABS(1:i),'color','blue','Linewidth',2);
    pause(dt)
end
legend([hABS,hnoABS],{'ABS','noABS'})

fprintf('Press enter to continue.\n')
pause
%% Plot braking maneuver
figure('Name','Braking Maneuver Comparison','NumberTitle','off'),
hold on
xlim([0 max(x_ABS(end),x_noABS(end))+length])
axis equal
grid on
yL = get(gca,'YLim');
x_brake = x_ABS(time==brake_time);
plot([x_brake x_brake],[0 50],'--k')
text(x_brake+1,5,'$\leftarrow$ Braking starting point')
ylim([0 50])
xlabel('$x$ [m]')
ylabel('$y$ [m]')
title(strcat('Braking maneuver on',road_condition_names{auxdata.road_condition}))
h1 = plot([x_ABS(1) x_ABS(1)],[15 15],'--r','LineWidth',2);
h2 = plot([x_noABS(1) x_noABS(1)],[30 30],'--b','LineWidth',2);
% legend([h1,h2],{'with ABS','without ABS'})
for i = 1:10:numel(time)
    plot([x_ABS(1) x_ABS(i)],[15 15],'--r','LineWidth',2)
    plot([x_noABS(1) x_noABS(i)],[30 30],'--b','LineWidth',2)
    hABS = rectangle('Position',[x_ABS(i)-length/2,15-width/2,length,width],'Curvature',[0.2,0.9],'FaceColor','red');
    hnoABS = rectangle('Position',[x_noABS(i)-length/2,30-width/2,length,width],'Curvature',[0.2,0.9],'FaceColor','blue');
    pause(dt)
    set(hABS,'Visible','off')
    set(hnoABS,'Visible','off')
end
hABS = rectangle('Position',[x_ABS(i)-length/2,15-width/2,length,width],'Curvature',[0.2,0.9],'FaceColor','red');
hnoABS = rectangle('Position',[x_noABS(i)-length/2,30-width/2,length,width],'Curvature',[0.2,0.9],'FaceColor','blue');




