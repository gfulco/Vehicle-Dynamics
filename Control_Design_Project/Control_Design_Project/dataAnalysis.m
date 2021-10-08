%% PostProcessing and data Analysis


%   ___       _
%  |   \ __ _| |_ __ _
%  | |) / _` |  _/ _` |
%  |___/\__,_|\__\__,_|
%


Lf =  vehicle_data.vehicle.Lf ;  % [m] Distance between vehicle CoG and front wheels axle
Lr =  vehicle_data.vehicle.Lr ;  % [m] Distance between vehicle CoG and front wheels axle
Lr =  vehicle_data.vehicle.L  ;  % [m] Vehicle length
hs =  vehicle_data.vehicle.hs ;  % [m] Position of chassis CoG from roll axis
hr =  vehicle_data.vehicle.hr ;  % [m] Heigth of vehicle roll axis
Wf =  vehicle_data.vehicle.Wf ;  % [m] Width of front wheels axle 
Wr =  vehicle_data.vehicle.Wr ;  % [m] Width of rear wheels axle 
                 
m =   vehicle_data.vehicle.m  ;   % [kg] Vehicle Mass
g =   vehicle_data.vehicle.g  ;   % [m/s^2] Gravitational acceleration


t0 = 0;
% tf = 10;
N = 1e3+1;
time = linspace(t0,t_f,N);
dt = time(2)-time(1);

% 
toDeg = 180/pi;
toRad = 180/pi;


%% Extract and interpolate data from simulink model

time_sim = input_sim.delta_D.time;

% Input
delta_D = interp1(time_sim,input_sim.delta_D.data,time);
Tw_rr   = interp1(time_sim,input_sim.Tw_rr.data,time);
Tw_rl   = interp1(time_sim,input_sim.Tw_rl.data,time);
Tw_fr   = interp1(time_sim,input_sim.Tw_fr.data,time);
Tw_fl   = interp1(time_sim,input_sim.Tw_fl.data,time);

figure('Name','Inputs','NumberTitle','off'),
% Steering angle delta_D
subplot(211)
plot(time,delta_D*180/pi)
grid on
title('Steering Angle','LineWidth',2)
ylabel('Angle [deg]')

% Torque at wheels
subplot(212)
hold on
plot(time,Tw_rr,'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,Tw_rl,'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,Tw_fr,'Color',color('orange'),'LineWidth',2)
plot(time,Tw_fl,'Color',color('red'),'LineStyle','--','LineWidth',2)
legend({'$Tw_{rr}$','$Tw_{rl}$','$Tw_{fr}$','$Tw_{fl}$'})
grid on
title('Torque at wheels')
xlabel('Time [s]')
ylabel('Torque [Nm]')

%% Filtered Input

delta_filt = interp1(time_sim,filtered_input.delta.data,time);
Tw_rr_filt = interp1(time_sim,filtered_input.Tw.data(:,1),time);
Tw_rl_filt = interp1(time_sim,filtered_input.Tw.data(:,2),time);
Tw_fr_filt = interp1(time_sim,filtered_input.Tw.data(:,3),time);
Tw_fl_filt = interp1(time_sim,filtered_input.Tw.data(:,4),time);

figure('Name','Filtered Inputs','NumberTitle','off'),
% Steering angle delta_D
subplot(211)
plot(time,delta_filt*180/pi)
grid on
title('Steering Angle','LineWidth',2)
ylabel('Angle [deg]')

% Torque at wheels
subplot(212)
hold on
plot(time,Tw_rr_filt,'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,Tw_rl_filt,'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,Tw_fr_filt,'Color',color('orange'),'LineWidth',2)
plot(time,Tw_fl_filt,'Color',color('red'),'LineStyle','--','LineWidth',2)
legend({'$Tw_{rr}$','$Tw_{rl}$','$Tw_{fr}$','$Tw_{fl}$'})
grid on
title('Torque at wheels')
xlabel('Time [s]')
ylabel('Torque [Nm]')

% Vehicle Measurements
u       = interp1(time_sim,vehicle.u_t___m_s_.data,time);
v       = interp1(time_sim,vehicle.v_t___m_s_.data,time);
Omega   = interp1(time_sim,vehicle.Omega_t___rad_s_.data,time);
phi     = interp1(time_sim,vehicle.phi_t___rad_.data,time);
delta   = interp1(time_sim,vehicle.delta_t___rad_.data,time);
alpha   = interp1(time_sim,vehicle.alpha.data,time);
kappa   = interp1(time_sim,vehicle.kappa.data,time);
Fz      = interp1(time_sim,vehicle.Fz.data,time);
phi_dot = interp1(time_sim,vehicle.phi_dot_t___rad_s_.data,time);
omega   = interp1(time_sim,vehicle.omega.data,time);

% % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
figure('Name','Vehicle Meas.','NumberTitle','off'),
% --- u --- %
ax(1) = subplot(321);
plot(time,u,'LineWidth',2)
grid on
title('$u$ [m/s]')
% --- v --- %
ax(2) = subplot(322);
plot(time,v,'LineWidth',2)
grid on
title('$v$ [m/s]')
% --- beta --- %
ax(3) = subplot(322);
plot(time,atan2(v,u)*toDeg,'LineWidth',2)
grid on
title('$\beta$ [deg]')
% --- Omega --- %
ax(4) = subplot(323);
plot(time,Omega*toDeg,'LineWidth',2)
grid on
title('$\Omega$ [deg/s]')
% --- phi --- %
ax(5) = subplot(324);
plot(time,phi*toDeg,'LineWidth',2)
grid on
title('$\phi$ [deg]')
% --- delta --- %
ax(6) = subplot(325);
plot(time,delta*toDeg,'LineWidth',2)
grid on
title('$\delta$ [deg]')
% --- phi_dot --- %
ax(7) = subplot(326);
plot(time,phi_dot*toDeg,'LineWidth',2)
grid on
title('$\dot{\phi}$ [deg/s]')

linkaxes(ax,'x')
clear ax


%% Steering behaviour: fixed steering angle
VG  = sqrt(v.^2+u.^2); % absolute velocity
R   = VG./Omega;       % curvature radius
ay  = VG.*Omega;       % lateral acceleration
idx = find(time > Tstart_US_test);

figure('Name','Steering behaviour','NumberTitle','off'),
subplot(3,1,1)
plot(ay(idx),R(idx),'LineWidth',2)
hold on
plot(ay(idx),1/curv_akermann*ones(size(idx)),'--r','LineWidth',2)
xlabel('$a_y [\frac{m}{s^2}]$')
ylabel('$R [m]$')

subplot(3,1,2)
plot(VG(idx),R(idx),'LineWidth',2)
hold on
plot(VG(idx),1/curv_akermann*ones(size(idx)),'--r','LineWidth',2)
xlabel('$V_G [\frac{m}{s}]$')
ylabel('$R [m]$') 

subplot(3,1,3)
plot(VG(idx),Omega(idx)./delta(idx),'LineWidth',2)
xlabel('$V_G [\frac{m}{s}]$')
ylabel('$\frac{\Omega}{\delta} [\frac{1}{s}]$') 
title('Yaw gain')


%%
% % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
figure('Name','Wheel Meas.','NumberTitle','off'),
% --- alpha --- %
ax(1) = subplot(221);
hold on
plot(time,alpha(:,1),'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,alpha(:,2),'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,alpha(:,3),'Color',color('orange'),'LineWidth',2)
plot(time,alpha(:,4),'Color',color('red'),'LineStyle','--','LineWidth',2)
grid on
title('Lateral Slip')
legend({'$\alpha_{rr}$','$\alpha_{rl}$','$\alpha_{fr}$','$\alpha_{fl}$'})
% --- kappa --- %
ax(2) = subplot(222);
hold on
plot(time,kappa(:,1),'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,kappa(:,2),'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,kappa(:,3),'Color',color('orange'),'LineWidth',2)
plot(time,kappa(:,4),'Color',color('red'),'LineStyle','--','LineWidth',2)
% ylim([-2 +2])
YL = get(gca,'YLim');
Vlow = 15;
Vth = 5;
indx0 = find(abs(u)<=0.1);
indxVlow = find(abs(u)>=Vlow);
indxVth  = find(abs(u)>=Vth);
% plot([time(indx0(1)) time(indx0(1))],YL,'Color','m')
% plot([time(indxVlow(1)) time(indxVlow(1))],YL,'Color','r')
% plot([time(indxVth(1)) time(indxVth(1))],YL,'Color','b')
grid on
title('Longitudinal Slip')
legend({'$\kappa_{rr}$','$\kappa_{rl}$','$\kappa_{fr}$','$\kappa_{fl}$'})
% --- Fz --- %
ax(3) = subplot(223);
hold on
plot(time,Fz(:,1),'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,Fz(:,2),'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,Fz(:,3),'Color',color('orange'),'LineWidth',2)
plot(time,Fz(:,4),'Color',color('red'),'LineStyle','--','LineWidth',2)
grid on
title('Vertical Load')
legend({'$Fz_{rr}$','$Fz_{rl}$','$Fz_{fr}$','$Fz_{fl}$'})
% --- omega --- %
ax(4) = subplot(224);
hold on
plot(time,omega(:,1),'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,omega(:,2),'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,omega(:,3),'Color',color('orange'),'LineWidth',2)
plot(time,omega(:,4),'Color',color('red'),'LineStyle','--','LineWidth',2)
grid on
title('Wheel Speed')
legend({'$\omega_{rr}$','$\omega_{rl}$','$\omega_{fr}$','$\omega_{fl}$'})

linkaxes(ax,'x')

clear ax

%% Forces
Fx_rr = interp1(time_sim,force_rr.data(:,1),time);
Fx_rl = interp1(time_sim,force_rl.data(:,1),time);
Fx_fr = interp1(time_sim,force_fr.data(:,1),time);
Fx_fl = interp1(time_sim,force_fl.data(:,1),time);

Fy_rr = interp1(time_sim,force_rr.data(:,2),time);
Fy_rl = interp1(time_sim,force_rl.data(:,2),time);
Fy_fr = interp1(time_sim,force_fr.data(:,2),time);
Fy_fl = interp1(time_sim,force_fl.data(:,2),time);

figure('Name','Forces','NumberTitle','off')
ax(1) = subplot(211);
hold on
plot(time,Fx_rr,'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,Fx_rl,'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,Fx_fr,'Color',color('orange'),'LineWidth',2)
plot(time,Fx_fl,'Color',color('red'),'LineStyle','--','LineWidth',2)
grid on
legend({'$Fx_{rr}$','$Fx_{rl}$','$Fx_{fr}$','$Fx_{fl}$'})
title('Longitudinal Forces')
xlim([t0 t_f])

ax(2) = subplot(212);
hold on
plot(time,Fy_rr,'Color',color('deepsky_blue'),'LineWidth',2)
plot(time,Fy_rl,'Color',color('blue'),'LineStyle','--','LineWidth',2)
plot(time,Fy_fr,'Color',color('orange'),'LineWidth',2)
plot(time,Fy_fl,'Color',color('red'),'LineStyle','--','LineWidth',2)
grid on
legend({'$Fy_{rr}$','$Fy_{rl}$','$Fy_{fr}$','$Fy_{fl}$'})
title('Lateral Forces')

linkaxes(ax,'x')
clear ax



%% Plot trajectory with animation
% Trajectory
theta = interp1(time_sim,trajectory.data(:,1),time);
x     = interp1(time_sim,trajectory.data(:,2),time);
y     = interp1(time_sim,trajectory.data(:,3),time);

fprintf('Press enter to continue...\n')
pause

figure,
xlim([min(x)-5 max(x)+5])
ylim([min(y)-5 max(y)+5])
axis equal
grid on
for i = 1:N
  hold on
  rot_mat = [cos(theta(i)) -sin(theta(i)) ; sin(theta(i)) cos(theta(i))];
  pos_rr = rot_mat*[-Lr -Wr]';
  pos_rl = rot_mat*[-Lr +Wr]';
  pos_fr = rot_mat*[+Lf -Wf]';
  pos_fl = rot_mat*[+Lf +Wf]';
  pos = [pos_rr pos_rl pos_fl pos_fr];
  p = patch(x(i) + pos(1,:),y(i) + pos(2,:),'blue');
  plot(x(1:i),y(1:i),'Color',color('red'),'LineWidth',2,'LineStyle',':')
  pause(dt)
  if(i ~= N)
    delete(p)
  end
end

%%
figure('Name','Trajectory','NumberTitle','off'),
plot(x,y)
grid on
xlabel('x-coord [m]')
ylabel('y-coord [m]')
title('Trajectory')
axis equal







