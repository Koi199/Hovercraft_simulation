% MECH 223 Hovercraft- Simulation II
% Team 14 - BEST TEAM BABY

% This script uses a for-loop to obtain position, velocity, and
% acceleration for the hovercraft as well as the lift and plenum 
%One lap is 17.34 m
%Turn radius is 1.2 m
%end to end -> 7.3 m 

%% Defining Constants

rho = 1.225; %still Air
Cl = 0.3; 
std_Cl = 0.1*Cl;
C_fd = 0.15; %Friction drag
C_pd = 1.48;%Cd from FLow Analysis (Pressure)
Cd = C_fd+C_pd;
std_Cd = 0.1*Cd;
c = 10/1000;
std_c = 0.0015/2;
g = 9.81;
delta_t = 0.1;
final_t = 15;
% Diameter of the propeller in inches
D = 6.5;
R = (D*0.0254)/2;
%Pitch is the forward distance that a propeller would theoretically travel 
% in a single rotation if there were no slip present – imagine a screw being
% driven into a piece of wood (inches)
pitch = 2;

%Frontal Area
L = 0.35;
W = 0.25;
H = 0.04;
Area = pi/4*((D*0.0254/2)^2) + (W*H);

%Number of blades for the propeller
N = 2;

%% Mass Parameters

%Number of Batteries
NUMbat = 12;
m_batteries = NUMbat*0.02;
m_electronics = 0.15;
m_materials = 0.31;
m = m_batteries + m_electronics+m_materials;

%% lift and Plenum Calculations

% Other useful constants
K = 0.6;
hp = 0.0001;
% Underneath Area 
% Dimensions in inches
L_p = 12;
W_p = 6.25;
Hyp_p = 4.25;
H_p = sqrt(Hyp_p^2-(W_p/2)^2);
L_pm = L_p*0.0254;
W_pm = W_p*0.0254;
H_pm = H_p*0.0254;
Area_Plenum = L_pm*W_pm + H_pm*W_pm*0.5;

% Blower Rating
Qfmax = 0.453069552/60;
Pfmax = 1*g/Area_Plenum; %Estimated value for the maximum pressure, 
% We stacked a few wooden blocks to check lift against Pressure. 

FUN1 = @(Qf) Pfmax*(1-Qf/Qfmax)-(Qf/(K*L*hp))^2*rho/2;
Qf_op = fzero(FUN1, 0.0001);
Pressure_needed = (m*g)/Area_Plenum;
% Qf = Qfmax*(1 - Pf/Pfmax);
h_surface = Qf_op/(K*L*sqrt(2*Pressure_needed/rho));

%% Determining Ideal RPM
% Motor_curve eqn
% 11.4/1000-((11.4/1000)/5190).*RPM
% Torque = (1/8)*rho.*Rad_s.^2.*Cd*c*R^4;
% Torque = (1/8)*rho*(RPM*2*pi/60)^2*Cd*c*R^4;
% Thrust = (1/6)*rho*Rad_s^2*Cl*c*R^3;
% Pressure_drag = (1/2)*Vel^2*rho*Area*C_pd;
% Motor_curve = 11.4/1000-((11.4/1000)/5190)*RPM;

% plot(RPM, Torque)
% hold on
% plot(RPM, Motor_curve)

% d(Torque)/d(C_fd) = (1/8)*rho*(RPM*2*pi/60)^2*c*R^4;
% d(Torque)/d(c) = (1/8)*rho*(RPM*2*pi/60)^2*C_fd*R^4;

FUN = @(RPM) (N*(1/8)*rho*(RPM*2*pi/60)^2*C_fd*c*R^4)-(11.4/1000-((11.4/1000)/5190)*RPM);
RPM_op = fzero(FUN, 3000);
%RPM of the Propeller
RPM = RPM_op;

%% Iteration for Velocity, Acceleration and Distance

% Defining arrays
t = 0:delta_t:final_t;
a = zeros(size(t));
v = zeros(size(t));
x = zeros(size(t));

% Thrust = 1.225*pi*((0.0254*D)^2/4)*(((RPM*0.0254*pitch/60)^2)-(RPM*0.0254*pitch ...
%     /60).*VelData)*(D/(pitch*3.29546))^1.5;
% Drag = 0.5*Cd*rho*VelData.^2*A;

% Use a for loop to find corresponding values for each quantity with respect
% to time
for t_int = 1:length(t)-1
   
   % Determines velocity
   v(t_int+1) = v(t_int) + a(t_int)*delta_t;
   % Determines displacement
   x(t_int+1) = x(t_int) +v(t_int)*delta_t+0.5*a(t_int)*delta_t^2;
   % Determines acceleration
   a(t_int+1) = ((1.225*pi*((0.0254*D)^2/4)*(((RPM*0.0254*pitch/60)^2)-(RPM*0.0254*pitch ...
    /60)*v(t_int))*(D/(pitch*3.29546))^1.5)-(0.5*(Cd)*rho*v(t_int)^2*Area))/m;

end

% Plot relevant graphs
plot(t, v);
hold on
plot(t, x);
plot(t, a);
xlabel('Time (s)');
title('Simulation - Project II');
legend('Velocity (m/s)','Position (m)','Acceleration (m/s^2)');
grid on  

% Check the final values
disp(v(end));
disp(x(end));
disp(a(end));

%% Uncertainty in Thrust and Drag
Thrust = N*(1/6)*rho*(RPM*2*pi/60)^2*Cl*c*R^3;
Drag = N*(1/6)*rho*(RPM*2*pi/60)^2*Cd*c*R^3;

% At 95% confidence interval
std_Thrust = (1/6)*(rho)*(RPM*2*pi/60)^2*R^3*sqrt(c^2*(std_c)^2+Cl^2*(std_Cl)^2);
std_Drag = (1/6)*(rho)*(RPM*2*pi/60)^2*R^3*sqrt(c^2*(std_c)^2+Cl^2*(std_Cd)^2);
std_acc = a(end)*sqrt((std_Thrust/Thrust)^2+(std_Drag/Drag)^2);
std_Vel = v(end)*sqrt((std_acc/a(end))^2);
std_dis = x(end)*sqrt((std_acc/a(end))^2+(std_Vel/v(end))^2);

disp(['Thrust(N) = ' num2str(Thrust) ' +/- ' num2str(std_Thrust)])
disp(['Drag(N) = ' num2str(Drag) ' +/- ' num2str(std_Drag)])
disp(['Terminal Vel(m/s) = ' num2str(v(end)) ' +/- ' num2str(std_Vel)])
disp(['Distance (m) after ' num2str(final_t) 's = ' num2str(x(end)) ' +/- ' num2str(std_dis)])
disp(['Acceleration (m/s^2) after ' num2str(final_t) 's = ' num2str(a(end)) ' +/- ' num2str(std_acc)])