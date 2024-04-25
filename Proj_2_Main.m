%{
Author: Alden Roswell
Assignment: Project 2 indidvidual submission
Creation Date: 4/15/24
Inputs: getConst.m OdeFun.m Thrust.m
Outputs: Flight path of bottle rocket
Purpose: to predict different times of a bottle rocket with distinct
phases of flight
%}

%% Notes
% Any variable ending in "_i" is an initial condition

%% Common practices
clear;
clc;
close all;

%% Call Constant Function
const = getConst();

%% Determing initial conditions based on constants
Vol_air_i = const.Vol_bottle - const.Vol_w_i; %
m_air_i = (Vol_air_i * const.p_r_i)/(const.R_air * const.T_i); % mass of air at launch.
m_r_i = const.m_bottle + const.row_w * const.Vol_w_i + m_air_i; % mass of rocket at launch;
timespan = [0,5];

%% creating initial state vector From initial Conditions.
state_i = [const.x_i; 0; const.z_i; 0; m_r_i; Vol_air_i; m_air_i];

%% Running ODE45
[t,state] = ode45(@(t,state) OdeFun(t,state,const,m_air_i), timespan, state_i);

%% Finding the Force Values with Thrust Function that mimics OdeFun
% Thrust() is a modified OdeFun that returns thrust and phase rather than
% state changes.
F_Thrust = zeros(length(t),1);
for i = 1:length(t)
    [F_thrust(i),phase(i)] = Thrust(state(i,:),(const),m_air_i);
end
%% Finding answers
answers.MaxThrust = max(F_thrust);
answers.MaxAltitude = max(state(:,3));
answers.MaxDistance = max(state(:,1));
answers

%% Plot Thrust vs Time
figure()
title('Thrust Over Time')
hold on;
plot(t,F_thrust(:), 'LineWidth',1)
xline(t(length(phase(phase == 1))), 'color', 'green', 'LineWidth',1);
xline(t(length(phase(phase <= 2))), 'color', 'yellow', 'LineWidth',1);
xlim([0,.5])
grid("on")
xlabel('Time (s)')
ylabel('Thrust (N)')
legend('Thrust (N) Over Time','End of Phase 1', 'End of Phase 2')


%% Plot Trajectory
figure()
hold on;
title('Rocket Trajectory')
plot(state(:,1),state(:,3))
legend('Rocket Position')
xlabel('Distance (m)')
ylabel('Altitude (m)')
ylim([0, 20])
grid on;

