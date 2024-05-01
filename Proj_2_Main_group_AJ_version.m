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
format long;
%% scale values

%% Defining design space dimension for Gradient ascent.
Design_space_dimension = 3;

delta_theta = deg2rad(0.5); % radians
delta_p_i = 1 * 6894.76; % pascals
delta_w_i = 0.0001; % m^3

gradient_scalar = [1/10000, 10000000, 1/10000];

delta_vals = [delta_theta, delta_p_i, delta_w_i];

const = getConst();

for i = 1:200

    delta_matrix = zeros(Design_space_dimension);
    const_vector = [const.theta_i;const.p_r_i;const.Vol_w_i];

    for k = 1:Design_space_dimension
        delta_matrix(k,k) = delta_vals(k);
    end

    for j = 1:Design_space_dimension + 1
        %% Change Constant Function

        if j > 1
            const_vector = const_vector + delta_matrix(:,j-1);
        end

        test_const.theta_i = const.theta_i;
        test_const.p_r_i = const.p_r_i;
        test_const.Vol_w_i = const.Vol_bottle;
        test_const = const;
        temporary_delta_vals = delta_vals;
        %
        %  if const_vector(1) > 0 && const_vector(1) < deg2rad(80)
        %  test_const.theta_i = const_vector(1);
        % % temporary_delta_vals(1,i) = delta_vals(1);
        %  end
        %
        %  if const_vector(2) > (const.p_amb) && const_vector(2) < (85 * 6894.76)
        %  test_const.p_r_i = const_vector(2);
        %  %temporary_delta_vals(2,i) = delta_vals(2);
        %  end
        %
        %  if const_vector(3) > 0 && const_vector(3) < const.Vol_bottle
        %  test_const.Vol_w_i = const_vector(3);
        %  %temporary_delta_vals(3,i) = delta_vals(3);
        %  end
        %
        step.theta(i) = test_const.theta_i;
        step.pressure(i) = test_const.p_r_i;
        step.water(i) = test_const.Vol_w_i;

        test_const.theta_i = const_vector(1);
        test_const.p_r_i = const_vector(2);
        test_const.Vol_w_i = const_vector(3);


        %% Determing initial conditions based on constants
        Vol_air_i = test_const.Vol_bottle - test_const.Vol_w_i; %
        m_air_i = (Vol_air_i * test_const.p_r_i)/(test_const.R_air * test_const.T_i); % mass of air at launch.
        m_r_i = test_const.m_bottle + test_const.row_w * test_const.Vol_w_i + m_air_i; % mass of rocket at launch;
        timespan = [0,5];
        %% creating initial state vector From initial Conditions.
        state_i = [test_const.x_i; 0; test_const.z_i; 0; m_r_i; Vol_air_i; m_air_i];
        %% Running ODE45
        state = [];
        [t,state] = ode45(@(t,state) OdeFun(t,state,test_const,m_air_i), timespan, state_i);
        distance_test(i,j) = max(state(:,1));

    end

    for l = 1:Design_space_dimension

        delta_dist(i,l) = ((distance_test(i,l+1) - distance_test(i,1)) / temporary_delta_vals(l));
    end

    update_vector = delta_dist(i,:) .* gradient_scalar;

    const.theta_i = max([min([const.theta_i + update_vector(1), deg2rad(80)]), 0]);
    const.p_r_i = max([min([const.p_r_i + update_vector(2), 85 * 6894.76]), const.p_amb + 1000 ]);
    const.Vol_w_i = max([min([const.Vol_w_i + .000001 * update_vector(3), const.Vol_bottle]), 0]);

end


for i = 1:100

    test_theta = linspace(deg2rad(1), deg2rad(80), 100);
    test_const = const;
    test_const.theta_i = test_theta(i);



    %% Determing initial conditions based on constants
    Vol_air_i = test_const.Vol_bottle - test_const.Vol_w_i; %
    m_air_i = (Vol_air_i * test_const.p_r_i)/(test_const.R_air * test_const.T_i); % mass of air at launch.
    m_r_i = test_const.m_bottle + test_const.row_w * test_const.Vol_w_i + m_air_i; % mass of rocket at launch;
    timespan = [0,5];
    %% creating initial state vector From initial Conditions.
    state_i = [test_const.x_i; 0; test_const.z_i; 0; m_r_i; Vol_air_i; m_air_i];
    %% Running ODE45
    state = [];
    [t,state] = ode45(@(t,state) OdeFun(t,state,test_const,m_air_i), timespan, state_i);
    theta_dist_chart(i) = max(state(:,1));

end
figure()
plot(rad2deg(test_theta),theta_dist_chart)
%% Commented out to test for imaginary numbers
for i = 1:100

    test_pressure = linspace(const.p_amb + 1, 80 * 6874, 100);
    test_const = const;
    test_const.p_r_i = test_pressure(i);



    %% Determing initial conditions based on constants
    Vol_air_i = test_const.Vol_bottle - test_const.Vol_w_i; %
    m_air_i = (Vol_air_i * test_const.p_r_i)/(test_const.R_air * test_const.T_i); % mass of air at launch.
    m_r_i = test_const.m_bottle + test_const.row_w * test_const.Vol_w_i + m_air_i; % mass of rocket at launch;
    timespan = [0,5];
    %% creating initial state vector From initial Conditions.
    state_i = [test_const.x_i; 0; test_const.z_i; 0; m_r_i; Vol_air_i; m_air_i];
    %% Running ODE45
    state = [];
    [t,state] = ode45(@(t,state) OdeFun(t,state,test_const,m_air_i), timespan, state_i);
    pressure_dist_chart(i) = max(state(:,1));

end
target_distance = 85;
target_pressure = interp1(real(pressure_dist_chart),test_pressure, target_distance);


figure()
plot(test_pressure,pressure_dist_chart)
%% Finds desired pressure
xline(target_pressure)
yline(target_distance)
for i = 1:100

    test_water = linspace(0+0.0000001, const.Vol_bottle, 100);
    test_const = const;
    test_const.Vol_w_i = test_water(i);



    %% Determing initial conditions based on constants
    Vol_air_i = test_const.Vol_bottle - test_const.Vol_w_i; %
    m_air_i = (Vol_air_i * test_const.p_r_i)/(test_const.R_air * test_const.T_i); % mass of air at launch.
    m_r_i = test_const.m_bottle + test_const.row_w * test_const.Vol_w_i + m_air_i; % mass of rocket at launch;
    timespan = [0,5];
    %% creating initial state vector From initial Conditions.
    state_i = [test_const.x_i; 0; test_const.z_i; 0; m_r_i; Vol_air_i; m_air_i];
    %% Running ODE45
    state = [];
    [t,state] = ode45(@(t,state) OdeFun(t,state,test_const,m_air_i), timespan, state_i);
    water_dist_chart(i) = max(state(:,1));

end
figure()
plot(test_water,water_dist_chart)


test_const = const;
test_const.p_r_i = target_pressure;

Vol_air_i = test_const.Vol_bottle - test_const.Vol_w_i; %
m_air_i = (Vol_air_i * test_const.p_r_i)/(test_const.R_air * test_const.T_i); % mass of air at launch.
m_r_i = test_const.m_bottle + test_const.row_w * test_const.Vol_w_i + m_air_i; % mass of rocket at launch;
timespan = [0,5];
%% creating initial state vector From initial Conditions.
state_i = [test_const.x_i; 0; test_const.z_i; 0; m_r_i; Vol_air_i; m_air_i];
%% Running ODE45
state = [];
[t,state] = ode45(@(t,state) OdeFun(t,state,test_const,m_air_i), timespan, state_i);
max_dist_land = max(state(:,1))


timevector = 1:length(distance_test);

figure()
hold on;
title('Delta Parameters')
subplot(3,1,1)
plot(timevector,delta_dist(:,1))
title('Launch Angle')

subplot(3,1,2)
plot(timevector,delta_dist(:,2))
title('Pressure of Air')

subplot(3,1,3)
plot(timevector,delta_dist(:,3))
title('Water Mass')

const.p_r_i / 6894.76;

figure()
plot(timevector, step.water)

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
answers;

%% Plot Thrust vs Time
figure()
title('Thrust Over Time')
hold on;
plot(t,F_thrust(:), 'LineWidth',1)
%xline(t(length(phase(phase == 1))), 'color', 'green', 'LineWidth',1);
%xline(t(length(phase(phase <= 2))), 'color', 'yellow', 'LineWidth',1);
grid("on")
xlabel('Time (s)')
ylabel('Thrust (N)')


%% Plot Trajectory
figure()
hold on;
title('Rocket Trajectory')
plot(state(:,1),state(:,3))

xlabel('Distance (m)')
ylabel('Altitude (m)')
ylim([0,40])
grid on;

