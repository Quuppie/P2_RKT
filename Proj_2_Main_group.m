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
% How can we 'mesh' together two variables to plot as a surf? That seems to
% be the last step, but it's complex.

%% Common practices
clear;
clc;
close all;
format long;
%% scale values

Design_space_dimension = 3;

delta_theta = deg2rad(0.5); % radians
delta_p_i = 3 * 6894.76; % pascals
delta_w_i = 0.0001; % m^3

gradient_scalar = [2/10000, 10000000, 2/10000];

delta_vals = [delta_theta, delta_p_i, delta_w_i];

const = getConst();
runs = 200;

%help input

for i = 1:runs

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

% Commented out to test for imaginary numbers
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

figure()
plot(test_pressure,pressure_dist_chart)
title 'Evolution of Pressure with Time'

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



%% NESTING TEST (For Surf)
%Theta is the outer (vertical), pressure inner (horizontal), nxn
n = 20

%Theta versus pressure
for i = 1:n
    test_theta_surf = linspace(deg2rad(1), deg2rad(80), n);
    test_pressure_surf = linspace(const.p_amb + 1, 80 * 6874, n);
        for j = 1:n
         test_const = const;
         test_const.theta_i = test_theta_surf(i);
         test_const.p_r_i = test_pressure_surf(j);

         %Have defined both, for i and j

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
        surf_dist_chart_tp(i,j) = max(state(:,1));


        end 
end

%Now water versus theta!
n = 20

for i = 1:n

    test_theta_surf = linspace(deg2rad(1), deg2rad(80), n);
    test_water_surf = linspace(0+0.0000001, const.Vol_bottle, n);
        for j = 1:n
         test_const = const;
         test_const.theta_i = test_theta_surf(i);
         test_const.Vol_w_i = test_water_surf(j);

         %Have defined both, for i and j

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
        surf_dist_chart_tw(i,j) = max(state(:,1));


        end 
end

%Water and Pressure
for i = 1:n

    test_water_surf = linspace(0+0.0000001, const.Vol_bottle, n);
    test_pressure_surf = linspace(const.p_amb + 1, 80 * 6874, n);

        for j = 1:n
         test_const = const;
         test_const.Vol_w_i = test_water_surf(i);
         test_const.p_r_i = test_pressure_surf(j);

         %Have defined both, for i and j

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
        surf_dist_chart_wp(i,j) = max(state(:,1));


        end 
end

%% Plotting the Surface!
%Theta vs Pressure
figure()
surf(real(test_theta_surf),real(test_pressure_surf),real(surf_dist_chart_tp));
title 'Testing Shenanigans'
xlabel 'Theta'
ylabel 'Pressure'
zlabel 'Distance of Rocket'

%Theta vs Water - THE COOL ONE!!
figure()
surf(rad2deg(real(test_theta_surf)),real(test_water_surf)*1000,real(surf_dist_chart_tw),'FaceAlpha',0.9);
view([45 22])
Z = ones(size(surf_dist_chart_tw))*85;
hold on
% patch([0,0,5,5],[-1,-1,1,1],[85,85,85,85],'w','FaceAlpha',0.7);
% surf(real(test_theta_surf),real(test_water_surf),Z);
title 'A 3D Slice of the Design Space: Water-Theta'
xlabel 'Theta (rad)'
ylabel 'Water Mass'
zlabel 'Distance of Rocket'

print ('topographic_theta_water', '-dpng', '-r300')

%Water vs Pressure
figure()
surf(real(test_water_surf),real(test_pressure_surf),real(surf_dist_chart_wp));
title 'Testing Shenanigans'
xlabel 'Water'
ylabel 'Pressure'
zlabel 'Distance of Rocket'



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

%Plotting Relationships

figure()
plot(test_water,water_dist_chart)
grid on
xlabel 'Mass of Water (kg?)'
ylabel 'Distance Traveled'
title 'Variation of Distance with Water Mass'

figure()
grid on
plot(rad2deg(test_theta),theta_dist_chart)
xlabel 'Launch Angle (deg)'
ylabel 'Distance Traveled'
title 'Variation of Distance with Launch Angle'

% surf_dist = [water_dist_chart',theta_dist_chart'];
% surf(test_water',test_theta',surf_dist);
% title 'Shenanigans'

const.p_r_i / 6894.76;

figure()
plot(timevector, step.water)
title 'Evolution of Water Mass over Time'

%% Finding the Force Values with Thrust Function that mimics OdeFun
% Thrust() is a modified OdeFun that returns thrust and phase rather than
% state changes.
F_Thrust = zeros(length(t),1);
for i = 1:length(t)
    [F_thrust(i),phase(i)] = Thrust(state(i,:),(const),m_air_i);
end
%% Finding answers
answers.MaxThrust = max(F_thrust)
answers.MaxAltitude = max(state(:,3));
answers.MaxDistance = max(state(:,1))
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

