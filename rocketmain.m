%{
Author: Cayden Stratford    
Assignment: ASEN 2012 - Project 02
Creation Date: 04.19.2024
Inputs: A Bottle Rocket's Governing Equations w/ Data
Outputs: Position and Velocity with Direction and Magnitude
Purpose: To explore the very boundaries of our numerical integration
capability

Current issue: Submitting!
%}

clear;
close all;
clc

load project2verification.mat
%V means velocity, v means volume.

%Inputs
c = getconst();
% 
% i.pump_P0 = 358527;     %Thermodynamics (Payload)
% i.water_v0 = 0.001; %Given
% i.water_m0 = c.rho_water*i.water_v0;   %Simple calc
% i.theta0 = 0.7*pi;   %Launch Parameters (Radians)
% i.Cd = 0.5; %Buh??

%Defining Realm
t = (0:0.1:5);
t0 = t(1);
tf = t(length(t));
Vz = 0;
Vx = 0;
Z = 0;
X = 0;
v0_air = c.i.v0_air;
M0_air = c.i.m0_air;
M0_rkt = c.m_empty+c.rho_water*(c.v_empty-v0_air)+((c.i.P0_air*v0_air)/(c.R_air*c.T0)); %I'm guessing that this flow is largely incompressible - rho same
M_rkt = M0_rkt; M_air = M0_air; v_air = v0_air;

%Preallocating... X is Statevector containing V and H
X = [X Vx Z Vz M_rkt v_air M_air]';
X0 = [0 0 0.25 0.00000 M0_rkt v0_air M0_air]';    %V0 = 0; x0 = 0; z0 = 0.25

%% ODE45 Call and Outputs
[time, ode_out] = ode45(@(t,X) dy_fun(t,X,c), [t0 tf], X0);

%Very Bad Magic Numbers
% for i = 1:height(ode_out)
% if ode_out(i,3) < 0
%     ode_out(i,3) = 0;
% end
% end

%CALCULATING THRUST
for i = 1:numel(time)
     [~,thrust(i,:)] = dy_fun(time(i),ode_out(i,:),c);
end
for i = 1:numel(time)
     site_Distance = sqrt(ode_out(i,1)^2+ode_out(i,3)^2);
end


maxThrust = max(thrust)
maxHeight = max(ode_out(:,3))
maxDistance = max(site_Distance)


max_dist_index = find(ode_out(:,1)>maxDistance,1);

fthrust = figure();
plot(time,thrust(:,1))
hold on
plot(verification.time,verification.thrust,'Color',"#A2142F")
title 'Thrust as a function of time'
xlabel 'Time (s)'
ylabel 'Thrust (n)'

fvol = figure();
plot(time,ode_out(:,6))
hold on
plot(verification.time,verification.volume_air, ':k')
title 'Volume of Air as a function of time'
xlabel 'Time (s)'
ylabel 'Thrust (n)'
legend ('Calculated Thrust', 'Verification Thrust');

%plots - Trajectory (x vs z), and Thrust curve (T, t)
ftest = figure();
plot(ode_out(:,1),ode_out(:,3));
hold on
plot(verification.distance,verification.height,'Color',"#A2142F")
legend 'Distance'
title 'Trajectory'
xlabel 'Drift (m)'
ylabel 'Altitude (m)'
legend ('Calculated Trajectory', 'Verification Trajectory');
print('Plot_Trajectory','-dpng','-r300')

fplots = figure();
subplot(2,2,1)
plot(time, ode_out(:,1))
hold on
plot(verification.time,verification.distance,'Color',"#A2142F")
title 'X Position'
xlabel 'Time (s)'
ylabel 'Drift (m)'
text(time(max_dist_index),ode_out(max_dist_index,1)-2,'|','Color','k','FontSize',5);
text(time(max_dist_index),ode_out(max_dist_index,1)+2,'Oof','Color','k','FontSize',6);
hold off
subplot(2,2,2)
plot(time, ode_out(:,2))
hold on
plot(verification.time,verification.velocity_x,'Color',"#A2142F")
title 'X Velocity'
xlabel 'Time (s)'
ylabel 'v (m/s)'
legend('Calculated','Verification','Location','Northeast')
hold off
subplot(2,2,3)
plot(time, ode_out(:,3))
hold on
plot(verification.time,verification.height,'Color',"#A2142F")
title 'Z Position'
xlabel 'Time (s)'
ylabel 'altitude (m)'
hold off
subplot(2,2,4)
plot(time, ode_out(:,4));
hold on
plot(verification.time,verification.velocity_y,'Color',"#A2142F")
title 'Z Velocity'
xlabel 'Time (s)'
ylabel 'v (m/s)'
hold off
print('Plot_Kinematic','-dpng','-r300')

fpropellant = figure;
subplot(3,1,1)
plot(time, ode_out(:,5));
hold on
title 'Overall Mass of Rocket'
xlabel 'Time (s)'
ylabel 'Mass (kg)'
subplot(3,1,2)
plot(time, ode_out(:,6));
title 'Volume of air Remaining'
xlabel 'Time (s)'
ylabel 'Volume (m^3)'
subplot(3,1,3)
plot(time, ode_out(:,7));
title 'Mass of Air Remaining'
xlabel 'Time (s)'
ylabel 'mass (kg)'
print('Plot_Propulsion','-dpng','-r300')

%% Functions
function [dy_out, telem] = dy_fun(t,X,c)
%Function to calculate net force for each phase of flight
%All are column vectors, [x, z]
%For notation, refer to initial state vector call   X = [X Vx Z Vz M_rkt v_air M_air]';
x = X(1);   %X value? 'x1'
Vx = X(2);    
z = X(3);
Vz = X(4);
M_rkt = X(5);
v_air = X(6);
M_air = X(7);

%Direction and Velocity
velocity = [Vx; Vz]; %[Vx, Vz]
%mag_v = sqrt((velocity(1)^2+velocity(2)^2));
mag_v = sqrt(velocity(1)^2+velocity(2)^2);
hbar = [(Vx/mag_v);(Vz/mag_v)];

%Remaining Resources
M_water = M_rkt - M_air - c.m_empty;
P_end = 0;

A_throat = (pi/4)*(c.d_throat^2);
A_rkt = (pi/4)*(c.d_bottle^2);

%Relevant equations - All exclude velocity (?)
%fexhaustvel = @(pressure,ambient,rhowater) sqrt((2*(pressure-ambient))/rhowater);    %Equation 8
%fpressure = @(P0_air,v0_air,v_air,gamma_air) P0_air*((v0_air/v_air)^gamma_air);  %Equation 5
%fpressure_end = @(P_end,m0_air,m_air,gamma_air) P_end*((m_air/m0_air)^gamma_air);
%fthrust = @(c_dis,a_throat,p_air,p_ambient) 2*c_dis*a_throat*(p_air-p_ambient);    %Equation 9
%fdrag = @(Cd,Ac,rho) Cd*Ac*0.5*rho;
%fvol_dt = @(c_dis,a_throat,rho_water,p_zero,p_ambient,v_zero,v_air,gamma_air) c_dis*a_throat*sqrt((2/rho_water)*((p_zero*(v_zero/v_air)^gamma_air)-p_ambient)); %oof - only for S1, at S2 vair = vb
%fmdot_water = @(c_dis,rho_water,a_throat,v_exhaust) c_dis*rho_water*a_throat*v_exhaust;
%magnitude = @(input) sqrt(input(1)^2+input(2)^2);

pressure_1 = c.i.P0_air*((c.i.v0_air/v_air)^c.gamma_air); %Thermodynamics
P_end = c.i.P0_air*(c.i.v0_air/c.v_empty)^(c.gamma_air); 
pressure_2 = P_end*((M_air/c.i.m0_air)^c.gamma_air); %Thermodynamics

%State 1 - Commencing the Booleans
if c.v_empty > v_air %((M_rkt-c.i.m0_air) > (c.m_empty))
     if z <= c.l_stand*sin(c.i.theta0)   %ON STAND??! State 0...Fixing the heading here...
        hbar = [cos(c.i.theta0); sin(c.i.theta0)];
       % disp 'On Stand'
     else
        hbar = hbar;
     end

exhaustvel = sqrt((2*(pressure_1-c.P_atm))/c.rho_water);

gravity = M_rkt*c.g;    %Kommencing Kinetics
thrust = 2*c.dis*A_throat*(pressure_1-c.P_atm);  %Scalar*heading unit vector
drag = c.i.Cd*A_rkt*0.5*c.rho_air*(mag_v)^2;   %Above  

mdot_water = -c.dis*c.rho_water*A_throat*exhaustvel; %5    ... From equation 9. This begins the other parameters, but is mdot(water)
vdot = c.dis*A_throat*sqrt((2/c.rho_water)*((c.i.P0_air*(c.i.v0_air/v_air)^c.gamma_air)-c.P_atm)); %6
mdot_air = 0;   %7

mass = M_rkt;

%disp 'Climbing' %, imaginary and velocity as follows'
%disp (num2str(pressure_1));

%State 2 - Thrust from exhausting air still
elseif (pressure_2 >= c.P_atm)   %No mo' water - 'b' handle
%P_end = pressure_1;
P_end = c.i.P0_air*(c.i.v0_air/c.v_empty)^(c.gamma_air);    %Constant value - At least it should be

pressure_2 = P_end*((M_air/c.i.m0_air)^c.gamma_air); %Thermodynamics 

rho_bottle = M_air/v_air;  %More thermodynamic definitions
temperature_bottle = pressure_2 / (rho_bottle*c.R_air);
pressure_crit = pressure_2 * (2 / (c.gamma_air + 1))^(c.gamma_air/(c.gamma_air-1)); %This is where the fun begins
vdot = 0; %6

if pressure_crit > c.P_atm  %Choked Flow
    mach = 1;
    exhausttemp = temperature_bottle*(2/(c.gamma_air+1));
    exhaustvel = sqrt(c.gamma_air*c.R_air*exhausttemp);
    exhaustpressure = pressure_crit;
    rho_exit = exhaustpressure / (c.R_air*exhausttemp);
    
  %  disp 'choked'
else %Non Choked
    mach = sqrt(((2*((pressure_2/c.P_atm)^((c.gamma_air-1)/c.gamma_air))-1)/(c.gamma_air-1)));
    exhausttemp = temperature_bottle / (1 + ((c.gamma_air - 1)/2)*mach^2);
    exhaustvel = mach * sqrt(c.gamma_air*c.R_air*exhausttemp);
    rho_exit = c.P_atm / (c.R_air * exhausttemp);
    exhaustpressure = c.P_atm;
  %  disp 'unchoked'
end

mdot_air = -c.dis*rho_exit*A_throat*exhaustvel;
mdot_water = 0; %5    ... From equation 9. This begins the other parameters, 
thrust = (-mdot_air*exhaustvel) + (exhaustpressure-c.P_atm)*A_throat;
drag = c.i.Cd*A_rkt*0.5*c.rho_air*(mag_v)^2;                        %0

gravity = M_rkt*c.g;    %Kommencing Kinetics

mass = M_rkt;
%disp 'On Air'

%Phase 3
elseif pressure_2 < c.P_atm && z > 0
    thrust = 0;
    gravity = c.g*M_rkt;

    drag = c.i.Cd*A_rkt*0.5*c.rho_air.*((mag_v).^2);
    mdot_water = 0;
    vdot = 0;
    mdot_air = 0;
    mass = M_rkt;
   % disp 'Water Exhausted'
%disp (num2str(hbar));
else
   %  disp 'On Ground'
     gravity = 0;
     thrust = 0;
     drag = 0;
     Vx = 0;
     Vz = 0;
     vdot = 0;
     mdot_water = 0;
     mdot_air = 0;
     mass = M_rkt;

end
f_net_x = thrust*(hbar(1)) - drag*(hbar(1));
f_net_z = thrust*(hbar(2)) - drag*(hbar(2)) - gravity;        %All set!                          %Drag is set to be inverse within states, and gravity assumed to go down
Ax = f_net_x/mass;
Az = f_net_z/mass;

dy_out = [Vx; Ax; Vz; Az; (mdot_water+mdot_air); vdot; mdot_air];
%X     = [X   Vx  Z   Vz     M_rkt    v_air    M_air]';

%Perhaps a simple convention could be good
telem = thrust;
imcheck = abs(imag(v_air));
if imcheck > 0 
    disp 'oof'
end
end

function c = getconst()
%Physical Constants
c.g = 9.81;
c.dis = 0.8;

%Thermodynamic Properties
c.R_air = 287;  %J/kg*k
c.gamma_air = 1.4;
c.rho_air = 0.961;  %kg/m^3
c.rho_water = 1000;

%Bottle Parameters
c.d_bottle = 0.105;
c.d_throat = 0.021; %m    - Need area from these
c.v_empty = 0.002;  %m^3
c.m_empty = 0.15; %kg

%Environmental Parameters
c.T0 = 300; %K
c.P_atm = 12.1*6894.757; %12.1 psig - 1 psi = 6.89476*10^3 Pa
c.l_stand = 0.5;    %length of stand

%All variable inputs are placed inside main...
i.P0_air = 52*6894.757 + c.P_atm;     %Thermodynamics (Payload)
i.water_v0 = 0.001;
i.water_m0 = c.rho_water*i.water_v0;
i.v0_air = c.v_empty - i.water_v0;
i.m0_air = (i.P0_air*i.v0_air)/(c.R_air*c.T0);
i.theta0 = deg2rad(42);   %Launch Parameters
i.Cd = 0.5; 
c.i = i;

c = c;
end

%{
Remarks.

Rocketry, even within such a small scale, is surprisingly complex.
The fact that these governing equations are so closely related makes this a
very complicated problem - One for which numerical integration is maybe
best suited. 
Thinking of other factors of an actual flight - Even something as simple as 
winds at altitude, not to mention the fluctuations in exhaust pressure or
instabilities in heading - would make this exponentially more complex.

Perhaps this is the merit of flight test; Yet, in spite of the difficulty
of manual troubleshooting, the ardure of aligning vectors to physical
reality, and the shenanigans of ODE45, this simulation is a powerful thing.
And indeed, ODE45, despite its shenanigans, is a quality integrator.
%}
