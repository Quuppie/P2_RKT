function [dstate_dt] = OdeFun(t,state,const,m_air_i)

%% Common Notation
% x: X position
% z: Z position
% lowercase v: velocity
% lowercase m: mass
% Vol: volume
% row: density
% uppercase A: Area;
% lowercase a: accel;
% uppercase R: Universal gas Constant
%% Second Part designates its Use
% For example, A_throat is the area of the throat or v_x is velocity in x
% direction.

%% Taking in state values
x = state(1);
v_x = state(2);
z = state(3);
v_z = state(4);
m_r = state(5);
Vol_air = state(6);
m_air = state(7);

%% determing Heading vector

v_air = [v_x, v_z];
airspeed = norm(v_air);

if norm([x-const.x_i, z-const.z_i]) > 0.5
    h = v_air/airspeed;
else
    h = [cos(const.theta_i), sin(const.theta_i)];
end


%% Pressure for ifelse test cases
p_end = const.p_r_i * ((const.Vol_bottle-const.Vol_w_i)/const.Vol_bottle)^const.gam; % I might be wrong
p = p_end * (m_air/m_air_i)^const.gam;

%% Phase 1 Water Remaining
if Vol_air < const.Vol_bottle
    
    % Definition of pressure in first phase
    p = const.p_r_i * ((const.Vol_bottle-const.Vol_w_i)/Vol_air)^const.gam; % I might be wrong
    % General Calculations
    A_throat = pi * (const.dia_throat/2)^2;
    v_exhaust = sqrt(((2 * (p - const.p_amb))/const.row_w));
    mass_flow_w = const.c_dis * const.row_w * A_throat * v_exhaust;

    %% Important Declarations
    % Forces
    F_thrust = h * abs(2 * const.c_dis * A_throat * (p - const.p_amb));
    % State Vector
    dstate_dt(5,1) = -mass_flow_w;
    dstate_dt(6,1) = const.c_dis*A_throat*v_exhaust;
    dstate_dt(7,1) = 0;

    %% Phase 2 No Water, Remaining pressure differential
elseif (p - const.p_amb) > 0.001

    % Common calculations not dependent on choked flow
    row = m_air/const.Vol_bottle;
    p_crit = p*(2/(const.gam + 1))^(const.gam/(const.gam-1));
    p_end = const.p_r_i * ((const.Vol_bottle-const.Vol_w_i)/const.Vol_bottle)^const.gam;
    T = p/(row*const.R_air);

    % Choked Flow
    if p_crit > const.p_amb

        T_e = T*(2/(const.gam + 1));
        v_e = sqrt(const.gam*const.R_air*T_e);
        p_e = p_crit;
        row_e = p_crit/(const.R_air*T_e);

    % Unchoked Flow
    else

        M_e = sqrt((2*(-1+(p/const.p_amb)^((const.gam-1)/const.gam)))/(const.gam-1));
        T_e = T/(1 + M_e^(2) * (const.gam-1)/2);
        row_e = const.p_amb / (const.R_air * T_e);
        p_e = const.p_amb;
        v_e = M_e*sqrt(const.gam*const.R_air*T_e);

    end

    A_throat = pi * (const.dia_throat/2)^2;
    dot_m_air = const.c_dis*row_e*A_throat*v_e;

    %% Important Declarations
    % Force
    F_thrust = h * abs(dot_m_air * v_e + (p_e - const.p_amb) * A_throat);
    % State Vector
    dstate_dt(5,1) = -dot_m_air;
    dstate_dt(6,1) = 0;
    dstate_dt(7,1) = -dot_m_air;

    %% Phase 3 / Ballistic Trajectory
else

    F_thrust = h * 0; %% water and air pressure is exhausted so thrust = 0.
    dstate_dt(5,1) = 0; %% mass and air volume are unchanging.
    dstate_dt(6,1) = 0;
    dstate_dt(7,1) = 0;
    
end
%% Compiling the forces into a 2d vector
Force_drag = - h * 0.5 * const.row_air * airspeed^2 *const.C_D * pi * (const.dia_bottle/2)^2;
F_grav = - (m_r * const.g)*[0,1];
F_net = F_thrust + F_grav + Force_drag;

%% If else statement to determine if and when the rocket hits the ground.
if z > 0
    dstate_dt(3,1) = v_z; % change in z position is z velocity
    dstate_dt(1,1) = v_x; % change in in x position is x velocity
    dstate_dt(2,1) = F_net(1)/m_r; % change in velocity or acceleration is F/ma
    dstate_dt(4,1) = F_net(2)/m_r; % change in velocity or acceleration is F/ma
else
    dstate_dt(3,1) = 0;
    dstate_dt(1,1) = 0;
    dstate_dt(2,1) = 0;
    dstate_dt(4,1) = 0;
end


end
