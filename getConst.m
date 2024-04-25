function const = getConst()
const.g = 9.81; % m/s^2
const.c_dis = 0.8; % coeff
const.row_air = 0.961; % kg/m^3
const.Vol_bottle = .002; % meters^3
const.p_amb =  12.1 * 6894.76; % pascals
const.gam = 1.4; % specific heat ratio;
const.row_w = 1000; %kg/m^3
const.dia_throat = 0.021; % bottle throat diameter
const.dia_bottle = .105; % bottle diameter
const.R_air = 287; % J/(kg*K) specific gas constant of air
const.m_bottle = 0.15; % kg
const.C_D = 0.5; % ratio
const.p_r_i = 52 * 6894.76 + const.p_amb; % initial guage pressure in bottle
const.Vol_w_i = 0.001; % meters^3
const.T_i = 300; % K temperature of air
const.v_i = 0; % initial velocity;
const.theta_i = deg2rad(42); % radians initial rocket angle
const.x_i = 0; % initial horizontal distance
const.z_i = 0.25; % meters initial vertical distance
const.l_s = 0.5; % length of test stand.

end