%   ORBIT OF A SATELLITE 
%
%   This script generates a plot showing the orbit of a satellite in
%   Cartesian space from it's Hamiltonian expressed in Delaunay variables
%   (equations from HW3)
%
%   Author: Sarah Wittenauer
%   Date: April 2023

clc
clear

% INPUTS
a = 1;
e = 0.5;
i = 45*pi/180;
t = linspace(0,99,1000);

% CALCULATIONS
T = 2*pi*a^(3/2);
M = (2*pi/T)*t;
E = true_anomaly(e,M);

x = a*(1-e*cos(E))*cos(i);
[x,y]= pol2cart(E,x);
z = a*(1-e*cos(E))*sin(i);

% PLOT
plot3(x,y,z)
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
grid on

