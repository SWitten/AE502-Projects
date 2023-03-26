%   J2 ORBIT PERTURBATION EFFECT
%
%   This script generates plots showing the effect of J2 perturbations on
%   the elements of a given orbit.
%
%   Adapted from 'Orbital Mechanics for Engineering Students', Curtis, 2014
%   Example 12.6
%
%   Author: Sarah Wittenauer
%   Date: March 2023

close all
clc
clear

% ----------INPUTS----------

% Conversion Factors and Constants
deg = pi/180;       % Conversion factor from degrees to radians
RE = 6370;          % Radius of the earth [km]   
J2 = 0.00108;       % Earth's J2
mu = 3.986e5;       % Gravitational Parameter [km^3/s^2]

% Given Initial Orbital Elements
a0 = 26600;         % Initial Semimajor Axis [km]
i0 = 1.10654;       % Initial Inclination [rad]
e0 = 0.74;          % Initial Eccentricity [-]
w0 = 5*deg;         % Initial Argument of Perigee [rad]
RA0 = 90*deg;       % Initial Right Ascension of the Node [rad]
M0 = 10*deg;        % Initial Mean Anomaly Angle [rad]

% Calculated Orbit Parameters
rp0 = a0*(1-e0);                % Initial Perigee Radius [km]
h0 = sqrt(rp0*mu*(1+e0));       % Initial Angular Momentum [km^2/s]
T0 = 2*pi/sqrt(mu)*a0^(3/2);    % Initial Period [sec]
TA0 = true_anomaly(e0,M0);      % Intial True Anomaly [rad]

% Store the initial orbital elements in a vector
coe0 = [h0 e0 RA0 i0 w0 TA0];

%% ----------CALCULATIONS----------

% Set up the time vector and prepare for integration
t0 = 0;                         % Initial time [sec]
tf = 100*24*3600;               % Final Time [sec]
nout = 50000;                   % Number of data points [-]
tspan = linspace(t0,tf,nout);   % Time Vector
options = odeset('reltol',1e-8,'abstol',1e-8,'initialstep',T0/1000);
y0 = coe0;

% Run Integration on the rates subfunction using the parameters set up
% above
[t,y] = ode45(@rates,tspan,y0,options);

% Store solution in a variable for each parameter
h = y(:,1);
e = y(:,2);
RA = y(:,3);
i = y(:,4);
w = y(:,5);
TA = y(:,6);

%% ----------PLOT----------
% Plot the change in orbital elements (excluding the true anomaly) over the
% specified 100 days. Angles plotted as degrees and Time in days.

figure(1)
subplot(5,1,1)
plot(t/(24*3600),(RA - RA0)/deg)
title('Right Ascension (degrees)')
xlabel('days')
grid on
grid minor

subplot(5,1,2)
plot(t/(24*3600),(w - w0)/deg)
title('Argument of Perigee (degrees)')
xlabel('days')
grid on
grid minor

subplot(5,1,3)
plot(t/(24*3600),h - h0)
title('Angular Momentum (km^2/s)')
xlabel('days')
grid on
grid minor
ylim([-5,45])

subplot(5,1,4)
plot(t/(24*3600),e - e0)
title('Eccentricity')
xlabel('days')
grid on
grid minor
ylim([-5e-4,15e-4])

subplot(5,1,5)
plot(t/(24*3600),(i - i0)/deg)
title('Inclination (degrees)')
xlabel('days')
grid on
grid minor
ylim([-5e-3,20e-3])


%%
function dfdt = rates(t,f)

%RATES
%
%   Calculates d/dt of each orbital element
%   Adapted from 'Orbital Mechanics for Engineering Students', Curtis, 2014
%   Example 12.6 and Gauss' variational equations (Equations 12.89)
%
%   dfdt = rates(t,f)
%
%   INPUTS:     t, time vector
%               f, initial argument
%
%   OUTPUTS:    dfdt, d/dt of each element at each t

RE = 6370;          % Radius of the earth [km]   
J2 = 0.00108;       % Earth's J2
mu = 3.986e5;       % Gravitational Parameter [km^3/s^2]

% Store the Elements in Variables for ease of use
h = f(1);
e = f(2);
RA = f(3);
i = f(4);
w = f(5);
TA = f(6);

% Calculate the radius and argument of latitude
r = h^2/mu/(1+e*cos(TA));
u = w + TA;

% d/dt for each element (Equations 12.89)
hdot = -3/2*J2*mu*RE^2/r^3*sin(i)^2*sin(2*u);
edot = 3/2*J2*mu*RE^2/h/r^3 ...
*(h^2/mu/r*sin(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
-sin(2*u)*sin(i)^2*((2+e*cos(TA))*cos(TA)+e));
TAdot = h/r^2 + 3/2*J2*mu*RE^2/e/h/r^3 ...
*(h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
+ sin(2*u)*sin(i)^2*sin(TA)*(h^2/mu/r + 1));
RAdot = -3*J2*mu*RE^2/h/r^3*sin(u)^2*cos(i);
idot = -3/4*J2*mu*RE^2/h/r^3*sin(2*u)*sin(2*i);
wdot = 3/2*J2*mu*RE^2/e/h/r^3 ...
*(-h^2/mu/r*cos(TA)*(3*sin(i)^2*sin(u)^2 - 1) ...
- sin(2*u)*sin(i)^2*sin(TA)*(2 + e*cos(TA)) ...
+ 2*e*cos(i)^2*sin(u)^2);

% Store solution in single variable to pass back
dfdt = [hdot edot RAdot idot wdot TAdot]';

end