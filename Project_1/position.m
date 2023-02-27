%   TWO-BODY ORBIT PROPAGATOR
%
%   This script generates new state vectors based on initial state vectors
%   of an object and a given change in time
%
%   Author: Sarah Wittenauer
%   Date: February 2023

function [r_v,v_v] = position(mu,delta_t,r0_v,v0_v) 
%POSITION
%
%   Solves Time-of-Flight Problems
%   Valid for ellipse conics and hyperbola conics
%
%   [r_v,v_v]=position(mu,delta_t,r0_v,v0_v)
%
%   INPUTS:     mu, standard gravitational parameter
%               delta_t, change in time between initial and final state
%               vectors
%               r0_v, initial position vector
%               v0_v, initial velocity vector
%
%   OUTPUTS:    r_v, updated position vector
%               v_v, updated velocity vector

r0 = norm(r0_v);
v0 = norm(v0_v);
alpha = (2/r0)-(v0^2/mu);

x = universal_variable_f(mu,delta_t,r0_v,v0_v);

z = alpha*x^2;

f = 1-x^2/r0*stumpff_c(z);
g = delta_t - 1/sqrt(mu)*x^3*stumpff_s(z);
r_v = f*r0_v+g*v0_v;
r = norm(r_v);

fdot = sqrt(mu)/r/r0*(z*stumpff_s(z)-1)*x;
gdot = 1-x^2/r*stumpff_c(z);
v_v = fdot*r0_v+gdot*v0_v;
end

function c = stumpff_c(z)
%STUMPFF_C
%
%   Evaluates the Stumff series C(z) for a given z value
% 
%   INPUTS:     z, Input Argument
%   OUTPUTS:    c, value of C(z) at given input
 
    if z>0
        c = (1-cos(sqrt(z)))/z;
    elseif z<0
        c = (cosh(sqrt(-z))-1)/(-z);
    else
        c = 1/2;
    end
end

function s = stumpff_s(z)
%STUMPFF_S
%
%   Evaluates the Stumff series S(z) for a given z value
% 
%   INPUTS:     z, Input Argument
%   OUTPUTS:    s, value of S(z) at given input
 
    if z>0
        s = (sqrt(z)-sin(sqrt(z)))/(sqrt(z))^3;
    elseif z<0
        s = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z))^3;
    else
        s = 1/6;
    end
end

function x = universal_variable_f(mu,delta_t,r0_v,v0_v)
%UNIVERSAL_VARIABLE_F
%
%   Solves Universal Variable Equation using Newton's Method
%   Valid for ellipse conics and hyperbola conics
%
%   [x] = universal_variable_f(mu,delta_t,r0_v,v0_v)
%
%   INPUTS:     mu, standard gravitational parameter
%               delta_t, time to find universal variable
%               r0_v, position vector
%               v0_v, velocity vector
%
%   OUTPUTS:    x, universal variable, (km^0.5 or ft^0.5)

tol = 1.e-8;
r0 = norm(r0_v);
v0 = norm(v0_v);
alpha = (2/r0)-(v0^2/mu);
x = sqrt(mu)*abs(alpha)*delta_t;

n = 0;
ratio = 1;
while abs(ratio)> tol
    n = n+1;
    z = (x^2)*alpha;
    c_n = stumpff_c(z);
    s_n = stumpff_s(z);

    F = x^3*s_n + (dot(r0_v,v0_v)/sqrt(mu))*x^2*c_n+r0*x*(1-z*s_n)-delta_t*sqrt(mu);
    dFdx = x^2*c_n+(dot(r0_v,v0_v)/sqrt(mu))*x*(1-z*s_n)+r0*(1-z*c_n);
    ratio = F/dFdx;
    x = x-ratio;
end
end
