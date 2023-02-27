%   LAMBERT SOLVER
%
%   This script solves lambert's problem from initial and final position
%   vectors and a given time-of-flight
%   
%   Script adapted from Algorithim 5.2 in "Orbital Mechanics for
%   Engineering Students", H.D.Curtis, 2013.
%
%   Author: Sarah Wittenauer
%   Date: February 2023

function [v1_v,v2_v] = lambert(r1_v,r2_v,t,string,mu)

%LAMBERT
%
%   Solves Lambert's Problem
%
%   [v1_v,v2_v]=lambert(r1_v,r2_v,t,string,mu)
%
%   INPUTS:     r1_v, initial position vector of the transfer orbit
%               r2_v, final position vector of the transfer orbit
%               t, given time-of-flight
%               string, 'pro' for finding the prograde orbit
%                       'retro' for finding the retrograde orbit
%               mu, standard gravitational parameter
%
%   OUTPUTS:    v1_v, initial velocity vector of the transfer orbit
%               v2_v, final velocity vector of the transfer orbit

r1 = norm(r1_v);
r2 = norm(r2_v);
c12 = cross(r1_v,r2_v);
theta = acos(dot(r1_v,r2_v)/r1/r2);

if strcmp(string,'pro')
    if c12(3) <= 0 
        theta = 2*pi-theta;
    end
elseif strcmp(string,'retro')
    if c12(3) >= 0
        theta = 2*pi-theta;
    end
else
    if c12(3) <= 0 
        theta = 2*pi-theta;
    end
    fprintf('**Prograde trajectory assumed**')
end

A = sin(theta)*sqrt(r1*r2/(1-cos(theta)));
z = 0;
while F(z,t,mu,A,r1,r2) < 0
    z = z+0.1;
end

tol = 1.e-6;
nmax = 1000;

ratio = 1;
n = 0;
while (abs(ratio) > tol) && (n <= nmax)
    n = n+1;
    ratio = F(z,t,mu,A,r1,r2)/dFdz(z,A,r1,r2);
    z = z - ratio;
end

if n >= nmax
    fprintf('**The number of iterations exceeds the max**')
end

f = 1-y(z,r1,r2,A)/r1;
g = A*sqrt(y(z,r1,r2,A)/mu);
gdot = 1-y(z,r1,r2,A)/r2;
v1_v = 1/g*(r2_v-f*r1_v);
v2_v = 1/g*(gdot*r2_v-r1_v);

% SUBFUNCTIONS

    function dum = y(z,r1,r2,A)
        dum = r1+r2+A*(z*stumpff_s(z)-1)/sqrt(stumpff_c(z));
    end

    function dum = F(z,t,mu,A,r1,r2)
        dum = (y(z,r1,r2,A)/stumpff_c(z))^1.5*stumpff_s(z)+A*sqrt(y(z,r1,r2,A))-sqrt(mu)*t;
    end

    function dum = dFdz(z,A,r1,r2)
        if z == 0
            dum = sqrt(2)/40*y(0,r1,r2,A)^1.5 + A/8*(sqrt(y(0,r1,r2,A))+A*sqrt(1/2/y(0,r1,r2,A)));
        else
            dum = (y(z,r1,r2,A)/stumpff_c(z))^1.5*(1/2/z*(stumpff_c(z)-3*stumpff_s(z))...
                +3*stumpff_s(z)^2/4/stumpff_c(z))+A/8*(3*stumpff_s(z)/stumpff_c(z)*sqrt(y(z,r1,r2,A))...
                +A*sqrt(stumpff_c(z)/y(z,r1,r2,A)));
        end
    end

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
