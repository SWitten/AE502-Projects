%   ORBITAL ELEMENT CALCULATOR
%
%   This script generates the orbital elements of a body given it's state
%   vectors.
%
%   Author: Sarah Wittenauer
%   Date: February 2023

function coe = coe_from_sv(r_v,v_v,mu)

%COE_FROM_SV
%
%   Calculates Orbital Elements from the state vectors
%
%   coe =coe_from_sv(r_v,v_v,mu)
%
%   INPUTS:     r_v, position vector
%               v_v, velocity vector
%               mu, standard gravitational parameter
%
%   OUTPUTS:    coe, vector containing orbital elements:
%                   (1): h, magnitude of the angular momentum (km^2/sec)
%                   (2): e, eccentricity
%                   (3): RA, right ascension of the ascending node (rad)
%                   (4): incl, inclination (rad)
%                   (5): w, argument of perigee (rad)
%                   (6): TA, true anomaly (rad)
%                   (7): a, semimajor axis (km/nmi)

eps = 1.e-10;

r = norm(r_v);
v = norm(v_v);
vr = dot(r_v,v_v)/r;

h_v = cross(r_v,v_v);
h = norm(h_v);

incl = acos(h_v(3)/h);

n_v = cross([0 0 1],h_v);
n = norm(n_v);
if n ~=0
    RA = acos(n_v(1)/n);
    if n_v(2) < 0
        RA = 2*pi - RA;
    end
else
    RA = 0;
end

e_v = 1/mu*((v^2 - mu/r)*r_v - r*vr*v_v);
e = norm (e_v);

if n~= 0
    if e > eps
        w = acos(dot(n_v,e_v)/n/e);
        if e_v(3) < 0
            w = 2*pi - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end

if e>eps
    TA = acos(dot(e_v,r_v)/e/r);
    if vr < 0
        TA = 2*pi - TA;
    end
else
    cp = cross(n_v,r_v);
    if cp(3) >= 0
        TA = acos(dot(n_v,r_v)/n/r);
    end
end
a = h^2/mu/(1-e^2);
coe = [h e RA incl w TA a];
end
