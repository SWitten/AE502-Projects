%   TRUE ANOMALY CALCULATOR
%
%   This script generates the true anomaly from the mean anomaly and
%   eccentricity of an orbit.
%
%   Author: Sarah Wittenauer
%   Date: March 2023

function f = true_anomaly(e,M)

%TRUE_ANOMALY
%
%   Calculates true anomaly from mean anomaly and eccentricity
%
%   f = true_anomaly(e,M)
%
%   INPUTS:     e, eccentricity
%               M, mean anomaly [rad]
%
%   OUTPUTS:    f, true anomaly [rad]

tol = 1e-8;

if M < pi
    E = M + e/2;
else
    E = M - e/2;
end

ratio = 1;
while abs(ratio) > tol
    ratio = (E - e*sin(E) - M)/(1-e*cos(E));
    E = E - ratio;
end

B = e/(1+sqrt(1-e^2));
f = E + 2*atan(B*sin(E)/(1-B*cos(E)));
end
