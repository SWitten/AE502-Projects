%   PORKCHOP PLOT GENERATOR
%
%   Author: Sarah Wittenauer
%   Date: February 2023

close all
clear
clc

%----------INPUTS----------

%Conversion Factors and Constants
mu = 1.327e11;      %Standard Gravitational Parameter for the Sun [km^3/s^2]
au_km = 1.496e8;    %Conversion from au to km
aud_kms = 1731.46;  %Conversion from au/day to km/s
d_s = 86400;        %Conversion from days to seconds

% Initial State vectors of earth and the planetary body
% Epoch is 01-01-2017 00:00:00 UTC 
re_v = [-1.796136509111975e-1, 9.667949206859814e-1, 3.668681017942158e-5]*au_km; %km
ve_v = [-1.720038360888334e-2, -3.211186197806460e-3, 7.927736735960840e-7]*aud_kms; %km/s
r1_v = [7.249472033259724, 14.61063037906177, 14.24274452216359]*au_km; %km
v1_v = [-8.241709369476881e-3, -1.156219024581502e-2, -1.317135977481448e-2]*aud_kms; %km/s

% Departure and Arrival Date Inputs
departure_start = juliandate(2017,1,1);
departure_end = juliandate(2020,7,31);
arrival_start = juliandate(2019,6,1);
arrival_end = juliandate(2022,1,31);
delta = arrival_start-departure_start;

%departure_delta = departure_end-departure_start;
step_dep = 2;
step_arr = 4;

% Create Arrays with the position and velocity vectors during departure and
% arrival
JD_dep = departure_start:step_dep:departure_end;
rArray_dep = zeros(length(JD_dep),3);
vArray_dep = zeros(length(JD_dep),3);

for i=1:length(JD_dep)
    dt=((step_dep*(i-1))*d_s);
    [rArray_dep(i,:),vArray_dep(i,:)]=position(mu,dt,re_v,ve_v);
end

JD_arr = arrival_start:step_arr:arrival_end;
rArray_arr = zeros(length(JD_arr),3);
vArray_arr = zeros(length(JD_arr),3);

for i=1:length(JD_arr)
    dt=(step_arr*(i-1)+delta)*d_s;
    [rArray_arr(i,:),vArray_arr(i,:)]=position(mu,dt,r1_v,v1_v);
end


%%
delta_dep = zeros(length(JD_dep),length(JD_arr));
delta_arr = zeros(length(JD_dep),length(JD_arr));
dv1_results = zeros(length(JD_dep),length(JD_arr));
dv2_results = zeros(length(JD_dep),length(JD_arr));

for i=1:length(JD_dep)
    JDi = JD_dep(i);
    for j=1:length(JD_arr)
        JDf = JD_arr(j);
        if JDi < JDf
            delta_dep(i,j) = JDi-departure_start;
            delta_arr(i,j) = JDf-arrival_start;
            TOF = (JDf-JDi)*d_s;
            [v1p,v2p] = lambert(rArray_dep(i,:),rArray_arr(j,:),TOF,'pro',mu);
            [v1r,v2r] = lambert(rArray_dep(i,:),rArray_arr(j,:),TOF,'retro',mu);
            dvd_p = norm(v1p-vArray_dep(i,:));
            dvd_r = norm(v1r-vArray_dep(i,:));
            dva_p = norm(v2p-vArray_arr(j,:));
            dva_r = norm(v2r-vArray_arr(j,:));
            dv1_results(i,j) = min(dvd_p,dvd_r);
            dv2_results(i,j) = min(dvd_p+dva_p,dvd_r+dva_r);
            clc
            disp(j)
            disp(i)
        else
            dv1_results(i,j) = nan;
            dv2_results(i,j) = nan;
        end
    end
end

%% PLOT

dv1_levels = [1 5 10 12.5 15 17.5 20];
dv2_levels = [37.5 40 42.5 45 47.5 50 52.5 55 57.5 60];

figure(1)
[c1,h1] = contourf(delta_dep,delta_arr,dv1_results,dv1_levels);
colormap jet
grid on
set(gca, 'layer', 'top');
colorbar;
xlabel('Departure (Days past 01-Jan-2017)')
ylabel('Arrival (Days past 01-Jun-2019)')
title('Earth to Borisov Fly-by Trajectories')

figure(2)
[c2,h2] = contourf(delta_dep,delta_arr,dv2_results,dv2_levels);
colormap jet
grid on
set(gca, 'layer', 'top');
colorbar;
xlabel('Departure (Days past 01-Jan-2017)')
ylabel('Arrival (Days past 01-Jun-2019)')
title('Earth to Borisov Rendezvous Trajectories')


