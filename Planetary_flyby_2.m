%% LAB5 ON POWERED FLYBY
clear all;
close all;
clc;

%% SET INTEGRATION OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');
steps = 1000;

%% DATA
muS= 132712*10^6; %km^3/s^2
muE = 398600; %km^3/s^2
AU = 149597870.7; %km
Vm = [31.5, 4.69, 0]; %km/s
Vp = [38.58, 0, 0]; %km/s
r_earth = [0,-1,0]*AU; %km
radius_earth = 6371; %km
h_atm_earth = 100; %km
% Compute the planet's linear speed
n_earth = sqrt(muS/norm(r_earth)^3); % Earth's rotation speed around the Sun
n_earth_v = n_earth*[0, 0, 1]; % % Earth's rotation speed vector [s^-1]
V_earth = cross(n_earth_v,r_earth); % Linear velocity of Earth

%% 1.Compute vp_inf and vm_inf
DeltaV = Vp - Vm;
vm_inf = Vm - V_earth;
vp_inf = DeltaV + vm_inf;

%% 2.Compute the turning angle delta
delta = acos(dot(vp_inf,vm_inf)/(norm(vm_inf)*norm(vp_inf)));

%% 3.Solve the non-linear system for rp

x0 = radius_earth;
opts = optimoptions('fsolve','OptimalityTolerance',1e-13);
rp = fsolve(@(x) root(x,norm(vm_inf),norm(vp_inf),delta,muE),x0,opts)
em = (1+rp*norm(vm_inf)^2/muE);
ep = (1+rp*norm(vp_inf)^2/muE);
deltam = 2*asin(1/em);
deltap = 2*asin(1/ep);

%% 3.Check rp validity
validity = rp > radius_earth + h_atm_earth

%% 4.Compute velocities at pericenter and Deltavp
vmp_inf = sqrt(norm(vm_inf)^2 + 2*muE/rp); %eq.8.58 from Curtis, v at perigee
vpp_inf = sqrt(norm(vp_inf)^2 + 2*muE/rp);
Deltavp = vpp_inf - vmp_inf
h_ga = rp - radius_earth

%% DEFINING THE FSOLVE FUNCTION
function F = root(x,vm_inf,vp_inf,delta,muP) %rp,em,ep,deltam,deltap

    rp = x(1);
    em = (1+rp*(vm_inf)^2/muP);
    ep = (1+rp*(vp_inf)^2/muP);
    deltam = 2*asin(1/em);
    deltap = 2*asin(1/ep);
    
    F(1) = delta - deltam/2-deltap/2;
end

