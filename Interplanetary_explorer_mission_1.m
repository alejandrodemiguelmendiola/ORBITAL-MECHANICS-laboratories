%% ASSIGNMENT 1.NEPTUNE DEPARTURE - URANUS FLYBY - JUPITER ARRIVAL
clc;
clear all;
close all;

%% DATA
%gravitational constant [km^3/s^2]
[muS] = astroConstants(4); 
[muN] = astroConstants(18);
[muU] = astroConstants(17);
[muJ] = astroConstants(15);
%mean radius [km]
[meanrN] = astroConstants(28);
[meanrU] = astroConstants(27);
[meanrJ] = astroConstants(25);
%mean velocity [km/s]
[meanvN] = sqrt(muN/meanrN); 
[meanvU] = sqrt(muU/meanrU);
[meanvJ] = sqrt(muJ/meanrJ);
%atmosphere height [km]
atm_U = 500000; 

%% INTEGRATION OPTIONS
options_int = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');
step_int = 1000; %steps for the integration

%% LAMBERT OPTIONS
orbitType = 0; %OJO!!! ASK IF WE SHOULD CONSIDER RETROGRADE TRANSFERS
Nrev = 0; % ALSO ASK
Ncase = 0;
optionsLMR = 1;

%% LAUNCH WINDOWS
first_dep = date2mjd2000([2020,1,1,0,0,0]); %1/1/2020 in julian days
last_dep = date2mjd2000([2022,1,1,0,0,0]); %1/1/2030 in julian days

%first_arr = date2mjd2000([2050,1,1,0,0,0]);
%last_arr = date2mjd2000([2070,1,1,0,0,0]);

step_dis = 200; %steps for the discretization
window = linspace(first_dep,last_dep,step_dis);
TOF = linspace(1e8,1e9,step_dis);
%TOFmax = first_arr - last_dep;
%TOF = linspace(3e8, TOFmax, step_dis);

%% CALCULATING ORBITAL ELEMENTS
counter = 0;
counter2 = 0;
ef_dp = zeros(step_dis,6); %efemerides (r,v) of the departure planet
for v = window
    counter = counter + 1; 
    [kep1,~] = uplanet(v,8);
    [r1,v1] = kep2car(kep1(1),kep1(2),kep1(3),kep1(4),kep1(5),kep1(6),muS);
    ef_dp(counter,1) = r1(1);
    ef_dp(counter,2) = r1(2);
    ef_dp(counter,3) = r1(3);
    ef_dp(counter,4) = v1(1);
    ef_dp(counter,5) = v1(2);
    ef_dp(counter,6) = v1(3);
end  
    

%% LAMBERT ARC NEPTUNE-URANUS
Deltav1 = zeros(step_dis);
for i = 1:step_dis
    RI = [ef_dp(i,1),ef_dp(i,2),ef_dp(i,3)];
    v1 = [ef_dp(i,4),ef_dp(i,5),ef_dp(i,6)];
    for j = 1:step_dis
        [kep2,~] = uplanet(window(i) + TOF(j)/(24*3600),7);
        [r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muS);
        RF = r2';
        [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF(j),muS,0,0,0);
        Deltav1(i,j) = norm(VI - v1) + norm(v2' - VF); 
    end
end

%find the minimum Deltav1
[M,I] = min(Deltav1(:));
M
[I_row, I_col] = ind2sub(size(Deltav1),I);
%Lambert 1 chosen
dep_chosen = window(I_row); %day chosen for departure
arrGA_chosen = dep_chosen + TOF(I_col)/(3600*24); %arrival to the GA planet chosen
[kep2,~] = uplanet(arrGA_chosen,7); %efemerides (r,v) of the GA planet
[r2,v2] = kep2car(kep2(1),kep2(2),kep2(3),kep2(4),kep2(5),kep2(6),muS);
RI = [ef_dp(I_row,1),ef_dp(I_row,2),ef_dp(I_row,3)];
RF = r2';
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF(I_col),muS,0,0,0);

%% POWERED GA
Vm = VF; %OJOOOO, I DONT KNOW IF THIS IS OKAY

%% ORBIT INTEGRATION (for plotting)
%initial conditions
r0_N = [ef_dp(I_row,1),ef_dp(I_row,2),ef_dp(I_row,3)];
v0_N = [ef_dp(I_row,4),ef_dp(I_row,5),ef_dp(I_row,6)];
r0_U = r2;
v0_U = v2;
r0_lambert = r0_N;
v0_lambert = VI;
%integration time
T_N = getT(r0_N,v0_N,muS);
T_U = getT(r0_U,v0_U,muS);
t_N = linspace(0,T_N,step_int);
t_U = linspace(0,T_U,step_int);
t_lambert = linspace(0,TOF(I_col),step_int);
%integration
[~,dy1] = ode113(@ode_2bodyproblem,t_N,[r0_N,v0_N],options_int,muS);
[~,dy2] = ode113(@ode_2bodyproblem,t_U,[r0_U,v0_U],options_int,muS);
[~,dylambert] = ode113(@ode_2bodyproblem,t_lambert,[r0_lambert,v0_lambert],options_int,muS);

%% PLOT
figure('Name','Mission')
plot3(0,0,0,'y*','markersize',8);
hold on
plot3(dy1(:,1),dy1(:,2),dy1(:,3),'b');
plot3(dy2(:,1),dy2(:,2),dy2(:,3),'r');
plot3(r0_N(1),r0_N(2),r0_N(3),'bo');
plot3(r0_U(1),r0_U(2),r0_U(3),'ro');
plot3(dylambert(:,1),dylambert(:,2),dylambert(:,3),'g');

xlabel('x[km]');
ylabel('y[km]');
zlabel('z[km]');
legend('Sun','Neptune orbit','Uranus orbit','Dep1','Arr1');

%% Degrees of freedom --> departure time, tof, flyby time
