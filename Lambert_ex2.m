clc;
clear;

% DATA
mu = 398600; %km^3/s^2;
tof = 3300; %s

a1 = 12000; %km
e1 = 0;
i1 = 0; %rad
O1 = 0; %rad
o1 = 0; %rad
nu1 = 120*pi/180; %rad

a2 = 9500; %km
e2 = 0.3;
i2 = 0; %rad
O2 = 0; %rad
o2 = 0; %rad
nu2 = 250*pi/180; %rad

% FROM KEPLERIAN TO CARTESIAN 
[r1,v1] = kep2car(a1,e1,i1,O1,o1,nu1, mu);
y1 = [r1,v1];

[r2,v2] = kep2car(a2,e2,i2,O2,o2,nu2, mu);
y2 = [r2,v2];

% SET TIMESPAN
T1 = 2*pi*sqrt(a1^3/mu);
step = 1000;
t1 = linspace(0,T1,step);

T2 = 2*pi*sqrt(a2^3/mu);
t2 = linspace(0,T2,step);

% INITIAL CONDITIONS
y0_1 = [r1,v1];
y0_2 = [r2,v2];

% SET OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');

% ORBIT INTEGRATION
[~,dy1] = ode113(@ode_2bodyproblem,t1,y0_1,options,mu);
[~,dy2] = ode113(@ode_2bodyproblem,t2,y0_2,options,mu);

% LAMBERT PROBLEM
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1,r2,tof,mu,0,0,0);

% TRANSFER COST
deltav1 = norm(VI - v1);
deltav2 = norm(v2 - VF);
deltav = deltav1 + deltav2; %is the result right?

% TRANSFER ARC
y0_3 = [r1,VI'];
t3 = linspace(0,tof,step);
[~,dy3] = ode113(@ode_2bodyproblem,t3,y0_3,options,mu);

plot3(dy1(:,1),dy1(:,2),dy1(:,3),'b'); %orbit1
hold on
plot3(dy2(:,1),dy2(:,2),dy2(:,3),'g'); %orbit2
plot3(r1(1),r1(2),r1(3),'bo'); %point1
plot3(r2(1),r2(2),r2(3),'go'); %point2
plot3(dy3(:,1),dy3(:,2),dy3(:,3),'r'); %transfer arc
plot3(0,0,0,'ko');
xlabel('[km]');
ylabel('[km]');
zlabel('[km]');
legend('orbit 1','orbit 2','initial point', 'final point', 'transfer arc','Earth');




