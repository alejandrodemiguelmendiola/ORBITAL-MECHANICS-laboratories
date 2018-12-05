clc;
clear;

%% DATA
mu = 398600; %km^3/s^2;

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

%% FROM KEPLERIAN TO CARTESIAN 
[r1,v1] = kep2car(a1,e1,i1,O1,o1,nu1, mu);

[r2,v2] = kep2car(a2,e2,i2,O2,o2,nu2, mu);

%% SET TIMESPAN
step = 1000;

T1 = 2*pi*sqrt(a1^3/mu);
t1 = linspace(0,T1,step);

T2 = 2*pi*sqrt(a2^3/mu);
t2 = linspace(0,T2,step);

tof = 3*T2;
t = [0,1,T1/8]; %chosen departure time is the last component of the vector, the rest is for making a vector for ode113 integration
%probably this can be done better

%% INITIAL CONDITIONS
y0_1 = [r1,v1];
y0_2 = [r2,v2];

%% SET OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');

%% ORBIT INTEGRATION
[~,dy1] = ode113(@ode_2bodyproblem,t1,y0_1,options,mu); %just to draw orbit 1
[~,dy2] = ode113(@ode_2bodyproblem,t2,y0_2,options,mu); %just to draw orbit 2

[~,aux1] = ode113(@ode_2bodyproblem,t,y0_1,options,mu); %to integrate the position in orbit 1 at any time t(end)
[~,aux2] = ode113(@ode_2bodyproblem,(t+tof),y0_2,options,mu); %to integrate the position in orbit 2 at any time t(end)

%% LAMBERT PROBLEM
P1 = [aux1(end,1),aux1(end,2),aux1(end,3)];
P2 = [aux2(end,1),aux2(end,2),aux2(end,3)];

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(P1,P2,tof,mu,0,0,0);

%% TRANSFER ARC INTEGRATION
t_arc = linspace(0,tof,step);
y0_3 = [P1,VI];
[~,dy3] = ode113(@ode_2bodyproblem,t_arc,y0_3,options,mu);

%% TRANSFER COSTS
v1 = [aux1(end,4),aux1(end,5),aux1(end,6)];
v2 = [aux2(end,4),aux2(end,5),aux2(end,6)];
deltav1 = norm(VI - v1);
deltav2 = norm(v2 - VF);
deltav = deltav1 + deltav2;

%% PLOT ORBITS
plottype = '3D'; %change between 'noplot', '2D' and '3D' to plot in different dimensions.
switch plottype
    case '3D'
        figure
        plot3(dy1(:,1),dy1(:,2),dy1(:,3),'b'); %orbit1
        hold on
        plot3(dy2(:,1),dy2(:,2),dy2(:,3),'g'); %orbit2
        plot3(r1(1),r1(2),r1(3),'b*');
        plot3(r2(1),r2(2),r2(3),'g*');
        plot3(dy3(:,1),dy3(:,2),dy3(:,3),'r'); %transfer arc
        plot3(0,0,0,'ko'); %Earth
        plot(P1(1),P1(2),P1(3),'bo'); %P1
        plot(P2(1),P2(2),P2(3),'go'); %P2
        legend('initial orbit','final orbit','initial dep point','initial fin point','transfer arc','Earth','P1','P2');
        xlabel('x[km]');
        ylabel('y[km]');
        zlabel('z[km]');
        
    case '2D'
        figure
        plot(dy1(:,1),dy1(:,2),'b'); %orbit1
        hold on
        plot(dy2(:,1),dy2(:,2),'g'); %orbit2
        plot(r1(1),r1(2),'b*');
        plot(r2(1),r2(2),'g*');
        plot(dy3(:,1),dy3(:,2),'r'); %transfer arc
        plot(0,0,'ko'); %earth
        plot(P1(1),P1(2),'bo'); %P1
        plot(P2(1),P2(2),'go'); %P2
        legend('initial orbit','final orbit','initial dep point','initial fin point','transfer arc','Earth','P1','P2');
        xlabel('x[km]');
        ylabel('y[km]');
        
    case 'noplot'
end

%% PORK-CHOP PLOTS

step_2 = 100;
vart = linspace(0.5*3600,5*3600,step_2); %different departure times
vartof = linspace(0,2.5*3600,step_2); %different variable times of flight

[~,aux3] = ode113(@ode_2bodyproblem,vart,y0_1,options,mu);
[~,aux4] = ode113(@ode_2bodyproblem,(vart+vartof),y0_2,options,mu);

deltav_2 = zeros(step_2,step_2);

for i = 1:step_2
    for j = 1:step_2
        
        P1_2 = [aux3(i,1),aux3(i,2),aux3(i,3)];
        v1_2 = [aux3(i,4),aux3(i,5),aux3(i,6)];
        P2_2 = [aux4(i,1),aux4(i,2),aux4(i,3)];
        v2_2 = [aux4(i,4),aux4(i,5),aux4(i,6)];
        [A,P,E,ERROR,VI_2,VF_2,TPAR,THETA] = lambertMR(P1_2,P2_2,vartof(j),mu,0,0,0);     

        deltav1_2 = norm(VI_2 - v1_2');
        deltav2_2 = norm(v2_2' - VF_2);
        deltav_2(i,j) = deltav1_2 + deltav2_2;
        end
end

[X,Y] = meshgrid(vart,vartof);
test = rand(step_2);

draw = 'n';
switch draw
    case 'y'
figure
contour(X,Y,deltav_2);
xlabel('t [s]');
ylabel('tof [s]');
legend('\Delta v');
figure
surf(X,Y,deltav_2);
    case 'n'
end
