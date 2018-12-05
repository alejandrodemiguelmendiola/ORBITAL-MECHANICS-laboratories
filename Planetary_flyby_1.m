%% LAB5 ON FLYBY
clear all;
close all;
clc;

%% DATA 
mu_earth = 398600; %km^3/s^2
mu_sun = 132712*10^6; %km^3/s^2
AU = 149597870.7; %km
vminus_inf = [15.1,0,0]; %km/s
impact_par = 9200; %km
r_earth = [1,0,0] * AU; %km

%% SET INTEGRATION OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');
steps = 10000;

%% 5.1 SOLVING THE 2D HYPERBOLA (turn_angle,r_p and Deltav)
a = - mu_earth/norm(vminus_inf)^2; %km
turn_angle = 2* atan(-a/impact_par); %rad
e = 1./sind(turn_angle/2);
r_p = a.*(1-e); %km
V_earth = [0,sqrt(mu_sun/norm(r_earth)),0]; %pressumed circular

%% 5.2 Compute vplus_inf for three locations of the incoming asymptote: 
caseofstudy = 'front'; % cases: 'front', 'behind', 'under'

switch caseofstudy
    case 'front'
        rotation = [0,0,-1];
    case 'behind'
        rotation = [0,0,1];
    case 'under'
        rotation = [0,-1,0];      
end

vplus_inf = rotvector(vminus_inf,rotation,turn_angle);

%% 5.3 COMPUTE V- and V+
Deltav = vplus_inf - vminus_inf;
Vminus = V_earth + vminus_inf; %velocity of entry s/c wrt Sun
Vplus = Vminus + Deltav; %exit velocity s/c wrt Sun

%% 5.4 PLOT TRAJECTORIES
T = getT(r_earth,Vminus,mu_sun);
t = linspace(0,T/3,steps);
[~,dy1] = ode113(@ode_2bodyproblem, t,[r_earth,-Vminus],options,mu_sun); %the - on Vminus is just so it plots "backwards" as in the asignment
[~,dy2] = ode113(@ode_2bodyproblem, t,[r_earth,Vplus],options,mu_sun);

switch caseofstudy
    case {'front','behind'}
        figure
        plot(0,0,'y*','markersize',30);
        hold on
        plot(r_earth(1),r_earth(2),'b.','markersize',30);
        plot(dy1(:,1),dy1(:,2),'b');
        plot(dy2(:,1),dy2(:,2),'r');
        xlabel('x[km]');
        ylabel('y[km]');
        legend('Sun','Earth','Before flyby', 'After flyby');
        grid on
        
    case 'under'
        figure
        plot3(0,0,0,'y*','markersize',30);
        hold on
        plot3(r_earth(1),r_earth(2),r_earth(3),'b.','markersize',30);
        plot3(dy1(:,1),dy1(:,2),dy1(:,3),'b');
        plot3(dy2(:,1),dy2(:,2),dy2(:,3),'r');
        xlabel('x[km]');
        ylabel('y[km]');
        zlabel('z[km]');
        legend('Sun','Earth','Before flyby', 'After flyby');
        grid on
end



