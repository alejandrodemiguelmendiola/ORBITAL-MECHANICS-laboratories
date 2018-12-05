%% LAB5 ON FLYBY
clear all;
close all;
clc;

%% DATA 
mu_earth = 398600; %km^3/s^2
mu_sun = 132712*10^6; %km^3/s^2
AU = 149597870.7; %km
vminus_inf = [15.1,0,0]; %km/s
r_earth = [1,0,0] * AU; %km
V_earth = [0,sqrt(mu_sun/norm(r_earth)),0]; %pressumed circular

%% SET INTEGRATION OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');
steps = 10000;

%% 1-Choose a location for the incoming asymptote
location = [1,0,0]* AU;

%% 2-Solve the hyperbola
impact_par = [5200,9200,15200];
a = - mu_earth/norm(vminus_inf)^2; %km
turn_angle = zeros(1,numel(impact_par));
e = zeros(1,numel(impact_par));
r_p = zeros(1,numel(impact_par));

for i = 1:numel(impact_par)
    
    turn_angle(i) = 2* atan(-a/impact_par(i)); %rad
    e = 1./sind(turn_angle/2);
    r_p = a.*(1-e); %km
    
end

%% 3-Compute V-
Vminus = V_earth + vminus_inf;

%% 4-Compute vplus_inf and V+

for i = 1:numel(impact_par)
    vplus_inf = rotvector(vminus_inf,[0,0,1],turn_angle(i));
    Deltav = vplus_inf - vminus_inf;
    Vplus = Vminus + Deltav;
    
    T = getT(r_earth,Vminus,mu_sun);
    t = linspace(0,T/3,steps);
    [~,dy1] = ode113(@ode_2bodyproblem, t,[r_earth,-Vminus],options,mu_sun); %the - on Vminus is just so it plots "backwards" as in the asignment
    [~,dy2] = ode113(@ode_2bodyproblem, t,[r_earth,Vplus],options,mu_sun);
    
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









