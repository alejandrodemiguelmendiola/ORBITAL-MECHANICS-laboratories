%% PLANET EXPRESS
clc;
close all;
clear all;

%% MISSION INFO 
planetid_dep = 3; % Earth = #3
planetid_arr = 4; % Mars = #4
constraint_launcher = 11;

%% SET INTEGRATION OPTIONS
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');
steps = 100;

%% LAUNCH WINDOWS
day1_early = date2mjd2000([2003,4,1,12,0,0]); %1/4/2003 in days
day1_late = date2mjd2000([2003,8,1,12,0,0]); %1/8/2003 in days
day2_early = date2mjd2000([2003,9,1,12,0,0]);%1/9/2003 in days
day2_late = date2mjd2000([2004,3,1,12,0,0]); %1/3/2004 in days

%% EPHEMERIDES (planet position at a specific date)
[kep1_early, mu_sun] = uplanet(day1_early,planetid_dep);
[kep1_late, ~] = uplanet(day1_late,planetid_dep);

[kep2_early, ~] = uplanet(day2_early,planetid_arr);
[kep2_late, ~] = uplanet(day2_late,planetid_arr);

%% FROM KEPLERIAN ELEMENTS TO CARTESIAN COORDINATES
[r_orbit1_early,v_orbit1_early] = kep2car(kep1_early(1),kep1_early(2),kep1_early(3),kep1_early(4),kep1_early(5),kep1_early(6),mu_sun);
[r_orbit1_late,v_orbit1_late] = kep2car(kep1_late(1),kep1_late(2),kep1_late(3),kep1_late(4),kep1_late(5),kep1_late(6),mu_sun);

[r_orbit2_early, v_orbit2_early] = kep2car(kep2_early(1),kep2_early(2),kep2_early(3),kep2_early(4),kep2_early(5),kep2_early(6),mu_sun);
[r_orbit2_late, v_orbit2_late] = kep2car(kep2_late(1),kep2_late(2),kep2_late(3),kep2_late(4),kep2_late(5),kep2_late(6),mu_sun);

%% PLANET 1 (EARTH) ORBITS
T1 = 2*pi*sqrt(kep1_early(1)^3/mu_sun);

tspan1 = linspace(0,T1,steps);
[~,dy_1_early] = ode113(@ode_2bodyproblem,tspan1,[r_orbit1_early,v_orbit1_early],options,mu_sun);
[~,dy_1_late] = ode113(@ode_2bodyproblem,tspan1,[r_orbit1_late,v_orbit1_late],options,mu_sun);

% PLANET 2 (MARS) ORBITS
T2 = 2*pi*sqrt(kep2_early(1)^3/mu_sun);

tspan2 = linspace(0,T2,steps);
[~,dy_2_early] = ode113(@ode_2bodyproblem,tspan2,[r_orbit2_early,v_orbit2_early],options,mu_sun);
[~,dy_2_late] = ode113(@ode_2bodyproblem,tspan2,[r_orbit2_late,v_orbit2_late],options,mu_sun);

%% TRANSFER ARC
steps_pork = 100;
tspan_dep = linspace(day1_early,day1_late,steps_pork);
tspan_arr = linspace(day2_early,day2_late,steps_pork);
deltav = zeros(steps_pork,steps_pork);
deltav_total = deltav; %to initialize the matrix

for i = 1:steps_pork
    for j = 1:steps_pork
        
        [kep_dep,mu_1] = uplanet(tspan_dep(j),planetid_dep);
        [kep_arr,mu_2] = uplanet(tspan_arr(i),planetid_arr);
        
        [RI,v_dep] = kep2car(kep_dep(1),kep_dep(2),kep_dep(3),kep_dep(4),kep_dep(5),kep_dep(6),mu_1);
        [RF,v_arr] = kep2car(kep_arr(1),kep_arr(2),kep_arr(3),kep_arr(4),kep_arr(5),kep_arr(6),mu_2);
        TOF = (tspan_arr(i) - tspan_dep(j)) * 24 * 3600; %OJO, in seconds
        [A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,mu_sun,0,0,0,0);
        deltav_1 = norm(VI - v_dep');
        deltav_2 = norm(v_arr' - VF);
        deltav(i,j) = deltav_1 + deltav_2; 
        
        deltav_total(i,j) = deltav(i,j); %to save all the results before imposing the launcher constraint
        if deltav_1 > constraint_launcher %NOTE: the launcher constraint is for deltav_1 or deltav
            deltav(i,j) = NaN;
        end
    end
end

%% PORK-CHOP PLOTS
[X,Y] = meshgrid(tspan_dep,tspan_arr);
drawpork = 'no';
switch drawpork
    case 'yes'
        figure
        contour(X,Y,deltav,'ShowText','on');
    case 'no'
        
end

%% MINIMUM DELTAV WITHOUT THE LAUNCHER CONSTRAINT
totalmin = min(deltav_total,[],'all');

%% MINIMUM DELTAV WITH THE LAUNCHER CONSTRAINT
[mindeltav,I] = min(deltav(:)); %transforms the matrix into a single row matrix and finds the minimum element and its location
[dep_index,arr_index] = ind2sub(size(deltav),I); %using the previous location(index), it finds the row and the column the minimum element belongs for the original matrix

%% DATES FOR THE LAUNCH AND ARRIVAL
t_dep = tspan_dep(dep_index); %mjd2000
t_arr = tspan_arr(arr_index); %mdj2000
departure_date = mjd20002date(t_dep); %19/6/2003 08:50
arrival_date = mjd20002date(t_arr); %11/12/2003 14:39

%% CALCULATE MINIMUM TRANSFER ARC
days_tof = t_arr - t_dep; %OJO, in days
tof = days_tof * 24 * 3600;
t_arc = linspace(0,tof,steps);

[kep_dep,~] = uplanet(t_dep,planetid_dep);
[kep_arr,~] = uplanet(t_arr,planetid_arr);

[r_dep,v_dep] = kep2car(kep_dep(1),kep_dep(2),kep_dep(3),kep_dep(4),kep_dep(5),kep_dep(6),mu_sun);
[r_arr,v_arr] = kep2car(kep_arr(1),kep_arr(2),kep_arr(3),kep_arr(4),kep_arr(5),kep_arr(6),mu_sun);

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r_dep,r_arr,tof,mu_sun,0,0,0,0);

[~,dy_arc] = ode113(@ode_2bodyproblem, t_arc, [r_dep,VI'], options, mu_sun);


%% PLOTS
plotting = '3D';
switch plotting
    
    case '3D'
        figure
        plot3(0,0,0,'y*')
        hold on
        plot3(dy_1_early(:,1),dy_1_early(:,2),dy_1_early(:,3),'b');  
        plot3(dy_2_early(:,1),dy_2_early(:,2),dy_2_early(:,3),'r');
        plot3(r_dep(1),r_dep(2),r_dep(3),'b*');
        plot3(r_arr(1),r_arr(2),r_arr(3),'r*');
        plot3(dy_arc(:,1),dy_arc(:,2),dy_arc(:,3),'g');
    
        xlabel('x[km]');
        ylabel('y[km]');
        zlabel('z[km]');
        legend('Sun','Planet 1','Planet 2','Departure','Arrival','Transfer arc');
        
    case '2D'
        figure
        plot(0,0,'y*')
        hold on
        plot(dy_1_early(:,1),dy_1_early(:,2),'b');  
        plot(dy_2_early(:,1),dy_2_early(:,2),'r');
        plot(r_dep(1),r_dep(2),'b*');
        plot(r_arr(1),r_arr(2),'r*');
        plot(dy_arc(:,1),dy_arc(:,2),'g');
    
        xlabel('x[km]');
        ylabel('y[km]');
        zlabel('z[km]');
        legend('Sun','Planet 1','Planet 2','Departure','Arrival','Transfer arc');
        
    case 'no'
    
    %CASE 3D MOVING NOT WORKING YET    
    case '3D moving'
        % EPHEMERIDES AT THE TIME OF LAUNCH (for the '3D moving' plot only)
        timebeforeplotting = 50; %in days mjd2000
        t_plot = linspace(t_dep - timebeforeplotting,t_arr,steps);
        [kep1_dep,~] = uplanet(t_dep - timebeforeplotting,planetid_dep);
        [kep2_dep,~] = uplanet(t_dep - timebeforeplotting,planetid_arr);
        [r1_dep,v1_dep] = kep2car(kep1_dep(1),kep1_dep(2),kep1_dep(3),kep1_dep(4),kep1_dep(5),kep1_dep(6),mu_sun);
        [r2_dep,v2_dep] = kep2car(kep2_dep(1),kep2_dep(2),kep2_dep(3),kep2_dep(4),kep2_dep(5),kep2_dep(6),mu_sun);
        [~,dy1_dep] = ode113(@ode_2bodyproblem,t_plot,[r1_dep,v1_dep],options,mu_sun);
        [~,dy2_dep] = ode113(@ode_2bodyproblem,t_plot,[r2_dep,v2_dep],options,mu_sun);
        
        
        % PLOTTING CODE
        h1_1 = animatedline('color','b');
        h1_2 = animatedline('Color','b','Marker','*');
        h2_1 = animatedline('Color','r');
        h2_2 = animatedline('Color','r','Marker','*');
        h3 = animatedline('Color','g');
        hold on
        for i = 1:steps
            plot3(0,0,0,'y*');
            plot3(dy1_dep(1,1),dy1_dep(1,2),dy1_dep(1,3),'bo');
            addpoints(h1_1,dy1_dep(i,1),dy1_dep(i,2),dy1_dep(i,3));
            addpoints(h1_2,dy1_dep(i,1),dy1_dep(i,2),dy1_dep(i,3));
            addpoints(h2_1,dy2_dep(i,1),dy2_dep(i,2),dy2_dep(i,3));
            addpoints(h2_2,dy2_dep(i,1),dy2_dep(i,2),dy2_dep(i,3));
            addpoints(h3,dy_arc(i,1),dy_arc(i,2),dy_arc(i,3));
            xlabel('x[km]');
            ylabel('y[km]');
            zlabel('z[km]');
            legend('Sun','Planet 1','Planet 2','Departure','Arrival','Transfer arc');
            drawnow limitrate
            clearpoints(h1_2)
            clearpoints(h2_2)
        end
       
end



