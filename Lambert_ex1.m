% Data
mu = 398600; %km^3/s^2;
r1 = [-21800, 37900,0];
r2 = [27300, 27700, 0];
deltat = 15*3600 + 6*60 + 40;

% Lambert solver, orbit reconstruction
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1, r2, deltat, mu,0,0,0);

% Set timespan
T = 2*pi*sqrt(A^3/mu);
step = 1000;
t = linspace(0,T,step);

% Initial conditions
y0 = [r1,VI];

% Set options
options = odeset('RelTol',1e-13, 'AbsTol',1e-14,'Stats','on');

% Orbit integration
[t,dy] = ode113(@ode_2bodyproblem, t, y0, options, mu);

% Plots
plot3(dy(:,1),dy(:,2),dy(:,3));
hold on
plot3(r1(1),r1(2),r1(3),'ro');
plot3(r2(1),r2(2),r2(3),'go');






