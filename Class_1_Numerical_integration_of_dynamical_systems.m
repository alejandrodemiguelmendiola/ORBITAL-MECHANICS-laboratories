%juanluis.gonzalo@polimi.it%
%Numerical integration of dynamical systems%
%ode solvers for orbital mechanics--> ode113%
%SYNTAX: [t,y]= odeXX(odefun,tspan,y0,options)%
%OPTIONS:options=odeset('RelTol', 1e-13, 'AbsTol',1e-14)%

%Oscillator parameters%
%c = 1;
%m = 1;
%k = 1;
%gamma = c/(2*m);
%w0 = sqrt(k/m);
%w = sqrt(w0^2-gamma^2);
gamma = 1;
w0 = 2;
w = sqrt(w0^2-gamma^2);


%Initial conditions%
x0 = 1;
v0 = 1; %xdot%

%Timespan%
t = linspace(0,10,100);

%Analytic equations%
x = exp(-gamma.*t).*(x0*cos(w.*t)+(v0+gamma*x0)/w * sin(w.*t));

v = exp(-gamma.*t).* (v0*cos(w.*t)- (w0^2*x0+gamma*v0)/w * sin(w.*t));

%Set options%
options = odeset('RelTol',1e-13, 'AbsTol',1e-14)

%ODE solver%
odefun = @(t,y) [y(2); -2*gamma*y(2)-w0^2*y(1)]; 

[t,y] = ode113(odefun, t, [x0;v0],options);

%Results%
figure

subplot(2,1,1);
plot(time,x);
hold
plot(t,y(:,1),'rx');
xlabel('t');
ylabel('x');
legend('Analytic Position','Numeric Position');

subplot(2,1,2);
plot(time,v);
hold
plot(t,y(:,2),'rx');
xlabel('t');
ylabel('v');
legend('Analytic Velocity','Numeric Velocity');

