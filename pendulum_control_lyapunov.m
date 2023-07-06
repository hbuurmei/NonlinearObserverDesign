% Code for Validating a Control Lyaponov Function 
% for a Pendulum with an Observer

clear; clc; close all; format compact;
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(groot, 'defaultAxesTickLabelInterpreter','latex');% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'defaultFigureUnits', 'pixels', 'defaultFigurePosition', [440   378   560   380]);
set(0,'defaultAxesFontSize',24);
set(0, 'DefaultLineLineWidth', 2);
%syms b m l g real % damping, mass, length, gravity
b = 0.1;m = 1;g = 9.8;l = 1;
syms theta theta_dot real % theta = output
syms theta_hat theta_dot_hat real
syms u real % u = control input
%syms k1 k2 real % k1, k2 = observer gains
syms z1 z2 z3 z4 real
z = [z1; z2; z3; z4];
k1 = 0.2006;k2 = 0.4985;
A = [0,1; 0, -b/(m*l^2)];
B = [0; (u - m*g*l*sin(z1))/(m*l^2)];
C = [1,0];
D = [C,0];
%K = [k1; k2];
K = [0.1;0];

F = [A,zeros(2,2);zeros(2,2),(A - K*C)];
G = [B;0;0];
z0 = [pi;0;0;0];
x = [theta; theta_dot];
x_hat = [theta_hat; theta_dot_hat];
e = x_hat - x;
%z = [x;e];
%ii = rand(2);
ii = ones(2,2);
Q = ii*ii.';

%V = 0.5*m*l^2*(theta_dot)^2 + m*g*l*(1+cos(theta)) + 0.5*e'*Q*e;
V = 0.5*m*l^2*(z2)^2 + m*g*l*(1+cos(z1)) + 0.5*[z3; z4]'*Q*[z3; z4];
gradient(V,z)
Vdot = ((dot(gradient(V,z),F*z + G)));
% Lg = ((dot(gradient(V,z),G)));
% Vdot = Lf + Lg;

VdotFun = matlabFunction(Vdot);
VFun = matlabFunction(V);

z1 = linspace(-pi,pi,100); % theta = z1
z2 = linspace(-100,100,100); % theta_dot = z2
%z3 = linspace(-pi,pi,100); % theta_hat = z3
z3 = zeros(100,1);
z4 = zeros(100,1);

[Z1,Z2] = meshgrid(z1,z2);
Vdot_ = VdotFun(-100*Z2,Z1,Z2,z3,z4);
figure
surf(Z1,Z2,Vdot_)
xlabel('$\theta$');
ylabel('$\dot{\theta}$');
zlabel('$\dot{V}(z)$');
title('Zero Error');
%zlim([-Inf 0])
% figure
% surf(THETA,THETA_DOT,VFun(THETA,THETA,THETA_DOT,THETA_DOT));
M = max(Vdot_,[],"all")





