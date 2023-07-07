clear; clc; close all;
b = 0.1;m = 1;g = 9.8;l = 1;
k1 = 0.2006;k2 = 0.4985;
A = [0,1; 0, -b/(m*l^2)];
%B = [0; (u - m*g*l*sin(z1))/(m*l^2)];
C = [1,0];
D = [C,0];
K = [k1;k2];
F = [A,zeros(2,2);zeros(2,2),(A - K*C)];
