%Main file to run shooting
%Parameter choosing
epsilon = 0.1; 
b = 0.05;
%b = 0.15;
p = 0.2;
%s = 0.0868;
s = 0.2;

[t, X, Y] = heteroclinicOrbit(epsilon, b, p, s);