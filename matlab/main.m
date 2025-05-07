%Main file to run shooting
%Parameter choosing
epsilon = 0.1; 
%b = 0.08;
%b = 0.04;
%b = 0.1;
b = 0.04;
p = 0.2;
%s = 0.0868;
s = 0.2;

[t, X, Y] = heteroclinicOrbit("X1",true,true,false,false,false,epsilon, b, p, s);