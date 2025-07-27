%Main file to run shootin
%Parameter choosing
epsilon = 0.1; 
p = 0.2;
%b = 2 * b;
%b = 3 * b;
%b = 0.05;
%b = 0.1;
b = 0.18;
%s = 0.08;
s = 0.16;
%s = 2 * s;
% = 2 * b;
tfinal = 100;
            %heteroclinicOrbit(initalCond,onlyS,plotFullSystem,fullSystem,fastTimeScale,plotSlowFlow,fullSlowFlow,vecFieldSlowFlow, epsilon, b, p, s)
[t, X, Y] = heteroclinicOrbit("nontrivial",false,true,false,false,true,true,false,epsilon, b, p, s, tfinal);


%FOR TESTING WITH FULL RESCALLING 
%p = 2 * epsilon;
%s = 2 * epsilon;