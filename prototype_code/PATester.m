clc;

E = 1; % Steel thermal conductivity W/(m.K)
A = 1; % Cross-sectional area in m^2
n_el = 5;
h = ones(n_el, 1);

f = 1; % No internal heat generation
T0 = 5; % Temperture at x=0 in K
P1 = 10; % Thermal flux at x=L in W/m^2

[K, F] = ProblemAssemblyPrototype(E, A, T0, f, P1, h);

fprintf('K:\n')
disp(K)
fprintf('F:\n')
disp(F)