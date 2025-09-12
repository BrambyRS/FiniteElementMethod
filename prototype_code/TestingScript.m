A = [2, -1, 0; 1, 2, 1; 0, 1, 3];
b = [1; 4; 2];

[Q, R] = QRPrototype(A);

x = QRSolvePrototype(A, b);
b_test = A*x;

E = 45; % Steel thermal conductivity W/(m.K)
A = 5e-5; % Cross-sectional area in m^2
L = 0.076; % Length in m

n_el = 10; % Number of elements
h = ones(n_el, 1) * (L / n_el);

f = -10; % No internal heat generation
T0 = 293.15; % Temperture at x=0 in K
P1 = 1000; % Thermal flux at x=L in W/m^2

[K, F] = ProblemAssemblyPrototype(E, A, T0, f, P1, h);
d = QRSolvePrototype(K, F);
T = [T0; d];

% Plot a figure
x = linspace(0, L, n_el+1);
hFig = figure(1);
clf();
hAx = defaultAxisSettings(axes(hFig));
plot(hAx, x, T - 273.15, '-o');
hAx.XLabel.String = 'Position [m]';
hAx.YLabel.String = 'Temperature [$^\circ$C]';
hAx.Title.String = 'Temperature Distribution';
