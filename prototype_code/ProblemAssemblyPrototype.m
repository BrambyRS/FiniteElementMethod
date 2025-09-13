function [K, F] = ProblemAssemblyPrototype(E, A, u0, f, t, h)
n_el = length(h);
K = zeros(n_el, n_el);
F = zeros(n_el, 1);

iN1 = 1; % Integral of N1 over the xi range [-1, 1]
iN2 = 1; % Integral of N2 over the xi range [-1, 1]
dN1 = -0.5; % Derivative of N1 with respect to xi
dN2 = 0.5; % Derivative of N2 with respect to xi

for i = 1:n_el-1
    K_local = 4*E*A/h(i) * [dN1*dN1, dN1*dN2; dN2*dN1, dN2*dN2];
    F_local = f*A*h(i) * [iN1; iN2];

    K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + K_local;
    F(i:i+1) = F(i:i+1) + F_local;
end
K(1, 1) = K(1, 1) + 4*E*A/h(1)*dN2*dN2;
F(1) = F(1) - 4*E*A/h(1)*dN1*dN2*u0;
F(end) = F(end) + t*A;
end
