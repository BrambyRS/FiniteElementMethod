function [K, F] = ProblemAssemblyPrototype(E, A, u0, f, t, h)
n_el = length(h);

iN1 = 1; % Integral of N1 over the xi range [-1, 1]
iN2 = 1; % Integral of N2 over the xi range [-1, 1]
dN1 = -0.5; % Derivative of N1 with respect to xi
dN2 = 0.5; % Derivative of N2 with respect to xi

left_bc = 'neumann';
right_bc = 'neumann';

n_dof = n_el+1;

if strcmp(left_bc, 'dirichlet')
    n_dof = n_dof - 1;
end
if strcmp(right_bc, 'dirichlet')
    n_dof = n_dof - 1;
end

K = zeros(n_dof, n_dof);
F = zeros(n_dof, 1);

switch left_bc
    case 'dirichlet'
        K(1, 1) = K(1, 1) + 4*E*A/h(1)*dN2*dN2;
        F(1) = f*A*h(1) - 4*E*A/h(1)*dN1*dN2*u0;
    case 'neumann'
        F(1) = F(1) + t*A;
    otherwise
        error('Left boundary condition type not recognized.');
end

for i = 1:n_dof-1
    K_local = 4*E*A/h(i) * [dN1*dN1, dN1*dN2; dN2*dN1, dN2*dN2];
    F_local = f*A*h(i) * [iN1; iN2];

    K(i:i+1, i:i+1) = K(i:i+1, i:i+1) + K_local;
    F(i:i+1) = F(i:i+1) + F_local;
end

switch right_bc
    case 'dirichlet'
        K(end, end) = K(end, end) + 4*E*A/h(end)*dN1*dN1;
        F(end) = F(end) + f*A*h(end) - 4*E*A/h(end)*dN1*dN2*u0;
    case 'neumann'
        F(end) = F(end) + t*A;
    otherwise
        error('Right boundary condition type not recognized.');
end

end
