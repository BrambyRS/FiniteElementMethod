function x = QRSolvePrototype(A, b)
% QRSOLVEPROTOTYPE Solve the linear system Ax = b using QR decomposition.
% This is a prototype function to work out the algorithm.
%
% Input:
%   A - An m x n matrix (m >= n).
%   b - An m x 1 vector.
% Output:
%   x - An n x 1 vector that solves Ax = b in the least squares sense.

    [m, n] = size(A);
    assert(m >= n, 'Matrix A of dimension m x n must have m >= n.');
    assert(length(b) == m, 'Vector b must have length m.');

    [Q, R] = QRPrototype(A);
    y = Q' * b;
    x = zeros(n, 1);
    for i = n:-1:1
        sum = 0;
        for j = i+1:n
            sum = sum + R(i, j) * x(j);
        end
        x(i) = (y(i) - sum) / R(i, i);
    end
end
