function [Q, R] = QRPrototype(A)
% QRPROTOTYPE Compute the QR decomposition of matrix A using Householder reflections.
% This is a prototype function to work out the algorithm.
%
% Input:
%   A - An m x n matrix to decompose (m >= n).
% Output:
%   Q - An m x m orthogonal matrix.
%   R - An m x n upper triangular matrix.

    [r, c] = size(A);
    assert(r >= c, 'Matrix A of dimension m x n must have m >= n.');

    Q = eye(r);
    R = A;

    for k = 1:c
        % Get column a
        a = R(k:r, k);
        e = zeros(length(a), 1);
        e(1) = 1;
        v = a + sign(a(1)) * norm(a) * e;
        H = eye(r);
        subH = eye(r-k+1) - 2 * (v*v')/(v'*v);
        H(k:r, k:r) = H(k:r, k:r) * subH;
        R = H * R;
        Q = Q * H';
    end
end
