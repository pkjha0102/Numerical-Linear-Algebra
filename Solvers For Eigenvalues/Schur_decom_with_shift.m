clc

A = [17 24 1 8 15; 23 5 7 14 16; 4 6 13 20 22; 10 12 19 21 3; 11 18 25 2 8];
maxIter = 20;

[my_eig, eig_inbuilt, iter_used, norm_error, error_mat] = SelfQRIter(A, maxIter)
%Compare eigenvectors from inbuilt commands

function [u, eig_inbuilt, k, norm_E, E] = SelfQRIter(A, maxIter)
[P, D] = eig(A);
eig_inbuilt = sort(diag(D));
[n, n] = size(A);
for i = 1:maxIter
    shift = A(n, n);
    [q, r] = qr(A - shift * eye(n));
    A = (r*q) + shift * eye(n);    
    %Eigenvalues lie on diagonals of upper triangular matrix
    u = sort(diag(A));
    k = i;
    E(:, i) = abs(eig_inbuilt - u);
end
norm_E = norm(E(:, end));
end