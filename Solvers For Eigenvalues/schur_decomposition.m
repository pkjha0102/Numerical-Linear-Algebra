clc
%n = input('size: ');
A = [17 24 1 8 15; 23 5 7 14 16; 4 6 13 20 22; 10 12 19 21 3; 11 18 25 2 8];

tol = 1e-8;
maxIter = 1000;

[eigenvalues, iter_used] = SelfQRIter(A, maxIter, tol)
%Compare eigenvectors from inbuilt commands
[P, D] = eig(A);
eig_inbuilt = diag(D)

function [u, k] = SelfQRIter(A, maxIter, tol)
[q, r] = qr(A);
A1 = r*q;
for i = 2:maxIter
    [q, r] = qr(A1);
    A = r*q;
    [q1, r1] = qr(A);
    A1 = r1*q1;
    %Checking convergence
    if norm(tril(A1) - tril(A), 'fro') < tol * norm(tril(A1), 'fro')
        break;
    end
    %Eigenvalues lie on diagonals of upper triangular matrix
    u = diag(A1);
    k = i;
end
end