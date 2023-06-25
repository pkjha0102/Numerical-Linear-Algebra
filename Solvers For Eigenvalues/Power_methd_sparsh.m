clc
%Timer start
tic

%n = input('size: ');
%sparse_A = rand(n, n);
A = [4 -1 0 0 0 0; -1 4 -1 0 0 0; 0 -1 4 -1 0 0; 0 0 -1 4 -1 0; 0 0 0 -1 4 -1; 0 0 0 0 -1 4];
[n, n] = size(A);
sparse_A = sparse(A);

tol = 1e-8;
maxIter = 100;
xo = sparse(ones(n,1));

[lambda_mat, eigenvalue, eigenvector, iter_used] = SelfPower(sparse_A, n, xo, maxIter, tol)
%Timer stop
toc
%Compare with inbuilt command
[V, D] = eig(A);
[eig_inbuilt, idx] = sort(diag(D), 'descend')

function [lambda_mat, lambda, x, k] = SelfPower(sparse_A, n, xo, maxIter, tol)
x = xo / norm(xo);
%Fill first values mannually to avoid indexing errors

% Matrix-vector multiplication step
p = sparse_A * x;

x = p / norm(p);
lambda = dot(x, p);
lambda_mat(1, 1) = 1;
lambda_mat(1, 2) = lambda;
for i = 2:maxIter
    % Matrix-vector multiplication step
    p = sparse_A*x;
    
    x = p / norm(p);
    lambda = dot(x, p);
    lambda_mat(i, 1) = i;
    lambda_mat(i, 2) = lambda;
    %Checking condition for convergence
    if abs(lambda_mat(i, 2) - lambda_mat(i-1, 2)) < tol * lambda_mat(i, 2)
        break
    end
    k = i;
end
end