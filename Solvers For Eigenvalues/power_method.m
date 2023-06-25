clc
%n = input('size: ');
%A = rand(n, n);

tol = 1e-8;
maxIter = 100;
xo = ones(n, 1);

[eigenvalue, lambda_mat, eigenvector, iter_used] = SelfPower(A, n, xo, maxIter, tol)
%Compare with inbuilt command
[V, D] = eig(A);

function [lambda, lambda_mat, x, k] = SelfPower(A, n, xo, maxIter, tol)
x = xo / norm(xo);
lambda_mat = zeros(n, 2);
%Fill first values mannually to avoid indexing errors
p = A*x;
x = p / norm(p);
lambda = dot(x, A*x);
lambda_mat(1, 1) = 1;
lambda_mat(1, 2) = lambda;
for i = 2:maxIter
    p = A*x;
    x = p / norm(p);
    lambda = dot(x, A*x);
    lambda_mat(i, 1) = i;
    lambda_mat(i, 2) = lambda;
    %Checking condition for convergence
    if abs(lambda_mat(i, 2) - lambda_mat(i-1, 2)) < tol * lambda_mat(i, 2)
        break
    end
    k = i;
end
lambda_mat = lambda_mat(1:k, :);
end