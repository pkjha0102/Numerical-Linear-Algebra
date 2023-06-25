clc

%n = input('Size: \n');
%A = rand(n, n);
%A = [4 -1 0 0 0 0; -1 4 -1 0 0 0; 0 -1 4 -1 0 0; 0 0 -1 4 -1 0; 0 0 0 -1 4 -1; 0 0 0 0 -1 4];
%[V, D] = eig(A);
%b = V(:, 4);
%b = [0; 5; 0; 6; -2; 6];
A = [1 -6 9; 6 2 3; 9 3 2];
b = [0; 5; 0];

maxIter = input('Maximum iterations: \n');
tol = 1e-8;

[solution, res_norm, iter_used] = SelfMR(A, b, maxIter, tol)
% Actual solution by backslash command
solution_backslash = A\b

function [x, z, t] = SelfMR(A, b, maxIter, tol)
[n, n] = size(A);
x = ones(n, 1);
r = b - A*x;
p = A*r;
for i = 2:maxIter
    a = dot(p, r) / dot(p, p);
    x = x + a*r;
    r = r - a*p;
    p = A*r;
    z = norm(r);
    t = i;
    if z < tol
        break;
    end
end
end