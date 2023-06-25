clc

%n = input('Size: \n');
%A = rand(n, n);
A = [1 -6 9; 6 2 3; 9 3 2];
b = [0; 5; 0];

maxIter = input('Maximum number of iterations: \n');
tol = 1e-8;

[solution, res_norm, iter_used] = SelfNSD(A, b, maxIter, tol)
% Actual solution by backslash command
solution_backslash = A\b

function [x, z, t] = SelfNSD(A, b, maxIter, tol)
[n, n] = size(A);
x = ones(n, 1);
r = b - A*x;
p = A'*r;
for i = 2:maxIter
    a = dot(p, p) / dot(A*p, A*p);
    x = x + a*p;
    r = r - a*A*p;
    p = A'*r;
    z = norm(r);
    t = i;
    if z < tol
        break;
    end
end
end