clc

%n = input('Size: \n');
%A = rand(n, n);
A = [5 7 6 5; 7 10 8 7; 6 8 10 9; 5 7 9 10];
b = ones(n, 1);

maxIter = input('Maximum number of iterations: \n');
tol = 1e-8;

[solution, residual, res_norm, iter_used] = SelfSD(A, b, maxIter, tol)
solution_backslash = A\b
norm_sol_bckslsh = norm(b - A*solution_backslash)

function [x, r, z, t] = SelfSD(A, b, maxIter, tol)
[n, n] = size(A);
x = ones(n, 1);
r = b - A*x;
p = A*r;
for i = 2:maxIter
    a = dot(r, r) / dot(r, p);
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