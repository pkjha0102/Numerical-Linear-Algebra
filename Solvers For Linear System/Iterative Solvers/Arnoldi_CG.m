clc

%n = input("Size: \n");
m = input("Size of subspace: \n");
max_iter = input("Maximum iterations: \n");
tol = 1e-8;
%A = rand(n, n);
A = [5 7 6 5; 7 10 8 7; 6 8 10 9; 5 7 9 10];
n = 4;
b = rand(n, 1);

[H, V] = ArnoldiSelf(A, b, n, m)
G = V * H * V'

[mysol, iter_used] = CGSelf(A, b, n, max_iter,tol)
% Compare with backslash solution
sol_inbuilt = A\b
norm_res = norm(mysol - sol_inbuilt)

% Function to implement Arnoldi
function [H, V] = ArnoldiSelf(A, b, n, m)
xo = ones(n,1);
ro = b - A*xo;
V(:, 1) = ro / norm(ro);
for j = 1:m
    w = A*V(:, j);
    for i = 1:j
        H(i, j) = dot(w, V(:, i));
        w = w - H(i, j) * V(:, i);
    end
    H(j+1, j) = norm(w);
    if H(j+1, j) == 0
        break;
    end
    V(:, j+1) = w / H(j+1, j);
end
H = H(1:m, :);
V = V(:, 1:m);
end

% Function to implement CG utilising Arnoldi
function [x, iter] = CGSelf(A, b, n, max_iter, tol)
x = ones(n, 1);
r = b - A*x;
p = r;
for j = 1:max_iter
    a = dot(r, r) / dot(A*p, p);
    x = x + a*p;
    r_new = r - a*A*p;
    beta = dot(r_new, r_new) / dot(r, r);
    p = r_new + beta * p;
    if norm(r_new) < tol
        iter = j;
        break;
    end
    r = r_new; %update 'r' if not converge
end
end
