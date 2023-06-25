clc

m = input("Size of subspace: \n")
A = [2 -1 0 0; -1 2 -1 0; 0 -1 2 -1; 0 0 -1 2];
[n, n] = size(A);
b = [1;0;0;0];
tol = 1e-6;
xo = zeros(n,1);

[V, H] = ArnoldiSelf(A, b, xo, n, m, tol)

mysol = SelfFOM(A, b, xo, m, V, H)
% Compare with backslash solution
sol_inbuilt = A\b
norm_res = norm(sol_inbuilt - mysol) %error in both solutions

% Function to implement Arnoldi
function [V, H] = ArnoldiSelf(A, b, xo, n, m, tol)
ro = b - A*xo;
V(:, 1) = ro / norm(ro);
for j = 1:m
    w = A*V(:, j);
    for i = 1:j
        H(i, j) = dot(w, V(:, i));
        w = w - H(i, j) * V(:, i);
    end
    H(j+1, j) = norm(w);
    V(:, j+1) = w / H(j+1, j);
    if H(j+1, j) < tol
        break;
    end
end
V = V(:, 1:m);
end

% Function to implement FOM utilising Arnoldi
function [x] = SelfFOM(A, b, xo, m, V, H)
ro = b - A*xo;
beta = norm(ro);
e = eye([m+1, 1]);
y = H\(beta * e);
x = xo + V*y;
r = b - A*x;
end