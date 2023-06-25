clc

%n = input("Size: \n");
m = input("Size of subspace: \n");
%A = rand(n, n);
A = [5 7 6 5; 7 10 8 7; 6 8 10 9; 5 7 9 10];
n = 4;
b = ones(n, 1);

[H, V] = ArnoldiSelf(A, b, n, m)
F = V' * A *V

mysol = FOMSelf(A, b, n, m, V, H)
% Compare with backslash solution
sol_inbuilt = A\b
norm_res = norm(sol_inbuilt - mysol) %error in both solutions

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

% Function to implement FOM utilising Arnoldi
function [x] = FOMSelf(A, b, n, m, V, H)
xo = ones(n, 1);
ro = b - A*xo;
beta = norm(ro);
e = eye([m, 1]);
y = H\(beta * e);
x = xo + V*y;
r = b - A*x;
end