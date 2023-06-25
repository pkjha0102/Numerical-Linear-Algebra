A = [6 1 1 1 1; 1 7 1 1 1; 1 1 8 1 1; 1 1 1 9 1; 1 1 1 1 10];
b = [-10;-6;0;8;18];
tol = 10^-8;
maxiter = 50;

[myans, error, iter] = GSIterSelf(A, b, tol, maxiter)

exact_ans = A\b

function [u, v, t] = GSIterSelf(A, b, tol, maxiter);
[n, n] = size(A);
%Assume initially the solution is x = [1;1;1;1;1]
x = zeros(n, 1);
p = diag(diag(A)) - (-1)*tril(A, -1)
sol_mat = zeros(n, maxiter); %sol_mat contains solution vector of each staep as columns
sol_mat(:, 1) = x;
err_mat = zeros(n, maxiter); %err_mat contains error as columns
err_mat(:, 1) = A*x - b;
for i = 2 : maxiter
    sol_mat(:, i) = (p^-1)*(p-A)*(sol_mat(:, i-1)) + (p^-1)*b;
    err_mat(:, i) = A*(sol_mat(:,i)) - b;
    if (norm(sol_mat(:,i) - sol_mat(:,i-1), inf) / norm(sol_mat(:,i-1), inf)) < tol
        u = sol_mat(:, i);
        v = err_mat(:, i);
        t = i;
        break;
    end
end
plot_mat = zeros(t, 2);
for j = 1 : t
    plot_mat(j, 1) = j;
    plot_mat(j, 2) = norm(err_mat(:, j), inf);
end
plot(plot_mat(:, 1), plot_mat(:, 2));
end