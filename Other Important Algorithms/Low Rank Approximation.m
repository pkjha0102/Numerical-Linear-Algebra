A = rand(8, 10);
k = input('k: ');

A1 = svds(A, k);

[U_k, S_k, V_k] = SelfLRA(A, k)

x_2 = norm(A - U_k * S_k * V_k')
x_f = norm(A - U_k * S_k * V_k', "fro")

function [U_k, S_k, V_k] = SelfLRA(A, k)
%V1 is eigenvector of A'*A
%D1 is diagonal matrix with eigenvalues on diagonal
[V, D] = eig(A' * A);
D_rooted = real(sqrt(D));
temp = diag(D_rooted);
[D_sorted, ind] = sort(temp,"descend")

%S1 is sorted diagonal matrix
S = diag(D_sorted);
%Choosing k largest singular values only
S_k = S(1:k, 1:k);
%V_k is first k columns of sorted V
V_sorted = V(:, ind);
V_k = V_sorted(:, 1:k);
for i = 1:k
    U_k(:, i) = (A * V_k(:, i)) / S_k(i, i);
end
end