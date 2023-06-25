A = [1 2 3; 4 5 6; 7 8 9];
tol = 1e-8
[U_svd, S_svd, V_svd] = svd(A);

[U, S, V] = SelfSVD(A, tol)

x = norm(A - U*S*V')

function [U_economic, S, V_economic] = SelfSVD(A, tol)
%V1 is eigenvector of A'*A
%D1 is diagonal matrix with eigenvalues on diagonal
[V, D] = eig(A' * A)
D_rooted = real(sqrt(D))
temp = diag(D_rooted);
[D_sorted, ind] = sort(temp,"descend");

%S is sorted diagonal matrix
S = diag(D_sorted);
% Remove zero rows
S(max(S)<tol,:)=[];
% Remove zero columns
S(:,max(S)<tol)=[];

r = size(S);
rank_svd = r(1,1);
%V is rearranged(sorted) V1
V_sorted = V(:, ind);
for i = 1:rank_svd
    V_economic(:, i) = V_sorted(:, i);
end
for i = 1:rank_svd
    U_economic(:, i) = (A * V_economic(:, i)) / S(i, i);
end
end