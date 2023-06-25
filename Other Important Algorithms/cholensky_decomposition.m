A = [2 3; 4 5];
b = [8; 14];

[x] = LSSNEqn(A, b)
d = A\b

% Function to solve AX = B using cholesky decomposition of A
function [x] = LSSNEqn(A, b);
b1 = A' * b;
A1 = A' * A
U = chol(A1);
L = U';
myobj = functionscontainer;
y = myobj.FdSubs(L, b1); % Solve LY = B
x = myobj.BdSubs(U, y);  % Solve UX = Y
end