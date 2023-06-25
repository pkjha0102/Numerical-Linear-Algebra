n = input("");
A = rand(10, n);

[Q, R] = HouseSelf(A)
norm(A - Q*R)

function [Q,R] = HouseSelf(A)
[m,n] = size(A);
Q = eye(m); % Orthogonal transform so far
R = A; % Transformed matrix so far
for j = 1:n
% -- Find H = I-tau*w*w’ to put zeros below R(j,j)
normx = norm(R(j:end,j));
s = -sign(R(j,j));
u1 = R(j,j) - s*normx;
w = R(j:end,j)/u1;
w(1) = 1;
tau = -s*u1/normx;
% -- R := HR, Q := QH
R(j:end,:) = R(j:end,:)-(tau*w)*(w'*R(j:end,:));
Q(:,j:end) = Q(:,j:end)-(Q(:,j:end)*w)*(tau*w)';
end
end