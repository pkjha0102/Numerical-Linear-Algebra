n = input("");
A = rand(n, n);

[Q1, R1] = qr(A);

[Q, R] = GramSelf(A)
norm(A - Q*R)

function [q, r] = GramSelf(A);
[n, m] = size(A);
A1 = A;
q = zeros(n, m);
r = zeros(m, m);
for i = 1:m
    for j = 1:i-1
        r(j, i) = dot(A1(:, i), q(:, j));
        A1(:, i) = A1(:, i) - r(j, i)*q(:, j); 
    end
    r(i, i) = norm(A1(:, i));
    q(:, i) = A1(:, i)/r(i, i);
end
end


% https://www.youtube.com/watch?v=58_31LW0S8E