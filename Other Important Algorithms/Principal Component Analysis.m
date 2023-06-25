N = 100;
x1 = rand(1, N);
x2 = 0.4*rand(1, N);

A = [x1;x2];
R = [cos(deg2rad(60)) -sin(deg2rad(60)); sin(deg2rad(60)), cos(deg2rad(60))];
D = R*A;

x = D(2,:);
y = D(1,:);
plot(x,  y, '.');
xlim([-4, 4])
ylim([-4, 4])

hold on

Q = (D*D')/(N-1);
[V, D] = eig(Q);
temp = diag(D);
[temp_sorted, idx] = sort(temp, "descend");
D_sorted = diag(temp_sorted)

V_sorted = V(:, idx)
u1 = V_sorted(:, 1)
plot(u1);