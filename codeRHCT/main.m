clear;
epsilon = 1;
syms x y real
u = (cos(2*pi*x)-1)*(cos(4*pi*y)-1);
ux = diff(u, x);
uy = diff(u, y);
uxx = diff(ux, x) + diff(uy, y);
f = epsilon^2 * (diff(diff(uxx, x), x) + diff(diff(uxx, y), y)) - uxx;
ux = matlabFunction([ux, uy]); %1-by-2, x, y
ux = @(lambda, v) ux(v(1, :) * lambda, v(2, :) * lambda); 
uxx = matlabFunction(uxx); %x, y
uxx = @(lambda, v) uxx(v(1, :) * lambda, v(2, :) * lambda);
f = matlabFunction(f); %x, y
f = @(lambda, v) f(v(1, :) * lambda, v(2, :) * lambda);
[phi, phix, phixx, T] = init_fespace();
phi = @(lambda, v) phi(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
phix = @(lambda, v) phix(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
phixx = @(lambda, v) phixx(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
T = @(v) T(v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
error = zeros(1, 4);
rate = zeros(1, 4);
parfor i = 1: 4
	error(i) = FEM(2^(-i-2), epsilon, ux, uxx, f, phi, phix, phixx, T);
end
for i = 2: 4
    rate(i) = log2(error(i-1)/error(i));
end