function [phix, T] = init_fespace() %3-by-2, 3-by-1, x_1, x_2, x_3, y_1, y_2, y_3
syms lambda1 lambda2 lambda3 x1 x2 x3 y1 y2 y3 real
xi1 = x2 - x3;
xi2 = x3 - x1;
xi3 = x1 - x2;
eta1 = y2 - y3;
eta2 = y3 - y1;
eta3 = y1 - y2;
T2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
lambda1x = eta1 / T2;
lambda2x = eta2 / T2;
lambda3x = eta3 / T2;
lambda1y = -xi1 / T2;
lambda2y = -xi2 / T2;
lambda3y = -xi3 / T2;
phi = [lambda1; lambda2; lambda3];
phix = lambda1x * diff(phi, lambda1) + lambda2x * diff(phi, lambda2) + lambda3x * diff(phi, lambda3);
phiy = lambda1y * diff(phi, lambda1) + lambda2y * diff(phi, lambda2) + lambda3y * diff(phi, lambda3);
phix = matlabFunction([phix, phiy]);
T = matlabFunction(T2 / 2);
end