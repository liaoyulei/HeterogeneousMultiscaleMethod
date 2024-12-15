function [phix, phixx, T] = init_fespace() %9-by-2, 9-by-4
syms lambda1 lambda2 lambda3 x1 x2 x3 y1 y2 y3 real
xi1 = x2 - x3;
xi2 = x3 - x1;
xi3 = x1 - x2;
eta1 = y2 - y3;
eta2 = y3 - y1;
eta3 = y1 - y2;
inner11 = eta1^2 + xi1^2;
inner12 = eta1 * eta2 + xi1 * xi2;
inner13 = eta1 * eta3 + xi1 * xi3;
inner22 = eta2^2 + xi2^2;
inner23 = eta2 * eta3 + xi2 * xi3;
inner33 = eta3^2 + xi3^2;
T2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
lambda1x = eta1 / T2;
lambda2x = eta2 / T2;
lambda3x = eta3 / T2;
lambda1y = -xi1 / T2;
lambda2y = -xi2 / T2;
lambda3y = -xi3 / T2;
bK = lambda1 * lambda2 * lambda3;
p = bK * [2 * lambda1 - 1; 2 * lambda2 - 1; 2 * lambda3 - 1];
phi = [
    lambda1^2 * (3 - 2 * lambda1) + 6 * inner12 / inner22 * p(2) + 6 * inner13 / inner33 * p(3);
    lambda2^2 * (3 - 2 * lambda2) + 6 * inner23 / inner33 * p(3) + 6 * inner12 / inner11 * p(1);
    lambda3^2 * (3 - 2 * lambda3) + 6 * inner13 / inner11 * p(1) + 6 * inner23 / inner22 * p(2);
    lambda1^2 * (xi2 * lambda3 - xi3 * lambda2) + 2* (xi2 + xi3) * bK * (lambda2 - lambda3) + 3 * xi2 * inner12 / inner22 * p(2) - 3 * xi3 * inner13 / inner33 * p(3);
    lambda2^2 * (xi3 * lambda1 - xi1 * lambda3) + 2 * (xi3 + xi1) * bK * (lambda3 - lambda1) + 3 * xi3 * inner23 / inner33 * p(3) - 3 * xi1 * inner12 / inner11 * p(1);
    lambda3^2 * (xi1 * lambda2 - xi2 * lambda1) + 2 * (xi1 + xi2) * bK * (lambda1 - lambda2) + 3 * xi1 * inner13 / inner11 * p(1) - 3 * xi2 * inner23 / inner22 * p(2);
    lambda1^2 * (eta2 * lambda3 - eta3 * lambda2) + 2 * (eta2 + eta3) * bK * (lambda2 - lambda3) + 3 * eta2 * inner12 / inner22 * p(2) - 3 * eta3 * inner13 / inner33 * p(3);
    lambda2^2 * (eta3 * lambda1 - eta1 * lambda3) + 2 * (eta3 + eta1) * bK * (lambda3 - lambda1) + 3 * eta3 * inner23 / inner33 * p(3) - 3 * eta1 * inner12 / inner11 * p(1);
    lambda3^2 * (eta1 * lambda2 - eta2 * lambda1) + 2 * (eta1 + eta2) * bK * (lambda1 - lambda2) + 3 * eta1 * inner13 / inner11 * p(1) - 3 * eta2 * inner23 / inner22 * p(2);
];
phix = lambda1x * diff(phi, lambda1) + lambda2x * diff(phi, lambda2) + lambda3x * diff(phi, lambda3);
phiy = lambda1y * diff(phi, lambda1) + lambda2y * diff(phi, lambda2) + lambda3y * diff(phi, lambda3);
phixx = lambda1x * diff(phix, lambda1) + lambda2x * diff(phix, lambda2) + lambda3x * diff(phix, lambda3);
phixy = lambda1y * diff(phix, lambda1) + lambda2y * diff(phix, lambda2) + lambda3y * diff(phix, lambda3);
phiyy = lambda1y * diff(phiy, lambda1) + lambda2y * diff(phiy, lambda2) + lambda3y * diff(phiy, lambda3);
phix = matlabFunction([phix, phiy]);
phixx = matlabFunction([phixx, phixy, phixy, phiyy]);
T = matlabFunction(T2 / 2);
end