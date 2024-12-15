function [phi, phix, phixx, T] = init_fespace() %9-by-2, 9-by-1, lambda1, lambda2, lambda3, x_1, x_2, x_3, y_1, y_2, y_3
syms lambda1 lambda2 lambda3 x1 x2 x3 y1 y2 y3 real
xi1 = x2 - x3;
xi2 = x3 - x1;
xi3 = x1 - x2;
eta1 = y2 - y3;
eta2 = y3 - y1;
eta3 = y1 - y2;
l1 = eta1^2 + xi1^2;
l2 = eta2^2 + xi2^2;
l3 = eta3^2 + xi3^2;
E1 = (l3 - l2) / l1;
E2 = (l1 - l3) / l2;
E3 = (l2 - l1) / l3;
T2 = det([1, x1, y1; 1, x2, y2; 1, x3, y3]);
lambda1x = eta1 / T2;
lambda2x = eta2 / T2;
lambda3x = eta3 / T2;
lambda1y = -xi1 / T2;
lambda2y = -xi2 / T2;
lambda3y = -xi3 / T2;
phi = [
    1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, xi2, -xi3, 0, 0, 0, 0;
    0, 0, 0, 0, 0, xi3, -xi1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, xi1, -xi2;
    0, 0, 0, eta2, -eta3, 0, 0, 0, 0;
    0, 0, 0, 0, 0, eta3, -eta1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, eta1, -eta2
] * [
    -1/2*(E2-E3), 0, 0, 3/2*(3+E2), 3/2*(3-E3), 0, 0, 0, 0, 0;
    1/2*(1-2*E1-E3), 1, 0, -3/2*(1-E1), 3/2*(E1+E3), 3, 3, 0, 0, 3*(1-E1);
    1/2*(1+2*E1+E2), 0, 1, -3/2*(E1+E2), -3/2*(1+E1), 0, 0, 3, 3, 3*(1+E1);
    -1/4*(1+E2), 0, 0, 1/4*(5+3*E2), 1/2, 0, 0, 0, 0, 0;
    -1/4*(1-E3), 0, 0, 1/2, 1/4*(5-3*E3), 0, 0, 0, 0, 0;
    1/4*(1-E3), 0, 0, -1/2, -1/4*(1-3*E3), 1, 0, 0, 0, 1;
    -1/2*E1, 0, 0, -1/4*(1-3*E1), 1/4*(1+3*E1), 0, 1, 0, 0, 1/2*(1-3*E1);
    1/2*E1, 0, 0, 1/4*(1-3*E1), -1/4*(1+3*E1), 0, 0, 1, 0, 1/2*(1+3*E1);
    1/4*(1+E2), 0, 0, -1/4*(1+3*E2), -1/2, 0, 0, 0, 1, 1
] * [lambda1^3; lambda2^3; lambda3^3; lambda1^2*lambda3; lambda1^2*lambda2; lambda2^2*lambda1; lambda2^2*lambda3; lambda3^2*lambda2; lambda3^2*lambda1; lambda1*lambda2*lambda3];
phix = lambda1x * diff(phi, lambda1) + lambda2x * diff(phi, lambda2) + lambda3x * diff(phi, lambda3);
phiy = lambda1y * diff(phi, lambda1) + lambda2y * diff(phi, lambda2) + lambda3y * diff(phi, lambda3);
phixx = lambda1x * diff(phix, lambda1) + lambda2x * diff(phix, lambda2) + lambda3x * diff(phix, lambda3);
phiyy = lambda1y * diff(phiy, lambda1) + lambda2y * diff(phiy, lambda2) + lambda3y * diff(phiy, lambda3);
phi = matlabFunction(phi); 
phix = matlabFunction([phix, phiy]);
phixx = matlabFunction(phixx + phiyy);
T = matlabFunction(det([1, (x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3; 1, x2, y2; 1, x3, y3]) / 2);
end