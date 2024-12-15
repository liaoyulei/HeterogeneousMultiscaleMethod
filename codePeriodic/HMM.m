function a_H= HMM(a, iota, delta, center, p, meshes, trans, phix, phixx, T)
omega = @(x, y) (1 + cos(2*pi*(x-center(1))/delta)) * (1 + cos(2*pi*(y-center(2))/delta));
omega = @(lambda, v) omega(v(1, :) * lambda, v(2, :) * lambda);
%omega = @(lambda, v) 1;
wg = [0.144315607677787, 0.103217370534718, 0.103217370534718, 0.103217370534718, 0.032458497623198, 0.032458497623198, 0.032458497623198, 0.095091634267285, 0.095091634267285, 0.095091634267285, 0.027230314174435, 0.027230314174435, 0.027230314174435, 0.027230314174435, 0.027230314174435, 0.027230314174435];
ld = [
    1/3, 1/3, 1/3; 
    0.170569307751760, 0.170569307751760, 0.658861384496480;
    0.170569307751760, 0.658861384496480, 0.170569307751760;
    0.658861384496480, 0.170569307751760, 0.170569307751760;
    0.050547228317031, 0.050547228317031, 0.898905543365938;
    0.050547228317031, 0.898905543365938, 0.050547228317031;
    0.898905543365938, 0.050547228317031, 0.050547228317031;
    0.459292588292723, 0.459292588292723, 0.081414823414554;
    0.459292588292723, 0.081414823414554, 0.459292588292723;
    0.081414823414554, 0.459292588292723, 0.459292588292723;
    0.263112829634638, 0.728492392955404, 0.008394777409958;
    0.263112829634638, 0.008394777409958, 0.728492392955404;
    0.728492392955404, 0.263112829634638, 0.008394777409958;
    0.728492392955404, 0.008394777409958, 0.263112829634638;
    0.008394777409958, 0.728492392955404, 0.263112829634638;
    0.008394777409958, 0.263112829634638, 0.728492392955404
]';
vertices = [delta * p + center, [ones(1, size(p, 2)), zeros(1, size(p, 2)); zeros(1, size(p, 2)), ones(1, size(p, 2))]];
A = zeros(size(vertices, 2)); %Nv-by-Nv
b = zeros(size(vertices, 2), 2); %Nv-by-2
c = zeros(size(b)); %Nv-by-2
for i = 1: size(meshes, 2)
	v = vertices(:, meshes(1: 3, i));
	jacobi = T(v);
    local = zeros(9, 9);
	for j = 1: size(wg, 2)
		local = local + wg(j) * jacobi * (iota^2 * phixx(ld(:, j), v) * phixx(ld(:, j), v)' + phix(ld(:, j), v) * a(ld(:, j), v) * phix(ld(:, j), v)'); 
	end
    A(trans(meshes(:, i)), trans(meshes(:, i))) = A(trans(meshes(:, i)), trans(meshes(:, i))) + local;
    b(trans(meshes(:, i)), :) = b(trans(meshes(:, i)), :) - local *  vertices(:, meshes(:, i))';
end
A = sparse(A);
free = trans == 1: size(vertices, 2);
free(1) = 0;
c(free, :) = A(free, free) \ b(free, :);
c = c(trans, :) + vertices';
a_H = zeros(2);
v_ave = zeros(2);
for i = 1: size(meshes, 2)
	v = vertices(:, meshes(1: 3, i));
	jacobi = T(v);
	for j = 1: size(wg, 2)
		a_H = a_H + wg(j) * jacobi * omega(ld(:, j), v) * a(ld(:, j), v) * phix(ld(:, j), v)' * c(meshes(:, i), :);
        v_ave = v_ave + wg(j) * jacobi * omega(ld(:, j), v) * phix(ld(:, j), v)' * c(meshes(:, i), :);
	end
end
a_H = a_H * inv(v_ave);
end
