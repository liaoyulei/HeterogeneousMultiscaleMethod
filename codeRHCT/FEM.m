function error = FEM(h, epsilon, ux, uxx, f, phi, phix, phixx, T)
wg = [0.109951743655322, 0.109951743655322, 0.109951743655322, 0.223381589678011, 0.223381589678011, 0.223381589678011];
ld = [
    0.816847572980459, 0.091576213509771, 0.091576213509771;
    0.091576213509771, 0.816847572980459, 0.091576213509771;
    0.091576213509771, 0.091576213509771, 0.816847572980459;
    0.108103018168070, 0.445948490915965, 0.445948490915965;
    0.445948490915965, 0.108103018168070, 0.445948490915965;
    0.445948490915965, 0.445948490915965, 0.108103018168070;
]';
lambda = [ld(1, :) / 3; ld(1, :) / 3 + ld(2, :); ld(1, :) / 3 + ld(3, :)];
[vertices, meshes, free] = init_mesh(h);
A = zeros(size(vertices, 2)); %Nv-by-Nv
b = zeros(size(vertices, 2), 1); %Nv-by-1
c = zeros(size(vertices, 2), 1); %Nv-by-1
for i = 1: size(meshes, 2)
	for j = 0: 2
		idx = [mod([j, j+1, j+2], 3) + 1, mod([j, j+1, j+2], 3) + 4, mod([j, j+1, j+2], 3) + 7];
		v = vertices(:, meshes(idx(1: 3), i));
		jacobi = T(v);
		for k = 1: size(wg, 2)
			A(meshes(idx, i), meshes(idx, i)) = A(meshes(idx, i), meshes(idx, i)) + wg(k) * jacobi * (epsilon^2 * phixx(lambda(:, k), v) * phixx(lambda(:, k), v)' + phix(lambda(:, k), v) * phix(lambda(:, k), v)'); 
			b(meshes(idx, i)) = b(meshes(idx, i)) + wg(k) * jacobi * f(lambda(:, k), v) * phi(lambda(:, k), v);	
		end
	end
end
A = sparse(A); %Nv-by-Nv
c(free) = A(free, free) \ b(free);
error = 0;
normu = 0;
for i = 1: size(meshes, 2)
	for j = 0: 2
		idx = [mod([j, j+1, j+2], 3) + 1, mod([j, j+1, j+2], 3) + 4, mod([j, j+1, j+2], 3) + 7];
		v = vertices(:, meshes(idx(1: 3), i));
		jacobi = T(v);
		for k = 1: size(wg, 2)
			ex =  ux(lambda(:, k), v) - c(meshes(idx, i))' * phix(lambda(:, k), v); %1-by-2
	    	exx = uxx(lambda(:, k), v) - c(meshes(idx, i))' * phixx(lambda(:, k), v); %1
	    	error = error + wg(k) * jacobi * (epsilon^2 * exx^2 + ex * ex');
	    	normu = normu + wg(k) * jacobi * (epsilon^2 * uxx(lambda(:, k), v)^2 + ux(lambda(:, k), v) * ux(lambda(:, k), v)');
		end
	end
end
error = (error / normu)^(1/2);
end