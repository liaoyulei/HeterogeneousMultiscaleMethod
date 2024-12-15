function error = FEM(h, ux, f, a, vertices, meshes, free, phi, phix, T) %3-by-1, lambda; 3-by-2, v; 1-by-1, v
wg = [1/3, 1/3, 1/3];
ld = [
    1, 0, 0;
    0, 1, 0;
    0, 0, 1
]';
A = zeros(size(vertices, 2)); %Nv-by-Nv
b = zeros(size(vertices, 2), 1); %Nv-by-1
c = zeros(size(vertices, 2), 1); %Nv-by-1
for i = 1: size(meshes, 2)
	v = vertices(:, meshes(:, i));
	jacobi = T(v);
	for j = 1: size(wg, 2)
		A(meshes(:, i), meshes(:, i)) = A(meshes(:, i), meshes(:, i)) + wg(j) * jacobi * phix(v) * a(:, :, meshes(j, i)) * phix(v)'; 
		b(meshes(:, i)) = b(meshes(:, i)) + wg(j) * jacobi * f(ld(:, j), v) * phi(ld(:, j));	
	end
end
A = sparse(A); %Nv-by-Nv
c(free) = A(free, free) \ b(free);
error = 0;
normu = 0;
for i = 1: size(meshes, 2)
	v = vertices(:, meshes(:, i));
	jacobi = T(v);
	for j = 1: size(wg, 2)
		ex =  ux(ld(:, j), v) - c(meshes(:, i))' * phix(v); %1-by-2
    	error = error + wg(j) * jacobi * ex * ex';
    	normu = normu + wg(j) * jacobi * ux(ld(:, j), v) * ux(ld(:, j), v)';
	end
end
error = (error / normu)^(1/2);
end