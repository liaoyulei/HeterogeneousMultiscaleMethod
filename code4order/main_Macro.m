clear;
my_cluster=parcluster('local');
mypar = parpool(my_cluster, my_cluster.NumWorkers);
epsilon = 2^(-5);
delta = 2^(-2);
h = 2^(-6);
syms x y real
a = matlabFunction([(20*pi+4*pi^2*cos(2*pi*x/epsilon)) / (2*pi-cos(2*pi*x/epsilon)), 2 + sin(2*pi*x) + cos(2*pi*y/epsilon); 3 + cos(2*pi*y) + sin(2*pi*x/epsilon), (22*pi+4*pi^2*sin(2*pi*y/epsilon)) / (2*pi-sin(2*pi*y/epsilon))]);
a = @(lambda, v) a(v(1, :) * lambda, v(2, :) * lambda);
u = sin(pi*x) * sin(pi*y);
gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
[p_h, meshes_h, free_h] = init_mesh(h); %2-by-Nv; 1-by-Nv, 9-by-Nt
[phix_h, phixx_h, T_h] = init_fespace();
phix_h = @(lambda, v) phix_h(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
phixx_h = @(lambda, v) phixx_h(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
T_h = @(v) T_h(v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
phi_H = @(lambda) lambda;
phix_H = @(v) [
	v(2, 2) - v(2, 3), v(1, 3) - v(1, 2); 
	v(2, 3) - v(2, 1), v(1, 1) - v(1, 3); 
	v(2, 1) - v(2, 2), v(1, 2) - v(1, 1)
] / det([1, 1, 1; v]);
T_H = @(v) det([1, 1, 1; v])/2;
error = zeros(1, 5);
rate = zeros(size(error));
for gamma = [1/4, 1, 4]
	fprintf("%.2f", gamma);
	if gamma < 1
		a11 = 4*pi*(2*pi^2 + 5) / sqrt(4*pi^2 - 1) - 4*pi^2;
		a22 = 2*pi*(4*pi^2+11) / sqrt(4*pi^2 - 1) - 4*pi^2;
	elseif gamma ==1
		a11 = 10;
		a22 = 11;
	else
		a11 = 4*pi^2 / ((5 + 2*pi^2) / sqrt(25 - pi^2) - 1);
		a22 = 4*pi^2 / ((11 + 4*pi^2) / sqrt(121 - 4*pi^2) - 1);
	end
	syms x y real
	a_bar = [a11, 2+sin(2*pi*x); 3+cos(2*pi*y), a22];
	ux = diff(u, x);
	uy = diff(u, y);
	f = -diff(a_bar(1, 1)*ux + a_bar(1, 2)*uy, x) - diff(a_bar(2, 1)*ux + a_bar(2, 2)*uy, y);
	ux = matlabFunction([ux, uy]); %1-by-2
	ux = @(lambda, v) ux(v(1, :) * lambda, v(2, :) * lambda);
	f = matlabFunction(f);
	f = @(lambda, v) f(v(1, :) * lambda, v(2, :) * lambda);
	for i = 1: size(error, 2)
		H = 2^(-i);
		iota = epsilon^gamma;
		meshes_H = generateMesh(model, 'Hmax', H, 'GeometricOrder', 'linear');
		[vertices_H, ~, meshes_H] = meshToPet(meshes_H); %2-by-Nv
		meshes_H = meshes_H(1: 3, :);
		free_H = prod(vertices_H .* (1 - vertices_H), 1) ~= 0;
		a_H = zeros(2, 2, size(vertices_H, 2));
		parfor i = 1: size(vertices_H, 2)
			a_H(:, :,  i) = HMM(a, iota, delta, vertices_H(:, i), p_h, meshes_h, free_h, phix_h, phixx_h, T_h);
		end
		error(i) = FEM(H, ux, f, a_H, vertices_H, meshes_H, free_H, phi_H, phix_H, T_H);
		fprintf(" & %.4e", error(i));
	end
	fprintf("\\\\\nrate &");
	for i = 2: size(error, 2)
		rate(i) = log2(error(i-1)/error(i));
		fprintf(" & %.2f", rate(i));
	end
	fprintf("\\\\ \\hline\n");
end


