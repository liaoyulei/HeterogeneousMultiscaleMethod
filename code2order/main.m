clear;
epsilon = 1/100;
delta = 10 * epsilon;
a = @(x) [5+3/2*cos(2*pi*x(1)/epsilon), 2+sin(2*pi*x(1)/epsilon); 3+2*cos(2*pi*x(1)/epsilon), 6+2*sin(2*pi*x(1)/epsilon)];
a_0 = [sqrt(91)/2, 2; 2/3*sqrt(91)-11/3, 6];
H = 1;
h = 1/128;
gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes_H = generateMesh(model, 'Hmax', H, 'GeometricOrder', 'linear');
vertices_H = meshToPet(meshes_H); %2-by-Nv
[vertices_h, meshes_h] = init_mesh(h); %2-by-Nv; 3-by-Nt
not_bdr = prod((vertices_h + 1/2) .* (1/2 - vertices_h), 1);
free_h = not_bdr ~= 0; %1-by-Nv
[phix, T] = init_fespace(); %3-by-2, 3-by-1, x_1, x_2, x_3, y_1, y_2, y_3
a_H = zeros(2, 2, size(vertices_H, 2));
parfor i = 1: size(vertices_H, 2)
	a_H(:, :, i) = HMM(a, delta * vertices_h + vertices_H(:, i), meshes_h, free_h, phix, T) / delta^2;
end