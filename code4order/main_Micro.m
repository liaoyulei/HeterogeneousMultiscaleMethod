clear;
H = 1/2;
gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes_H = generateMesh(model, 'Hmax', H, 'GeometricOrder', 'linear');
vertices_H = meshToPet(meshes_H); %2-by-Nv
a_H = zeros(2, 2, size(vertices_H, 2));
[phix, phixx, T] = init_fespace();
phix = @(lambda, v) phix(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
phixx = @(lambda, v) phixx(lambda(1), lambda(2), lambda(3), v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
T = @(v) T(v(1, 1), v(1, 2), v(1, 3), v(2, 1), v(2, 2), v(2, 3));
error = zeros(1, 5);
rate = zeros(size(error));
for gamma = [1/4, 1, 4]
	fprintf("%.2f", gamma);
	for i = 1: size(error, 2)
		epsilon = 2^(-5);
		iota = epsilon^gamma;
		N = 2^(i-1);
		delta = 2^(-1);
		h = 2^(-4)/N;
		syms x y real
%%{
		a = matlabFunction([50 + 4*pi^2 * cos(2*pi*x/epsilon), (4*pi^2 + 25/pi) * sin(2*pi*x/epsilon); 0, 4*pi^2-1 + sin(2*pi*x/epsilon)] / (2*pi - cos(2*pi*x/epsilon)));
		a = @(lambda, v) a(v(1, :) * lambda);
		if gamma < 1
			a11 = (50+8*pi^3)/sqrt(4*pi^2-1)-4*pi^2;
		elseif gamma == 1
			a11 = 25/pi;
		else
			a11 = 2*pi^2/((25+4*pi^3)/sqrt(2500-16*pi^4)-1/2);
		end
		a_bar = [a11, 0; 0, sqrt(4*pi^2-1)];
%}
%{
		a = matlabFunction([2+sign(ceil(x/epsilon)+floor(x/epsilon)-2*x/epsilon), 0; 0, 3+2*sign(ceil(y/epsilon)+floor(y/epsilon)-2*y/epsilon)]);
		a = @(lambda, v) a(v(1, :) * lambda, v(2, :) * lambda);
		a_bar = [3/2, 0; 0, 5/3];
%}
		[p_h, meshes_h, free_h] = init_mesh(h); %2-by-Nv; 1-by-Nv, 9-by-Nt
		parfor i = 1: size(vertices_H, 2)
			a_H(:, :,  i) = HMM(a, iota, delta, vertices_H(:, i), p_h, meshes_h, free_h, phix, phixx, T);
		end
		error(i) = sqrt(max(sum((a_H - a_bar).^2, [1, 2]))) / sqrt(max(sum(a_bar.^2, [1, 2])));
		fprintf(" & %.4e", error(i));
	end
	fprintf("\\\\\nrate & ");
	for i = 2: size(error, 2)
		rate(i) = log2(error(i-1)/error(i));
		fprintf(" & %.2f", rate(i));
	end
	fprintf("\\\\ \\hline\n");
end