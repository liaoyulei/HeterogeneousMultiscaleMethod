function [p, t, free] = init_mesh(h) %p: 2-by-Nv, t: 9-by-Nt
gm = [3; 4; -1/2; 1/2; 1/2; -1/2; -1/2; -1/2; 1/2; 1/2];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes = generateMesh(model, 'Hmax', h, 'GeometricOrder', 'linear');
[p, ~, t] = meshToPet(meshes);
t = [t(1: 3, :); t(1: 3, :) + size(p, 2); t(1: 3, :) + 2*size(p, 2)];
free = prod((1/2 + p) .* (1/2 - p), 1) ~= 0;
%free = [free, free, free]; %clamped
free = [free, (1/2 + p(2, :)) .* (1/2 - p(2, :)), (1/2 + p(1, :)) .* (1/2 - p(1, :))] ~= 0; %simply supported
end