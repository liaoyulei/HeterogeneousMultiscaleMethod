function [p, t] = init_mesh(h) %p: 2-by-Nv, t: 3-by-Nt
gm = [3; 4; -1/2; 1/2; 1/2; -1/2; -1/2; -1/2; 1/2; 1/2];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes = generateMesh(model, 'Hmax', h, 'GeometricOrder', 'linear');
[p, ~, t] = meshToPet(meshes);
t = t(1: 3, :);
end