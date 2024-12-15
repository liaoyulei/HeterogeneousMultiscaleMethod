function [p, t, trans] = init_mesh(h) %p: 2-by-Nv, t: 9-by-Nt
gm = [3; 4; -1/2; 1/2; 1/2; -1/2; -1/2; -1/2; 1/2; 1/2];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes = generateMesh(model, 'Hmax', h, 'GeometricOrder', 'linear');
[p, ~, t] = meshToPet(meshes);
t = [t(1: 3, :); t(1: 3, :) + size(p, 2); t(1: 3, :) + 2*size(p, 2)]; %3\x N
trans = 1: size(p, 2);
trans(2: 4) = 1;
trans(2/h + 3: 3/h + 1) = 1/h + 3: -1: 5;
trans(3/h + 2: 4/h) = 2/h + 2: -1: 1/h + 4;
trans = [trans, trans + size(p, 2), trans + 2*size(p, 2)];
end