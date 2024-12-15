gm = [3; 4; 0; 1; 1; 0; 0; 0; 1; 1];
model = createpde(1);
geometryFromEdges(model, decsg(gm));
meshes = generateMesh(model, 'Hmax', 1/2, 'GeometricOrder', 'linear');
[p, e, t] = meshToPet(meshes);
figure
hold on
axis equal
axis off
box off
pdemesh(p, e, t, NodeLabels="on"); 
%print -depsc2 meshes.eps
for i = [2,9,11,4]
	fprintf("(%.2f,%.2f)--", 3*p(1, i)+7, 3*p(2, i)+1);
end
print -dpng meshes.png