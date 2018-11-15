model = createpde();

voltage = 115;
gd = [3;4;0;1;1;0;1;1;0;0];
ns = 'R';
sf = 'R';

g = decsg(gd,sf,ns);
geometryFromEdges(model,g);
pdegplot(model, 'EdgeLabels', 'on');

applyBoundaryCondition(model, 'edge', 1, 'r', voltage);
applyBoundaryCondition(model, 'edge', 3, 'r', -voltage);
applyBoundaryCondition(model, 'edge', 2, 'g', 0);
applyBoundaryCondition(model, 'edge', 4, 'g', 0);

specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0,'face',1) 
generateMesh(model);
solution = solvepde(model);
u = solution.NodalSolution;

figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')