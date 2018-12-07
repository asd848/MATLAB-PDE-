model = createpde(); 

x_bottom = [0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = [0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = [0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = [20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

x_top = fliplr(x_top);
y_top = fliplr(y_top);

gd = [2, 22, x_bottom, x_top, y_bottom, y_top]';
dl = decsg(gd);
geom = geometryFromEdges(model,dl);

pdegplot(dl, 'EdgeLabels','on','VertexLabels','on', 'FaceLabels','on');
% 
% dcVoltage = 115;
% applyBoundaryCondition(model,'neumann','Edge',[11, 22],'q',0,'g',0);
% applyBoundaryCondition(model,'dirichlet','Edge',(12:1:21),'r',dcVoltage,'h',1);
% applyBoundaryCondition(model,'dirichlet','Edge',(1:1:10),'r',-dcVoltage,'h',1);
% 
% specifyCoefficients(model,'m',0,'d',0,'c',1e6*35,'a',0,'f',0);
% gM = generateMesh(model);
% solution = solvepde(model); % for stationary problems
% 
% u = solution.NodalSolution;
% figure(2);
% pdeplot(model,'XYData',u,'Mesh','on', 'ElementLabels','on')
% title('Voltage Map of the Window')
% xlabel('x')
% ylabel('y')