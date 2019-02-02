% Constants
conductivity = 1E6; %[S/m]
model = createpde(); 

% x and y coordinates of the Global Window
x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

% Assembling the 10x10 polygons into a geometry for MATLAB
gd = zeros(46, 1);
gd(1) = 2;
gd(2) = 22;
gd(3:13) = x_bottom';
gd(14:24) = fliplr(x_top)';
gd(25:35) = y_bottom';
gd(36:46) = fliplr(y_top)';

dl = decsg(gd);
% Plotting the geometry
%figure(1);
%pdegplot(dl);% 'FaceLabels', 'on', 'EdgeLabels', 'on');

gFM = geometryFromEdges(model,dl);

% Applying boundary conditions (no current flux)to left and right edges
% (hard coded for ease)
applyBoundaryCondition(model,'neumann','Edge',22,'q',0,'g',0);
applyBoundaryCondition(model,'neumann','Edge',11,'q',0,'g',0);

% Applying boundary conditions (voltage input )to top and bottom edges
% (hard coded for ease)
applyBoundaryCondition(model,'dirichlet','Edge',(12:1:21),'r',dcVoltage,'h',1);
applyBoundaryCondition(model,'dirichlet','Edge',(1:1:10),'r',-dcVoltage,'h',1);

coeffs = @(location, state) conductances(location.x, location.y)*(x_bottom(11)-x_bottom(1))/((x_bottom(11)-x_bottom(1))*(y_top(1)-y_bottom(1)));

%c = 1E-7*(x_bottom(11)-x_bottom(1))/((x_bottom(11)-x_bottom(1))*(y_top(1)-y_bottom(1)));

specifyCoefficients(model,'m',0,'d',0,'c',c,'a',0,'f',0,'face',1);
generateMesh(model, 'Hmax', 0.05);
solution = solvepde(model); % for steady state problems

u = solution.NodalSolution;

% Plotting the voltage drop across the geometry
figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')

function out = conductances(x,y)
    conductivity = 1e6;
    out = conductivity*y^2;
end
