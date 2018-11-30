function [outThick, qj] = busBarGeoSplitFull(dcVoltage, thickness, qjDes)
% busBarGeoSplitFull solves the electrostatics PDE and calculates thickness
% for a desired surface watt density on the actual geometry of the Global
% Window
% INPUTS:
% dcVoltage: voltage input (V)
% thickness: inital guess for ITO thickness
% qjDes: desired heating in watt/sq
% OUTPUTS:
% outThick: output thickness per square
% qj: heating per square

% Constants
conductivity = 6.41e5; %[S/m]
model = createpde(); 

% x and y coordinates of the Global Window
x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

% Creating 10x10 polygons to cover the surface area of the Global window
x = zeros(11,11);
y = zeros(11,11);

num_slope = y_top - y_bottom;
den_slope = x_top-x_bottom;

slope = num_slope./den_slope;

for i=1:11
    if slope(i) == inf
        slope(i) = 0;
    else
    end
end

b = y_bottom - slope.*x_bottom;

x_diff = x_bottom-x_top;

for i=1:11
    if x_diff(i) ~= 0
        x(:,i) = x_top(i):x_diff(i)/10:x_bottom(i)';
    else
        x(2:end-1,i) = x_top(i)*ones(9,1);
    end
    
    if slope~=0
        y(:,i) = slope(i)*x(:,i) + b(i);
    else
        y(:,i) = fliplr(y_bottom(i):(y_top(i)-y_bottom(i))/10:y_top(i))';
    end

end

% Assembling the 10x10 polygons into a geometry for MATLAB
gd = zeros(10, 10*10);
gd(1,:) = 2;
gd(2,:) = 4;

for i=1:10
    for j=1:10
        gd(3:end,10*(i-1)+j) = [x(i,j); x(i,j+1); x(i+1,j+1); x(i+1,j); y(i,j); y(i,j+1); y(i+1,j+1); y(i+1,j)];
    end
end

dl = decsg(gd);
% Plotting the geometry
figure(1);
pdegplot(dl) %'FaceLabels', 'on', 'EdgeLabels', 'on')

gFM = geometryFromEdges(model,dl);

% Applying boundary conditions (no current flux)to left and right edges
% (hard coded for ease)
applyBoundaryCondition(model,'neumann','Edge',(102:1:111),'q',0,'g',0);
applyBoundaryCondition(model,'neumann','Edge',(201:1:210),'q',0,'g',0);

% Applying boundary conditions (voltage input )to top and bottom edges
% (hard coded for ease)
applyBoundaryCondition(model,'dirichlet','Edge',(1:1:10),'r',dcVoltage,'h',1);
applyBoundaryCondition(model,'dirichlet','Edge',(11:1:20),'r',-dcVoltage,'h',1);

% Applying the conductance to each polygon
conductance = ones(10,10);

for i=1:10
    for  j=1:10
%         conductance(i,j) = conductivity*thickness(i,j)*sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
        conductance(i,j) = conductivity*sqrt((x(i,j)-x(i,j+1))^2+(y(i,j)-y(i,j+1))^2);
    end
end

for i=1:10
    for j=1:10
        if i==1
            c = conductance(i, j);
        elseif i==2
            c = conductance(end,j);
        else
            c = conductance(i-1, j);
        end
        specifyCoefficients(model,'m',0,'d',0,'c',c,'a',0,'f',0,'face',10*(i-1)+j);
    end
end

% Creating the mesh, 'Hmax' refines the mesh
generateMesh(model, 'Hmax', 1.2);
solution = solvepde(model); % for steady state problems

u = solution.NodalSolution;

% Plotting the voltage drop across the geometry
figure(2);
pdeplot(model,'XYData',u, 'mesh', 'on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')
drawnow

%% Evaluate the gradient
Qj = zeros(10,10);
for i=1:10
    for j=1:10
        % Creates a new model and meshes a single polygon in order to find (x,y) to
        % evaluate the gradient at. 
        temp_model = createpde();
        temp_geo = decsg(gd(:,10*(i-1)+j));
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        % Evaluating the gradient
        [gradx, grady] = evaluateGradient(solution, x_data, y_data);
        
        % Calculating average volumetric joule heating (Qj)
        Qj_temp = zeros(length(x_data), length(y_data));
        for k=1:length(x_data)
            for m=1:length(y_data)
                Qj_temp(k,m) = conductivity*(gradx(k).^2 + grady(m).^2);
            end
        end
        Qj(i,j) = mean(mean(Qj_temp));
    end
end

%% Getting output thickness
outThick  = qjDes./Qj;
qj = thickness.*Qj;
end
