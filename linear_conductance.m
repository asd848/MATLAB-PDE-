%% height and width represent how many squares should be along both the edges
height= 1;
width = 1;
sizeSquare = 1;
dcVoltage=115;

%% S/m according to http://www.mit.edu/~6.777/matprops/ito.htm
conductivity = 1e6; 
numMeasurements = 10; %relevant for sampling the Qj in each square
spacing = 1/numMeasurements;

%% create PDE model
model = createpde(); 

%% Setup the geometry of our complete obejct
gd = (1:10)';
ns = [""];
sf = "";
for i = 1:width
    for j = 1:height
        gd(1:end, end+1) = [3;4;i-1;i;i;i-1;j;j;j-1;j-1;];
        ns(end + 1) = "R" + num2str(i) + num2str(j);
        if sf == ""
            sf = ns(end);
        else
            sf = sf + " + " + ns(end); 
        end
    end
end
gd = gd(1:end, 2:end);
ns = ns(2:end);

g = decsg(gd, sf, ns);
geometryFromEdges(model,g);
figure(1)
pdegplot(g, 'EdgeLabels', 'on')


%% Applying boundary conditions to edges

applyBoundaryCondition(model,'neumann','Edge',4,'q',0,'g',0);
applyBoundaryCondition(model,'neumann','Edge',2,'q',0,'g',0);

applyBoundaryCondition(model,'dirichlet','Edge',1,'r',dcVoltage,'h',1); %bottom boundary
applyBoundaryCondition(model,'dirichlet','Edge',3,'r',-dcVoltage,'h',1); %top boundary

coeffs = @(location,state) conductances(location.x, location.y);

specifyCoefficients(model,'m',0,'d',0,'c',coeffs,'a',0,'f',0,'face',1) ;

generateMesh(model, 'Hmax', 0.05);
solution = solvepde(model); % for stationary problems

u = solution.NodalSolution;
figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')
set(gca,'FontSize',16)
drawnow

%% Finding Voltage Gradient (Electric Field) Code
[xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
mesh = [xData(:) yData(:)]; %% Not the mesh used in solving the PDE !!!!
xData = mesh(1:end,1);
yData = mesh(1:end,2);

[gradx, grady] = evaluateGradient(solution, xData, yData);

gradx = sizeSquare.*gradx;
grady = sizeSquare.*grady;

%% Calculate volumetric joule heating
for i = 1:length(xData)
    for j  = 1:length(yData)
        Qj(i,j) = conductivity* (gradx(i).^2 + grady(j).^2);
    end
end

function output = conductances(x, y)
conductivity = 1e6; 
output = conductivity*(0.1+y);
end

