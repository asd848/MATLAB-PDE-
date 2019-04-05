function [qj] = resisty_joule_heating(sizeSquare,dcVoltage, thickness)
%INPUTS:
%sizeSquare:
%dcVoltage: voltage input (V)
%thickness: inital guess for ITO thickness
%qjDes: desired heating in watt/sq
%OUTPUTS:
%outThick: output thickness per square
%qj: heating per square

%% height and width represent how many squares should be along both the edges
height= sizeSquare;
width = sizeSquare;

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
pdegplot(g,'FaceLabels','on')

%% Applying boundary conditions to edges
for i = 1:height
   %Left and Right Edges
   applyBoundaryCondition(model,'neumann','Edge',i,'q',0,'g',0);
   applyBoundaryCondition(model,'neumann','Edge',height*width*2-i,'q',0,'g',0);
end
for i = 1:width
   %Top and Bottom Edges (where we set our voltage)
   applyBoundaryCondition(model,'dirichlet','Edge',(height*width)-i+1,'r',-dcVoltage,'h',1); %bottom boundary
   applyBoundaryCondition(model,'dirichlet','Edge',(height*width)*2-sizeSquare-i+1,'r',dcVoltage,'h',1); %top boundary
end


%Interpolating thickness to make it a continuous function
%If it isn't continuous then it continuity boundary conditions between
%square interfaces won't be able to work.

%% Decompose thickness matrix into x-y positions with a z value
xLength = size(thickness,2);
yLength = size(thickness,1);
xSpacing = linspace(0,xLength,xLength);
ySpacing = linspace(yLength,0,yLength);
xPos = [];
yPos = [];
zThickness = [];
for k = 1:xLength
    for m = 1:yLength
        xPos = [xPos xSpacing(k)];
        yPos = [yPos ySpacing(m)];
        zThickness= [zThickness thickness(m,k)];
    end
end
xPos = xPos';
yPos = fliplr(yPos);
yPos = yPos';

zThickness = zThickness';
%surfaceFit is a function z(x,y), which takes in x and y position and
%outputs thickness. Fit can be adjusted with last parameter.
% WARNING!!! If program crashes it may be because zThickness did not have
% enough points for the fit to work. Try changing 'poly23' to 'poly12'.
disp('Surface fit');
surfaceFit = fit([xPos,yPos],zThickness,'linearinterp')
%surfaceFit = fit([xPos,yPos],zThickness,'poly33')


%% Apply coefficients to PDE
%In the for loop below, 'c' is the electric conductance of the material.
conductance = @(location,state) surfaceFit(location.x,location.y)*conductivity*width./sizeSquare; %sizeSquare is to scale our thing down to 1
% conductance = thickness.*conductivity.*width./sizeSquare;
specifyCoefficients(model,'m',0,'d',0,'c',conductance,'a',0,'f',0,'face',[1:100]);
% for i = 1:width
%     for j = 1:height
%        specifyCoefficients(model,'m',0,'d',0,'c',conductance(10-(j-1),i),'a',0,'f',0,'face',(i-1)*height+j) ;
%     end
% end
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,10),'a',0,'f',0,'face',99);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,10),'a',0,'f',0,'face',100);
% 
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,1),'a',0,'f',0,'face',5);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,2),'a',0,'f',0,'face',10);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,3),'a',0,'f',0,'face',15);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,4),'a',0,'f',0,'face',20);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,5),'a',0,'f',0,'face',24);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,1),'a',0,'f',0,'face',4);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,2),'a',0,'f',0,'face',9);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,3),'a',0,'f',0,'face',14);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,4),'a',0,'f',0,'face',19);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,5),'a',0,'f',0,'face',25);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(3,1),'a',0,'f',0,'face',3);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(3,2),'a',0,'f',0,'face',8);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(3,3),'a',0,'f',0,'face',13);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(3,4),'a',0,'f',0,'face',18);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(3,5),'a',0,'f',0,'face',23);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(4,1),'a',0,'f',0,'face',2);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(4,2),'a',0,'f',0,'face',7);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(4,3),'a',0,'f',0,'face',12);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(4,4),'a',0,'f',0,'face',17);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(4,5),'a',0,'f',0,'face',22);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(5,1),'a',0,'f',0,'face',1);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(5,2),'a',0,'f',0,'face',6);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(5,3),'a',0,'f',0,'face',11);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(5,4),'a',0,'f',0,'face',16);
% specifyCoefficients(model,'m',0,'d',0,'c',conductance(5,5),'a',0,'f',0,'face',21);

%% Get voltage solution from mesh of geometry
generateMesh(model,'hmax',0.3);
solution = solvepde(model); % for stationary problems

u = solution.NodalSolution;
% figure(2);
% pdeplot(model,'XYData',u,'Mesh','on')
% title('Voltage Map of the Window')
% xlabel('x')
% ylabel('y')
% set(gca,'FontSize',16)
% drawnow

%% Finding Voltage Gradient (Electric Field) Code
[xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
mesh = [xData(:) yData(:)]; %% Not the mesh used in solving the PDE !!!!
xData = mesh(1:end,1);
yData = mesh(1:end,2);

[gradx, grady] = evaluateGradient(solution, xData, yData);

gradx = sizeSquare.*gradx;
grady = sizeSquare.*grady;

gradv = sqrt(gradx.^2+grady.^2);

%% Calculate volumetric joule heating
for i = 1:length(xData)
    for j  = 1:length(yData)
        Qj(i,j) = conductivity* (gradx(i).^2 + grady(j).^2);
    end
end

%% Calculate average Qj for each square
avgQj = zeros(sizeSquare,sizeSquare);
for i = 1:sizeSquare
    for j = 1:sizeSquare
        xStartIndex = numMeasurements * (i -1) + 1;
        xEndIndex = numMeasurements * i;
        yStartIndex = numMeasurements * (j -1) + j;
        yEndIndex = numMeasurements * j;

        avgQj(i,j) = mean(mean(Qj(xStartIndex:xEndIndex,yStartIndex:yEndIndex)));
    end
end

%% Getting output qj
qj = avgQj.*thickness;
end

