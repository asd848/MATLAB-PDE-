function [outThick] = jouleHeater(sizeSquare,dcVoltage, thickness, qjDes)
height= sizeSquare;
width = sizeSquare;
conductivity = 1e6; %% made up number
numMeasurements = 10; %relevant for sampling the Qj in each square
spacing = 1/numMeasurements;

model = createpde(); 

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

% Applying boundary conditions to edges
for i = 1:height
    %Left and Right Edges
   applyBoundaryCondition(model,'neumann','Edge',i,'q',0,'g',0);
   applyBoundaryCondition(model,'neumann','Edge',height*width*2-i,'q',0,'g',0);
end
for i = 1:width
    %Top and Bottom Edges (where we set our voltage)
   applyBoundaryCondition(model,'dirichlet','Edge',(height*width)-i+1,'r',-dcVoltage,'h',1); %bottom boundary
   applyBoundaryCondition(model,'dirichlet','Edge',(height*width)*2-sizeSquare-i+1,'r',dcVoltage,'h',1); %top boundary
   
%    applyBoundaryCondition(model,'dirichlet','Edge',3,'r',-dcVoltage,'h',1); 
%    applyBoundaryCondition(model,'dirichlet','Edge',1,'r',dcVoltage,'h',1); 

end


%In the for loop below, 'c' is the electric conductance of the material.
%Eventually, it should be referencing a matrix as it iterates, but for now,
%it is set to 1.
conductance = zeros(sizeSquare, sizeSquare);

for i = 1:width
    for j = 1:height
       conductance(i,j) = conductivity*thickness(i,j)*width/sizeSquare;
       specifyCoefficients(model,'m',0,'d',0,'c',conductance(i,j),'a',0,'f',0,'face',(i-1)*height+j) ;
    end
end


% figure(4);
% pdegplot(model, 'EdgeLabels', 'on')

generateMesh(model);
solution = solvepde(model); % for stationary problems

u = solution.NodalSolution;
figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')
drawnow

% finalModel = model;


%% Finding Voltage Gradient (Electric Field) Code
[xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
mesh = [xData(:) yData(:)]; %% Not the mesh used in solving the PDE !!!!
xData = mesh(1:end,1);
yData = mesh(1:end,2);

% uintrp = interpolateSolution(results,xData, yData);

[gradx, grady] = evaluateGradient(solution, xData, yData);

gradx = sizeSquare.*gradx;
grady = sizeSquare.*grady;

xData = unique(xData);
yData = unique(yData);

Qj= zeros(sizeSquare*numMeasurements,sizeSquare*numMeasurements);

for i = 1:length(xData)
    for j  = 1:length(yData)
        Qj(i,j) = gradx(i).^2 + grady(j).^2;
    end
end

% Calculate average Qj for each square

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

%% Getting output thickness
outThick = zeros(sizeSquare, sizeSquare);

outThick  = qjDes./avgQj;

end

