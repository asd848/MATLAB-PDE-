function [finalModel, solution] = jouleHeater(sizeSquare,dcVoltage)
height= sizeSquare;
width = sizeSquare;


model = createpde(); 

gd = (1:10)';
%ns = ['';''];
ns = char('R11', 'R12', 'R21', 'R22')';
sf = 'R11+R12+R21+R22';
for i = 1:width
    for j = 1:height
        gd(1:end, end+1) = [3;4;i-1;i;i;i-1;j;j;j-1;j-1;];
%         %ns = [ns;strcat('R' , num2str(i) , num2str(j))]
%         if strcmp(sf, '')
%             sf = ns(end)
%         else
%             sf = strcat(sf , ' + ' , ns(end))
%         end
    end
end
gd = [3 3 3 3; 4 4 4 4; 0 0 1 1; 1 1 2 2; 1 1 2 2; 0 0 1 1; 1 2 1 2; 1 2 1 2; 0 1 0 1; 0 1 0 1]
% sf
disp(sf)
disp(ns)
% gd = gd(1:end, 2:end);
% % ns = ns(2:end);

%ns = [strcat(ns(1,1), ns(1,2), ns(1,3)) ; strcat(ns(2,1), ns(2,2), ns(2,3)); strcat(ns(3,1), ns(3,2), ns(3,3)); ...
            %strcat(ns(4,1), ns(4,2), ns(4,3))]

% % disp(ns)

g = decsg(gd, sf, ns)
geometryFromEdges(model,g);
% Applying boundary conditions to edges
for i = 1:height
    %Left and Right Edges
   applyBoundaryCondition(model,'edge',i,'q',0,'g',0)
   applyBoundaryCondition(model,'edge',height*width*2-i,'q',0,'g',0)
end
for i = 1:width
    %Top and Bottom Edges (where we set our voltage)
   applyBoundaryCondition(model,'edge',(height*width)-i+1,'r',-dcVoltage,'h',1) 
   applyBoundaryCondition(model,'edge',(height*width)*2-sizeSquare-i+1,'r',dcVoltage,'h',1) 
end


%In the for loop below, 'c' is the electric conductivity of the material.
%Eventually, it should be referencing a matrix as it iterates, but for now,
%it is set to 1.
sigma = [1 1 1 1];
for i = 1:width*height
   specifyCoefficients(model,'m',0,'d',0,'c',sigma(i),'a',0,'f',0,'face',i) 
end


figure(4);
pdegplot(model, 'FaceLabels', 'on')

generateMesh(model);
solution = solvepde(model); % for stationary problems

u = solution.NodalSolution;
figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')

finalModel = model;
end

