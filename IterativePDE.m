%% Attempt to code with more than one material
length = 2;
height= length;
width = length;


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
% sq1 = [3;4;0;1;1;0;1;1;0;0];
% sq2 = [3;4;1;2;2;1;1;1;0;0];
% gd = [sq1,sq2];
% sf = "R1 + R2";
% ns = ["R1", "R2"];
g = decsg(gd, sf, ns);
geometryFromEdges(model,g);
% Applying boundary conditions to edges
for i = 1:height
   applyBoundaryCondition(model,'dirichlet','Edge',i,'r',1,'h',1) 
%    applyBoundaryCondition(model,'dirichlet','Edge',height*width*2-i,'r',0,'h',1) 
%     applyBoundaryCondition(model,'neumann','Edge',i,'q',0,'g',0)
%     applyBoundaryCondition(model,'neumann','Edge',i,'q',0,'g',0)
   applyBoundaryCondition(model,'neumann','Edge',height*width*2-i,'q',0,'g',0)
end
for i = 1:width
   applyBoundaryCondition(model,'neumann','Edge',(height*width)-i,'q',0.00,'g',0) 
   applyBoundaryCondition(model,'neumann','Edge',(height*width)*2-length-i,'q',0.00,'g',0) 
end

for i = 1:width*height
   specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0,'face',i) 
end
%Specifying Coefficients


figure(1);
pdegplot(model, 'EdgeLabels', 'on')
% specifyCoefficients(model,'m',0,'d',0,'c',10,'a',0,'f',0,'face',1);
% specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0,'face',2);

% applyBoundaryCondition(model,'dirichlet','Edge',1,'r',1,'h',1);
% % applyBoundaryCondition(model,'neumann','Edge',1,'q',0,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',2,'q',0,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',3,'q',0.05,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',4,'q',0.05,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',5,'q',0.05,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',6,'q',0.05,'g',0);
% applyBoundaryCondition(model, 'neumann', 'Edge', 7, 'q', 0,'g',0);
% 
% 
% 
generateMesh(model);
results = solvepde(model); % for stationary problems

u = results.NodalSolution;
figure(2);
pdeplot(model,'XYData',u,'Mesh','on')
xlabel('x')
ylabel('y')

% g = decsg(gd, sf, ns);
% sf = "";
% gd = ;
% for i = 1:width
%     for j = 1:height
%         pdeObject(i,j) = createpde();
%         sq = [3;4;i-1;i+1;i+1;i-1;j+1;j+1;j-1;j-1;];
%         gd = [sq];
%         ns(i,j) = "R" + num2str(i) + num2str(j);
%         sf = sf + " + " + ns(i,j); 
%         geometryFromEdges(pdeObject(i,j),g);
%     end
% end

%% Just one material
% 
% 
% pdegplot(model, 'FaceLabels', 'on')
% pdegplot(model, 'EdgeLabels', 'on')
% xlim([0 10])
% ylim([0 10])
% 
% applyBoundaryCondition(model,'neumann','Edge',[1,3],...
%                        'q',0.05,'g',0);
% applyBoundaryCondition(model,'neumann','Edge',2,'q',0,'g',0);
% applyBoundaryCondition(model,'dirichlet', 'Edge',4,'h',1,'r',1);
% 
% specifyCoefficients(model,'m',0,'d',0,'c',1,'a',0,'f',0);
% generateMesh(model);
% results = solvepde(model); % for stationary problems
% 
% u = results.NodalSolution;
% pdeplot(model,'XYData',u,'Mesh','on')
% xlabel('x')
% ylabel('y')