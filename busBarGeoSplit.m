function [outThick, qj] = busBarGeoSplit(dcVoltage, thickness, qjDes)
close all;
conductivity = 1e6; %% made up number
model = createpde(); 

x_bottom = [0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = [0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = [0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = [20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

x_mid = [0:30.17/10:30.17];
y_mid = [8.958: (31.21-8.958)/10: 31.21];

x_top = fliplr(x_top);
y_top = fliplr(y_top);

gd_ll = [2, 12, x_bottom(1:6), fliplr(x_mid(1:6)), y_bottom(1:6),fliplr(y_mid(1:6))]';
gd_lr = [2, 12, x_bottom(6:11), fliplr(x_mid(6:11)), y_bottom(6:11), fliplr(y_mid(6:11))]';

gd_tl = [2, 12, x_mid(1:6), x_top(6:11), y_mid(1:6), y_top(6:11)]';
gd_tr = [2, 12, x_mid(6:11), x_top(1:6), y_mid(6:11), y_top(1:6)]';

gd = [gd_ll gd_lr gd_tl gd_tr];
ns = ["Rll" "Rlr" "Rtl" "Rtr"];
sf = "Rll + Rlr + Rtl + Rtr";

dl = decsg(gd, sf, ns);
figure;
pdegplot(dl)

gFM = geometryFromEdges(model,dl);

% pdegplot(dl, 'EdgeLabels','on','VertexLabels','on');

applyBoundaryCondition(model,'neumann','Edge',[11 17 24 25],'q',0,'g',0);
applyBoundaryCondition(model,'dirichlet','Edge',(12:1:16),'r',dcVoltage,'h',1);
applyBoundaryCondition(model,'dirichlet','Edge',(18:1:22),'r',dcVoltage,'h',1);
applyBoundaryCondition(model,'dirichlet','Edge',(1:1:10),'r',-dcVoltage,'h',1);

conductance = ones(2,2);

conductance(1,1) = conductivity*thickness(1,1)*15.09;
conductance(1,2) = conductivity*thickness(1,2)*(30.17-15.09);
conductance(2,1) = conductivity*thickness(2,1)*18.28;
conductance(2,2) = conductivity*thickness(2,2)*(35.46-18.28);

specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,1),'a',0,'f',0,'face',2);
specifyCoefficients(model,'m',0,'d',0,'c',conductance(1,2),'a',0,'f',0,'face',4);
specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,1),'a',0,'f',0,'face',1);
specifyCoefficients(model,'m',0,'d',0,'c',conductance(2,2),'a',0,'f',0,'face',3);

% specifyCoefficients(model,'m',0,'d',0,'c',1e6*35,'a',0,'f',0);
gM = generateMesh(model);
solution = solvepde(model); % for stationary problems

u = solution.NodalSolution;
figure;
pdeplot(model,'XYData',u)
title('Voltage Map of the Window')
xlabel('x')
ylabel('y')

%% temp mesh per square
temp_model_ll = createpde();

temp_geo = decsg(gd_ll);
geometryFromEdges(temp_model_ll,temp_geo);
mesh = generateMesh(temp_model_ll);

xData_ll = mesh.Nodes(1,:);
yData_ll = mesh.Nodes(2,:);

temp_model_lr = createpde();
temp_geo = decsg(gd_lr);
geometryFromEdges(temp_model_lr, temp_geo);
mesh = generateMesh(temp_model_lr);

xData_lr = mesh.Nodes(1,:);
yData_lr = mesh.Nodes(2,:);

temp_model_tl = createpde();
temp_geo = decsg(gd_tl);
geometryFromEdges(temp_model_tl, temp_geo);
mesh = generateMesh(temp_model_tl);

xData_tl = mesh.Nodes(1,:);
yData_tl = mesh.Nodes(2,:);

temp_model_tr = createpde();
temp_geo = decsg(gd_tr);
geometryFromEdges(temp_model_tr, temp_geo);
mesh = generateMesh(temp_model_tr);

xData_tr = mesh.Nodes(1,:);
yData_tr = mesh.Nodes(2,:);

%% finding the gradient

%lower left
[gradx_ll, grady_ll] = evaluateGradient(solution, xData_ll, yData_ll);

for i = 1:length(xData_ll)
    for j  = 1:length(yData_ll)
        Qj_ll(i,j) = conductivity*(gradx_ll(i).^2 + grady_ll(j).^2);
    end
end

%lower right
[gradx_lr, grady_lr] = evaluateGradient(solution, xData_lr, yData_lr);

for i = 1:length(xData_lr)
    for j  = 1:length(yData_lr)
        Qj_lr(i,j) = conductivity*(gradx_lr(i).^2 + grady_lr(j).^2);
    end
end

%top left
[gradx_tl, grady_tl] = evaluateGradient(solution, xData_tl, yData_tl);

for i = 1:length(xData_tl)
    for j  = 1:length(yData_tl)
        Qj_tl(i,j) = conductivity*(gradx_tl(i).^2 + grady_tl(j).^2);
    end
end

%top right
[gradx_tr, grady_tr] = evaluateGradient(solution, xData_tr, yData_tr);

for i = 1:length(xData_tr)
    for j  = 1:length(yData_tr)
        Qj_tr(i,j) = conductivity*(gradx_tr(i).^2 + grady_tr(j).^2);
    end
end

% Calculate average Qj for each square

avgQj = [mean(mean(Qj_tl)) mean(mean(Qj_tr)); mean(mean(Qj_ll)) mean(mean(Qj_lr))];

%% Getting output thickness
outThick  = qjDes./avgQj;
qj = outThick.*avgQj;
end