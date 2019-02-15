%% setting up nxn grid
sizeSquare = 10;
height= sizeSquare;
width = sizeSquare;
load('actual_watt_density_global.mat');
qj = actual_watt_density_global./0.0254^2;

load('actual_thickness_global.mat');
thick = actual_thickness_global.*1e-11;


model = createpde(); 

% Creating the PDE geometry
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
%pdegplot(g)

%% mesh each square individually
qj_xy = cell(10);
for i=1:10
    for j=1:10
        % Creates a new model and meshes a single polygon in order to find (x,y) to
        % evaluate the gradient at. 
        temp_model = createpde();
        temp_geo = decsg(gd(:,10*(i-1)+j));
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model, 'Hmin', 4);
%         pdeplot(temp_model, 'NodeLabels', 'on')
%         pause
        xy_points = mesh.Nodes;
        qj_xy{11-j,i} = xy_points;
    end
end

%% format data for fitting
x_data = [];
y_data = [];
qj_data = [];
thick_data = [];
for i=1:10
    for j=1:10
        x_data = [x_data; qj_xy{i,j}(1,:)'];
        y_data = [y_data; qj_xy{i,j}(2,:)'];
        qj_data = [qj_data; qj(i,j)*ones(length(qj_xy{i,j}(1,:)),1)]; 
        thick_data = [thick_data; thick(i,j)*ones(length(qj_xy{i,j}(1,:)),1)]; 
    end
end

%% fitting data
% f = fit([x_data, y_data], qj_data, 'poly55')
% plot(f, [x_data,y_data], qj_data)
