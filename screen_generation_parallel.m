run('global_window_comsol_iteration.m');

%%
delta_matrix = load('current_delta.csv');
x = delta_matrix(:,1);
y = delta_matrix(:,2);
delta = delta_matrix(:,3);
surface_fit = fit([x y], delta, 'linearinterp');

%%
% x and y coordinates of the Global Window
x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

top_bus_interp = fit(x_top', y_top', 'linearinterp');
bottom_bus_interp = fit(x_bottom', y_bottom', 'linearinterp');

rightm = (y_top(11)-y_bottom(11))/(x_top(11)-x_bottom(11));

x_bot_points = linspace(x_bottom(1), x_bottom(11), 11);
y_bot_points = bottom_bus_interp(x_bot_points);

x_top_points = zeros(1,11);
y_top_points = zeros(1,11);

window_grid_x = zeros(11,11);
window_grid_y = zeros(11,11);

for i=11:-1:1
    parallel_lines = @(x) (rightm*(x-x_bot_points(i))+y_bot_points(i));
    fun = @(x) top_bus_interp(x) - parallel_lines(x);
    x_zero = fzero(fun, x_top(i));
    y_zero = top_bus_interp(x_zero);
    x_top_points(i) = x_zero;
    y_top_points(i) = y_zero;
    x_query = linspace(x_bot_points(i),x_top_points(i),11);
    y_query = parallel_lines(x_query);
    window_grid_x(:,i) = x_query;
    window_grid_y(:,i) = y_query;
end

%% Assembling the 10x10 polygons into a geometry for MATLAB
gd = zeros(12, 10*10);

index = 1;
for i=1:10
    for j=1:10
        if (window_grid_x(i,j)>= 0 && window_grid_x(i,j+1)>=0 && window_grid_x(i+1,j) >=0 && window_grid_x(i+1,j+1)>=0)
            gd(1,index) = 2;
            gd(2,index) = 4;
            gd(3:10,index) = [window_grid_x(i,j); window_grid_x(i,j+1);...
                window_grid_x(i+1,j+1); window_grid_x(i+1,j); window_grid_y(i,j); ...
                window_grid_y(i,j+1); window_grid_y(i+1,j+1); window_grid_y(i+1,j)];
            index = index+1;
            disp('all positive');
            
        elseif(window_grid_x(i,j)>= 0 && window_grid_x(i,j+1)>= 0 && window_grid_x(i+1,j) <0 && window_grid_x(i+1,j+1)>=0)               
            m1 = (window_grid_y(i+1,j+1)-window_grid_y(i+1,j))/(window_grid_x(i+1,j+1)-window_grid_x(i+1,j));
            y1 = m1*-(window_grid_x(i+1,j+1)) + window_grid_y(i+1, j+1);
            m2 = (window_grid_y(i+1,j)-window_grid_y(i,j))/(window_grid_x(i+1,j)-window_grid_x(i,j));
            y2 = m2*-(window_grid_x(i+1,j)) + window_grid_y(i+1, j);
              
            if(window_grid_x(i,j)==0 && window_grid_y(i,j) ==0)
                gd(1,index) = 2;
                gd(2,index) = 4; 
                gd(3:10,index) = [window_grid_x(i,j); window_grid_x(i,j+1); window_grid_x(i+1,j+1); 0;...
                    window_grid_y(i,j); window_grid_y(i,j+1); window_grid_y(i+1,j+1); y1]
            else
                gd(1,index) = 2;
                gd(2,index) = 5; 
                gd(3:end,index) = [window_grid_x(i,j); window_grid_x(i,j+1); window_grid_x(i+1,j+1); 0; 0;...
                    window_grid_y(i,j); window_grid_y(i,j+1); window_grid_y(i+1,j+1); y1; y2]
            end
            
            index = index+1;
            disp('one negative');
            
        elseif(window_grid_x(i,j)< 0 && window_grid_x(i,j+1)>=0 && window_grid_x(i+1,j) <0 && window_grid_x(i+1,j+1)>=0)
            gd(1,index) = 2;
            gd(2,index) = 4;
            m1 = (window_grid_y(i+1,j+1)-window_grid_y(i+1,j))/(window_grid_x(i+1,j+1)-window_grid_x(i+1,j));
            y1 = m1*-(window_grid_x(i+1,j+1)) + window_grid_y(i+1, j+1);
            m2 = (window_grid_y(i,j+1)-window_grid_y(i,j))/(window_grid_x(i,j+1)-window_grid_x(i,j));
            y2 = m2*-(window_grid_x(i,j+1)) + window_grid_y(i+1,j);
                
            gd(3:10,index) = [0; window_grid_x(i,j+1); window_grid_x(i+1,j+1); 0;...
                y2; window_grid_y(i,j+1); window_grid_y(i+1,j+1); y1];
            index = index+1;
            disp('half negative');
            
        elseif(window_grid_x(i,j)< 0 && window_grid_x(i,j+1)>=0 && window_grid_x(i+1,j) <0 && window_grid_x(i+1,j+1)<0)
            gd(1,index) = 2;
            gd(2,index) = 3;
            m1 = (window_grid_y(i+1,j+1)-window_grid_y(i,j+1))/(window_grid_x(i+1,j+1)-window_grid_x(i,j+1));
            y1 = m1*-(window_grid_x(i+1,j+1)) + window_grid_y(i+1, j+1);
            m2 = (window_grid_y(i,j+1)-window_grid_y(i,j))/(window_grid_x(i,j+1)-window_grid_x(i,j));
            y2 = m2*-(window_grid_x(i,j+1)) + window_grid_y(i, j+1);
            
            gd(3:8,index) = [0; 0; window_grid_x(i,j+1);...
                y2; window_grid_y(i,j+1);y1] 
            index = index+1;
            disp('one positive');
        else
            disp('all negative');
        end
    end
end

end_index = find(gd(1,:)==0);
gd = gd(1:end, 1:(end_index(1)-1));


% R1 = [gd; zeros(4,100)];
% C1 = [2,6,0, fliplr(x_top_points(1:4)), 0, 0.5309, fliplr(y_top_points(1:4)),0]';
% gd = [R1, C1];
% 
% sf = 'R1-C1';
% ns = char('R1', 'C1')';
model = createpde(); 
% dl = decsg(gd,sf,ns);
% Plotting the geometry
dl = decsg(gd);
figure(1);
pdegplot(dl) %'FaceLabels', 'on', 'EdgeLabels', 'on')

%%
window_geom = csvread('perf_geom.csv');
delta_matrix = zeros(10,10);
idx = 1;
figure(5);
hold on;
for i=1:3
    for j=1:10
         temp_model = createpde();
        temp_geo = decsg(window_geom(:,idx));
        pdegplot(temp_geo)
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(11-i,j) = mean(fit_data);
        idx = idx+1;
        pause;
    end
end

for i=4:5
    for j=2:10
        temp_model = createpde();
        temp_geo = decsg(window_geom(:,idx));
        pdegplot(temp_geo)
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(11-i,j) = mean(fit_data);
        idx = idx+1;
        pause;

    end
end

for i=6:8
    for j=3:10
         temp_model = createpde();
        temp_geo = decsg(window_geom(:,idx));
        pdegplot(temp_geo)
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(11-i,j) = mean(fit_data);
        idx=idx+1;
        pause;

    end
end

for i=9:10
    for j=4:10
         temp_model = createpde();
        temp_geo = decsg(window_geom(:,idx));
        pdegplot(temp_geo)
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(11-i,j) = mean(fit_data);
        idx=idx+1;
        pause;
    end
end

%%
figure(6)
 temp_model = createpde();
        temp_geo = decsg(window_geom(:,39));
        pdegplot(temp_geo)
        geometryFromEdges(temp_model,temp_geo);
        mesh = generateMesh(temp_model);
        generateMesh(temp_model);
        x_data = mesh.Nodes(1,:);
        y_data = mesh.Nodes(2,:);
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(11-i,j) = mean(fit_data);
        idx = idx+1;
