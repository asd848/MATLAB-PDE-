run('global_window_comsol_iteration.m');

%%
delta_matrix = load('current_delta.csv');
x = delta_matrix(:,1);
y = delta_matrix(:,2);
delta = delta_matrix(:,3);
pause
surface_fit = fit([x y], delta, 'linearinterp');

delta_matrix = zeros(10,10);

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

%% Assembling the 10x10 polygons into a geometry for MATLAB
gd = zeros(10, 10*10);
gd(1,:) = 2;
gd(2,:) = 4;

for i=1:10
    for j=1:10
        gd(3:end,10*(i-1)+j) = [x(i,j); x(i,j+1); x(i+1,j+1); x(i+1,j); y(i,j); y(i,j+1); y(i+1,j+1); y(i+1,j)];
    end
end

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
        
        fit_data = surface_fit(x_data,y_data);
        fit_data = fit_data(~isnan(fit_data));
        
        delta_matrix(i,j) = mean(fit_data);
    end
end

image(delta_matrix)

