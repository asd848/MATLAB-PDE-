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
gd = zeros(10, 10*10);
gd(1,:) = 2;
gd(2,:) = 4;

for i=1:10
    for j=1:10
        gd(3:end,10*(i-1)+j) = [window_grid_x(i,j); window_grid_x(i,j+1);...
            window_grid_x(i+1,j+1); window_grid_x(i+1,j); window_grid_y(i,j); ...
            window_grid_y(i,j+1); window_grid_y(i+1,j+1); window_grid_y(i+1,j)];
    end
end
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


    