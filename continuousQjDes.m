% qj_des = 1.4/0.0254^2*ones(10,10);

qj_des = randi(10,10)

x_bottom = 0.0254*[0 3.622 7.411 11.107 14.724 18.279 21.772 25.214 28.611 32.037 35.463];
y_bottom = 0.0254*[0 0.815 3.101 5.532 8.081 10.715 13.431 16.211 19.046 22.421 25.795];

x_top = 0.0254*[0 2.712 5.378 7.994 10.566 13.09 15.567 18.005 20.404 22.764 25.09];
y_top = 0.0254*[20.902 22.085 23.37 24.754 26.217 27.76 29.38 31.057 32.79 34.576 36.404];

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

x_cen = [];
y_cen = [];

for i=1:10
    for j=1:10
%         x_cen(10*(i-1)+j,1) = x(i,j) - 0.5*(x(i,j)-x(i,j+1));
%         y_cen(10*(i-1)+j,1) = y(i,j) - 0.5*(y(i,j)-y(i+1,j));
        x_cen(10*(i-1)+j,1) = 0.25*(x(i,j) + x(i+1,j) + x(i+1,j+1) + x(i,j+1));
        y_cen(10*(i-1)+j,1) = 0.25*(y(i,j) + y(i+1,j) + y(i+1,j+1) + y(i,j+1));
    end
end
figure;
for i = 1:11
    hold on;
    plot(x(:,i), y(:,i))
    plot(x(i,:), y(i,:))
end

scatter(x_cen, y_cen)

qj_inter = [];
idx = 0;

% Loop over top row
for j=1:11
     if j == 1
        idx = idx + 1;
        qj_inter(idx,1) = qj_des(1,j);
    elseif j == 11
        idx = idx + 1;
        qj_inter(idx,1) = qj_des(1,j-1);    
    else
        idx = idx + 1;
        qj_inter(idx,1) = 0.5*(qj_des(1,j-1) + qj_des(1,j));
    end        
end

% Loop over bottom row
for j=1:11
     if j == 1
        idx = idx + 1;
        qj_inter(idx,1) = qj_des(10,j);
    elseif j == 11
        idx = idx + 1;
        qj_inter(idx,1) = qj_des(10,j-1);    
    else
        idx = idx + 1;
        qj_inter(idx,1) = 0.5*(qj_des(10,j-1) + qj_des(10,j));
    end        
end

% Loop over left column
for i=2:10
     idx = idx + 1;
     qj_inter(idx,1) = 0.5*(qj_des(i,1) + qj_des(i,1));
end

% Loop over right column
for i=2:10
     idx = idx + 1;
     qj_inter(idx,1) = 0.5*(qj_des(i,10) + qj_des(i,10));
end

qj_des_vec = reshape(qj_des, [length(x_cen), 1]);

x_fit = [x_cen; x(1,:)'; x(11,:)'; x(2:10,1); x(2:10, 11)];
y_fit = [y_cen; y(1,:)'; y(11,:)'; y(2:10,1); y(2:10, 11)];
qj_fit = [qj_des_vec; qj_inter];

fit_obj = fit([x_fit, y_fit], qj_fit, 'cubicinterp')
figure;
plot(fit_obj, [x_fit, y_fit], qj_fit)