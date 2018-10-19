%% Getting XY Data from a solution
voltage = 1;
sizeSquare = 2;

[xData, yData] = meshgrid(0.5:1:sizeSquare, 0.5:1:sizeSquare);
mesh = [xData(:) yData(:)];
xData = mesh(1:end,1);
yData = mesh(1:end,2);


[model, results] = jouleHeater(sizeSquare,voltage);



uintrp = interpolateSolution(results,xData, yData);
% uintrp(2,1:end) = interpolateSolution(results,xDataReverse,yDataReverse);

[gradx, grady] = evaluateGradient(results, xData, yData);
normV = gradx.^2 + grady.^2;

% [gradx(1:end),grady(2,1:end)] = evaluateGradient(results, xDataReverse, yDataReverse)';