%% Getting XY Data from a solution
voltage = 1;
sizeSquare = 2;

[xData, yData] = meshgrid(1:sizeSquare, 1:sizeSquare);
mesh = [xData(:) yData(:)];
xData = mesh(1:end,1);
yData = mesh(1:end,2);


[model, results] = jouleHeater(sizeSquare,voltage);



uintrp = interpolateSolution(results,xData, yData);
% uintrp(2,1:end) = interpolateSolution(results,xDataReverse,yDataReverse);

[gradx, grady] = evaluateGradient(results, xData, yData);
% [gradx(1:end),grady(2,1:end)] = evaluateGradient(results, xDataReverse, yDataReverse)';