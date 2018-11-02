%% Getting XY Data from a solution
voltage = 115;
sizeSquare = 2;


numMeasurements = 10;
%Assume that the length of each square is 1
spacing = 1/numMeasurements;

%% For each square
%{
0.1 0.2 0.3 0.4 0.5  -- 1
1.1 1.2 1.3 1.4 1.5  -- 2

%}

%% So far this taking 10 points for each square (for average voltage calculation)
[xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
mesh = [xData(:) yData(:)];
xData = mesh(1:end,1);
yData = mesh(1:end,2);



[model, results] = jouleHeater(sizeSquare,voltage);



uintrp = interpolateSolution(results,xData, yData);
% uintrp(2,1:end) = interpolateSolution(results,xDataReverse,yDataReverse);

[gradx, grady] = evaluateGradient(results, xData, yData);

xData = unique(xData);
yData = unique(yData);

normV = zeros(sizeSquare*numMeasurements,sizeSquare*numMeasurements);

for i = 1:length(xData)
    for j  = 1:length(yData)
        normV(i,j) = gradx(i).^2 + grady(j).^2;
    end
end





%% Find the average gradient for a specific square


% [gradx(1:end),grady(2,1:end)] = evaluateGradient(results, xDataReverse, yDataReverse)';