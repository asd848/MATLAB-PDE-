%% Getting XY Data from a solution
voltage = 115;
sizeSquare = 2;


numMeasurements = 10;
%Assume that the length of each square is 1/sizeSquare
spacing = 1/numMeasurements;

%% For each square
%{
0.1 0.2 0.3 0.4 0.5  -- 1
1.1 1.2 1.3 1.4 1.5  -- 2

%}

%% So far this taking 10 points for each square (for average voltage calculation)
% [xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
% mesh = [xData(:) yData(:)];
% xData = mesh(1:end,1);
% yData = mesh(1:end,2);



[model, results] = jouleHeater(sizeSquare,voltage);

[xData, yData] = meshgrid(spacing:spacing:sizeSquare, spacing:spacing:sizeSquare);
mesh = [xData(:) yData(:)];
xData = mesh(1:end,1);
yData = mesh(1:end,2);

uintrp = interpolateSolution(results,xData, yData);

[gradx, grady] = evaluateGradient(results, xData, yData);

gradx = sizeSquare.*gradx;
grady = sizeSquare.*grady;

xData = unique(xData);
yData = unique(yData);

Qj= zeros(sizeSquare*numMeasurements,sizeSquare*numMeasurements);

for i = 1:length(xData)
    for j  = 1:length(yData)
        Qj(i,j) = gradx(i).^2 + grady(j).^2;
    end
end

%% Calculate average Qj for each square

avgQj = zeros(sizeSquare,sizeSquare);
for i = 1:sizeSquare
    for j = 1:sizeSquare
        xStartIndex = numMeasurements * (i -1) + 1;
        xEndIndex = numMeasurements * i;
        yStartIndex = numMeasurements * (j -1) + j;
        yEndIndex = numMeasurements * j;

        avgQj(i,j) = mean(mean(Qj(xStartIndex:xEndIndex,yStartIndex:yEndIndex)));
    end
end




