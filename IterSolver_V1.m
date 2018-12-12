%% Getting XY Data from a solution
voltage = 115;
%% The number of squares along one edge
sizeSquare = 10;
%% Tolerance for iteration termination
tolerance = 1e-9;
%% Initial thickness guess, (uniform)
i_thickness = ones(sizeSquare, sizeSquare);
%% Desired surcace heating of the global window in matrix form
qjDes = 1.4/0.0254^2*ones(10,10);
%% initial call to jouleHeater to get a matrix of thicknesses and surface joule heating
[out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
 fprintf('The watt density is \n');
 disp(qj)
%% Initial error calculation
err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;

%% Iterations terminate if error drops below tolerance
while(err > tolerance)
   % count iterations
   i = i + 1;
   % update input thickness to be previous iterations output thickness
   i_thickness = out_thickness;
   % call to jouleHeater to solve for new thicknesses
   [out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
   % recalculate the error between thickness iterations currently looking 
   % at the sum of all the entries in the matrix as our metric
   err = abs(sum(sum(out_thickness - i_thickness)));
   fprintf('Iteration %d \n', i);
   fprintf('This is the error %d\n', err);
   fprintf('The thickness is \n');
   disp(out_thickness);
   fprintf('The watt density is \n');
   disp(qj)
end
