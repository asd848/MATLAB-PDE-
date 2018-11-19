% IterSolver_V1 that runs the iterative solver to find thickness for a desired
% watt density

% System parameters
voltage = 115;
sizeSquare = 5;
tolerance = 1e-9;

% Initial thickness guess
i_thickness = ones(sizeSquare, sizeSquare);

% Desired watt density
qjDes = randi(5,5);

% Initial solution
[out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);

% Error calculation
err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;

% Begin iterative method
while(err > tolerance)
   i = i + 1;
   i_thickness = out_thickness;
   [out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
   err = abs(sum(sum(out_thickness - i_thickness)));
   
   fprintf('Iteration %d \n', i);
   fprintf('This is the error %d\n', err);
   fprintf('The thickness is \n');
   disp(out_thickness);
   fprintf('The watt density is \n');
   disp(qj)
end
