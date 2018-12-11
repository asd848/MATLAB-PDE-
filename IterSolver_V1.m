%% Getting XY Data from a solution
voltage = 115;
sizeSquare = 10;
tolerance = 1e-9;

i_thickness = 1e-7*ones(sizeSquare, sizeSquare);
% qjDes = ones(5,5);
%qjDes = 1.4/0.0254^2*ones(10,10);
load('actual_watt_density_global.mat');
qjDes = actual_watt_density_global./0.0254^2;
[out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
 fprintf('The watt density is \n');
 disp(qj)

err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)

i = 0;

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
