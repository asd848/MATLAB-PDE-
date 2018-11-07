%% Getting XY Data from a solution
voltage = 115;
sizeSquare = 2;
tolerance = 1e-7;

i_thickness = ones(sizeSquare, sizeSquare);
qjDes = [10, 30; 20, 40];

out_thickness = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);

err = abs(sum(sum(out_thickness - i_thickness)));

while(err > tolerance)
   fprintf('This is the error %d\n', err);
   i_thickness = out_thickness;
   out_thickness = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
   err = abs(sum(sum(out_thickness - i_thickness)));
%    disp(out_thickness);
end
