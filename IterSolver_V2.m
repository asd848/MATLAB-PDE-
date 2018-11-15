voltage = 115;
tolerance = 1e-9;

i_thickness = ones(2, 2);
qjDes = 5*randi(2,2);
[out_thickness, qj_out] = busBarGeoSplit(voltage, i_thickness, qjDes);

err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;

while(err > tolerance)
   i = i +1;
   i_thickness = out_thickness;
   [out_thickness, qj_out] = busBarGeoSplit(voltage, i_thickness, qjDes);
   err = abs(sum(sum(out_thickness - i_thickness)));
   fprintf('Iteration %d \n', i);
   fprintf('This is the error %d\n', err);
   fprintf('The thickness is \n');
   disp(out_thickness);
   fprintf('The watt density is \n');
   disp(qj_out)
end
