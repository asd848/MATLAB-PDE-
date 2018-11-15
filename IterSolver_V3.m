voltage = 115;
tolerance = 1e-16;

i_thickness = 50e-9*ones(10,10);
qjDes = 3000*ones(10,10);
[out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);

%err = abs(sum(sum(out_thickness - i_thickness)));
err = norm(qj_out-qjDes);
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;

while(err > tolerance)
   i = i +1;
   i_thickness = out_thickness;
   [out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);
   err = norm(qj_out-qjDes);
%    err = abs(sum(sum(out_thickness - i_thickness)));
   fprintf('Iteration %d \n', i);
   fprintf('The error is: %d\n', err);
   fprintf('The thickness is: \n');
   disp(out_thickness);
   fprintf('The watt density is: \n');
   disp(qj_out)
end
