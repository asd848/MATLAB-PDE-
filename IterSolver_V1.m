%% Getting XY Data from a solution
voltage = 115;
sizeSquare = 10;
tolerance = 1e-9;

i_thickness = ones(sizeSquare, sizeSquare);
% qjDes = ones(5,5);
%qjDes = 1.4/0.0254^2*ones(10,10);
load('actual_watt_density_global.mat');
qjDes = actual_watt_density_global./0.0254^2;
[out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
 fprintf('The watt density is \n');
 disp(qj)

% Old error checking
% err = abs(sum(sum(out_thickness - i_thickness)));
% fprintf('The applied voltage is %d V \n', voltage)
% fprintf('The desired watt density is \n')
% disp(qjDes)

err = qj_out - qjDes;
percent_err = err./qjDes * 100;

i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10)./1000000);
i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10)./1000000);
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)

i = 0;

while(err > tolerance)
   i = i + 1;
%    i_thickness = out_thickness;
%    [out_thickness, qj] = jouleHeater(sizeSquare,voltage, i_thickness, qjDes);
%    err = abs(sum(sum(out_thickness - i_thickness)));
%    fprintf('Iteration %d \n', i);
%    fprintf('This is the error %d\n', err);
%    fprintf('The thickness is \n');
%    disp(out_thickness);
%    fprintf('The watt density is \n');
%    disp(qj)

   [out_thickness, qj_out] = jouleHeater(sizeSquare, voltage, i_thickness, qjDes);
   
   err = qj_out - qjDes;
   percent_err = err./qjDes * 100;

   i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10));
   i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10));
   
   i_thickness(i_thickness>8000e-10) = 8000e-10;
   i_thickness(i_thickness<300e-10) = 300e-10;
   fprintf('Iteration %d \n', i);
   fprintf('The error is: %d\n', abs(mean(mean(err))));
   fprintf('The thickness is: \n');
   disp(i_thickness);
   fprintf('The watt density is: \n');
   disp(qj_out)
end
