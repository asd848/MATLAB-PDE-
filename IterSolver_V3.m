voltage = 115;
tolerance = 1e-16;

i_thickness = 3000e-10*ones(10,10);
qjDes = 1.4/0.0254^2*ones(10,10);
% qjDes = global_heating./(0.0254^2);
initialRun = true;
[out_thickness, qj_out, voltageMap] = busBarGeoSplitFull(voltage, i_thickness, qjDes, initialRun, 0);
initialRun = false;
err = qj_out - qjDes;
percent_err = err./qjDes * 100;

% for i=1:10
%     for j=1:10
%         if percent_err(i,j) > 12
%             i_thickness(i,j) = i_thickness(i,j)./abs(percent_err(i,j));
%         elseif percent_err(i,j) < -10
%             i_thickness(i,j) = i_thickness(i,j).*abs(percent_err(i,j));
%         end
%     end
% end



i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10)./1000000);
% disp(percent_err(percent_err>10));
% pause
i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10)./1000000);
% disp(percent_err(percent_err<-10))
% pause
%err = abs(sum(sum(out_thickness - i_thickness)));
% err = norm(qj_out-qjDes);
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;
% pause
while(abs(mean(mean(err))) > tolerance)
   i = i +1;
%    i_thickness = out_thickness;
   [out_thickness, qj_out, voltageMap] = busBarGeoSplitFull(voltage, i_thickness, qjDes,initialRun, voltageMap);
   
   err = qj_out - qjDes;
   percent_err = err./qjDes * 100;

   i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10));
   i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10));
   
   i_thickness(i_thickness>8000e-10) = 8000e-10;
   i_thickness(i_thickness<300e-10) = 300e-10;
   
   
%    err = norm(qj_out-qjDes);
%    err = abs(sum(sum(out_thickness - i_thickness)));
   fprintf('Iteration %d \n', i);
   fprintf('The error is: %d\n', abs(mean(mean(err))));
   fprintf('The thickness is: \n');
   disp(i_thickness);
   fprintf('The watt density is: \n');
   disp(qj_out)
end