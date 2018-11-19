% IterSolver_V3 that runs the iterative solver to find thickness for a desired
% watt density

% System Parameters
voltage = 115;
tolerance = 1e-16;

% Initial thickness guess
i_thickness = 3000e-10*ones(10,10);

% Desired watt density
qjDes = 1.4/0.0254^2*ones(10,10);

% Intial solution
[out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);

% Error calculation
err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)
i = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was us trying a different approach of trying to vary thickness by a
% proportion of the error between qj and qjDesried
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% err = qj_out - qjDes;
% percent_err = err./qjDes * 100;
% 
% i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10)./1000000);
% disp(percent_err(percent_err>10));
% pause
% i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10)./1000000);
% disp(percent_err(percent_err<-10))
% pause
%err = abs(sum(sum(out_thickness - i_thickness)));
% err = norm(qj_out-qjDes);

% Begin iterative method
while(abs(mean(mean(err))) > tolerance)
    i = i +1;
    i_thickness = out_thickness;
    [out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);
    err = abs(sum(sum(out_thickness - i_thickness)));
    fprintf('Iteration %d \n', i);
    fprintf('The error is: %d\n', abs(mean(mean(err))));
    fprintf('The thickness is: \n');
    disp(i_thickness);
    fprintf('The watt density is: \n');
    disp(qj_out)
    
%%%%%%%%%%%%%%%%%%%%%    
%    Other Approach
%%%%%%%%%%%%%%%%%%%%%
%    err = qj_out - qjDes;
%    percent_err = err./qjDes * 100;
% 
%    i_thickness(percent_err>10) = i_thickness(percent_err>10) ./abs(percent_err(percent_err>10));
%    i_thickness(percent_err<-10) = i_thickness(percent_err<-10) .*abs(percent_err(percent_err<-10));
%    
%    i_thickness(i_thickness>8000e-10) = 8000e-10;
%    i_thickness(i_thickness<300e-10) = 300e-10; 


end
