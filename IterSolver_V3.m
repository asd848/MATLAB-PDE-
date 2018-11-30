% IterSolver_V3 that runs the iterative solver to find thickness for a desired
% watt density

% file = fopen('GlobalOutput.txt', 'wt');

% System Parameters
voltage = 115;
tolerance = 1;
factor = 50;

% Initial thickness guess
% i_thickness = 3000e-10*ones(10,10);
load('actual_thickness_global.mat');
i_thickness = actual_thickness_global*1e-11;
% Desired watt density
% qjDes = 1.4/0.0254^2*ones(10,10);
load('actual_watt_density_global.mat');
qjDes = actual_watt_density_global./0.0254^2;

% Intial solution
[out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);

% Error calculation
% err = abs(sum(sum(out_thickness - i_thickness)));

% fprintf(file, '%s %d \r\n', 'The applied voltage is', voltage)
% fprintf(file, '%s \r\n', 'The desired watt density is')
% dlmwrite('GlobalOutput.txt', qjDes, 'delimiter', ' ');

fprintf('The applied voltage is %d \n', voltage)
fprintf('The desired watt density is: \n')
disp(qjDes)
i = 0;

fprintf('The intial watt density is: \n')
disp(qj_out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This was us trying a different approach of trying to vary thickness by a
% proportion of the error between qj and qjDesried
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err = qj_out - qjDes;
percent_err = err./qjDes;
err = norm(percent_err);
% fprintf(file, 'The error is %d \n', err)
fprintf('The error is %d \n', err)
fprintf('The element wise error is \n')
disp(qj_out-qjDes);

i_thickness(percent_err>0.1) = i_thickness(percent_err>0.1).*(1-percent_err(percent_err>0.1)./factor);
% disp(percent_err(percent_err>10));
% pause
i_thickness(percent_err<-0.1) = i_thickness(percent_err<-0.1).*(1-percent_err(percent_err<-0.1)./factor);
% disp(percent_err(percent_err<-10))
% pause
% err = abs(sum(sum(out_thickness - i_thickness)));

% fprintf(file, '%s \r\n', 'The updated thickness is')
% dlmwrite('GlobalOutput.txt', i_thickness, 'delimiter', ' ');
fprintf('The updated thickness is: \n')
disp(i_thickness)

% pause
% Begin iterative method
while(err > tolerance)
    i = i +1;
%     i_thickness = out_thickness;
    [out_thickness, qj_out] = busBarGeoSplitFull(voltage, i_thickness, qjDes);
    err = qj_out - qjDes;
    percent_err = err./qjDes;

    i_thickness(percent_err>0.1) = i_thickness(percent_err>0.1).*(1-percent_err(percent_err>0.1)./factor);
    i_thickness(percent_err<-0.1) = i_thickness(percent_err<-0.1).*(1-percent_err(percent_err<-0.1)./factor);
    err = norm(percent_err);
%     err = abs(sum(sum(out_thickness - i_thickness)));
%     fprintf(file, 'Iteration %d \n', i)
%     fprintf(file, 'The error is %d \n', err)
%     fprintf(file, 'The thickness is: \n');
%     dlmwrite('GlobalOutput.txt', i_thickness, 'delimiter', ' ');
%     disp(i_thickness);
%     fprintf('The watt density is: \n');
%     dlmwrite('GlobalOutput.txt', qj_out, 'delimiter', ' ');
%     disp(qj_out)
    
    fprintf('Iteration %d \n', i)
    fprintf('The error is %d \n', err)
    fprintf('The element wise error is \n')
    disp(qj_out-qjDes);
    fprintf('The thickness is: \n')
    disp(i_thickness)
    fprintf('The watt density is: \n')
    disp(qj_out)
%      pause
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
% fclose(file)
