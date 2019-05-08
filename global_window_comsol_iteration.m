%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global_window_comsol_iteration.m
% Script for running COMSOL MATLAB iterative method via LiveLink for the
% Global window geometry
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main 

%{
Load the model using mphload command. NOTE: COMSOL model must be in the
same directory as MATLAB script or you must provide the full path to the
model
%}
model = mphload('window1.mph');

% System parameters
tol = 10.123;

%% Interpolating data with COMSOL from csv file


% If you have not created an interpolation function in COMSOL you must
% run the following command 
% model.component('comp1').func.create('int1', 'Interpolation');

% Interpolation of desired joule heating
model.component('comp1').func('int1').set('source', 'file');
model.component('comp1').func('int1').set('filename', 'qj_desired_global.csv');
model.component('comp1').func('int1').set('interp', 'linear');
model.component('comp1').func('int1').set('nargs', 2);
model.component('comp1').func('int1').set('defvars',false);
model.component('comp1').func('int1').set('argunit', 'W/m^2');


% Interpolation of initial thickness
% Again, if you have not created a second interpolation function in COMSOL
% you must run the following command
% model.component('comp1').func.create('int2', 'Interpolation');

model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('filename', 'initial_thickness_global.csv');
model.component('comp1').func('int2').set('interp', 'linear');
model.component('comp1').func('int2').set('nargs', 2);
model.component('comp1').func('int2').set('defvars',false);
model.component('comp1').func('int2').set('argunit', 'm');

% Assumes you have already created plot groups with a surface plot within

% Setting the plot result to voltage
model.result('pg1').feature('surf1').set('expr', 'V');

% Setting the plot result to surface joule heating
model.result('pg2').feature('surf1').set('expr', 'ecs.Qsrh');

% Set the ECS thickness to our interpolated thickness data
model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');

% Run the study
model.study('std1').run;

% Save the model
mphsave(model, 'window1.mph');

% Grabbing interpolation co-ordiantes to ensure we compare the same co-ords 
% when updating
qj_des = mphplot(model, 'pg3');
xy_coords = qj_des{1,1}{1,1}.p(1:2,:);
xy_coords = [xy_coords; 0.0127*ones(1, length(xy_coords))];
qj_des = qj_des{1,1}{1,1}.d;

% Writing desired qj data
qj_des_data = [xy_coords' qj_des];

% Removing any NaN values
indices = find(isnan(qj_des_data(:,4)));
j=0;
for i=1:length(qj_des_data)
    if ismember(i, indices) == 0
        j = j + 1;
        qj_des_data_copy(j,:) = qj_des_data(i,:);
    end
end

% Saving COMSOL's interpolation of desired joule heating to a csv
xy_coords = qj_des_data_copy(:,1:3);
qj_des = qj_des_data_copy(:,4);
dlmwrite('qj_des_data.csv', qj_des_data_copy, 'precision', 10);

% Grabbing interpolation of intial thickness estimate
current_delta = mphplot(model, 'pg4');
current_delta = current_delta{1,1}{1,1}.d;

% Interpolating joule heating from solution using interpolation co-ordinates
current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords', 'edim', 'boundary', 'selection',1);
current_qj = current_qj';

% Remove any NaN values
current_qj_data = [xy_coords current_qj];
indices = find(isnan(current_qj_data(:,4)));
j=0;
for i=1:length(current_qj_data)
    if ismember(i, indices) == 0
        j = j + 1;
        qj_des_data_copy_2(j,:) = qj_des_data_copy(i,:);
        current_qj_data_copy(j,:) = current_qj_data(i,:);
    end
end

% Saving current iteration's joule heating to a csv file
xy_coords = current_qj_data_copy(:,1:3);
current_qj = current_qj_data_copy(:,4);
current_qj_data = [xy_coords(:,1:2) current_qj];
qj_des =  qj_des_data_copy_2(:,4);
dlmwrite('current_qj.csv', current_qj_data, 'precision', 10);

% Saving the current iterations thickness (initial thickness) to a csv file
current_delta_data = [xy_coords(:,1:2) 3.18E-7*ones(length(xy_coords),1)];
dlmwrite('current_delta.csv', current_delta_data, 'precision', 10);
current_delta = current_delta_data(:,3);

% Calculating the error
sprintf('Iteration 1')
err = mean(abs(current_qj-qj_des));
sprintf('The error is %g', err)

% Update thickness and save to csv as the next iteration's thickness
updated_delta = (qj_des./current_qj).*current_delta;
updated_delta_data = [xy_coords(:,1:2) updated_delta];
dlmwrite('next_thickness.csv', updated_delta_data, 'precision', 10);

%% Begin iteration
i = 0;

while err >= tol
    i = i + 1;
    
    % Grab current delta and write to the current delta csv file
    current_delta = updated_delta;
    current_delta_data = [xy_coords(:,1:2) current_delta];
    dlmwrite('current_delta.csv', current_delta_data, 'precision', 10);
    
    % Upload new thickness to COMSOL interpolation
    model.component('comp1').func('int2').set('filename', 'next_thickness.csv');
    model.component('comp1').func('int2').set('extrap', 'const');
    
    % Set shell thickness to new thickness
    model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');
    
    % Run and svae the model
    model.study('std1').run;
    mphsave(model, 'window1.mph')
    
    % Interpolate solution for qj and save to a csv file
    current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords', 'edim', 'boundary', 'selection',1);
    current_qj = current_qj';
    current_qj_data = [xy_coords(:,1:2) current_qj];
    dlmwrite('current_qj.csv', current_qj_data, 'precision', 10);
    
    % Update delta using Resisty updating method
    updated_delta = resisty_updating(i, current_qj, qj_des, current_delta);
     
    % Pulls the x, y, thickness data from the results
    updated_delta_data = [xy_coords(:,1:2) updated_delta];
    
    % Writing udpated thickness data points to csv
    dlmwrite('next_thickness.csv', updated_delta_data, 'precision', 10);

    % Calculate error
    err = mean(abs(current_qj-qj_des));
    
    sprintf('%d Iteration', i)
    sprintf('The error is %d', err)
    disp('mean')
    disp(mean(current_qj-qj_des));
    disp('std')
    disp(std(current_qj-qj_des));
end
end

%% Resisty updating method
function [updated_delta] = resisty_updating(i, current_qj, qj_des, current_delta)
% INPUTS:
% i = iteration count
% current_qj = ith joule heating 
% qj_des = desired joule heating
% current_delta = ith delta profile
% OUTPUTS:
% updated_delta = new delta profile

   if mod(i,4) == 1
       updated_delta = current_delta.*(current_qj./qj_des);
   elseif mod(i,4) == 2
       updated_delta = current_delta.*(qj_des./current_qj).^2;
   elseif mod(i,4) == 3
       updated_delta = current_delta.*(current_qj./qj_des).^2; 
   else
       updated_delta = current_delta.*(qj_des./current_qj);
   end
end

