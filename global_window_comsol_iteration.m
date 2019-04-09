%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Iterative method in COMSOL using LiveLink
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
Load the model using mphload command. NOTE: COMSOL model must be in the
same directory as MATLAB script or you must provide the full path to the
model
%}
function main 
model = mphload('window1.mph');

% Global parameters for iterative method
alpha = 1e-2;
tol =1e-16;
err = 1;

%% Interpolating data with COMSOL from csv file

% Interpolation of qj
%{
If you have not created an interpolation function in COMSOL yet you must
run 
%}

% model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func('int1').set('source', 'file'); % 
model.component('comp1').func('int1').set('filename', 'qj_desired_global.csv');
model.component('comp1').func('int1').set('interp', 'linear');
model.component('comp1').func('int1').set('nargs', 2);
model.component('comp1').func('int1').set('defvars',false);
model.component('comp1').func('int1').set('argunit', 'W/m^2');

% Interpolation of initial thickness
% model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('filename', 'initial_thickness_global.csv');
model.component('comp1').func('int2').set('interp', 'linear');
model.component('comp1').func('int2').set('nargs', 2);
model.component('comp1').func('int2').set('defvars',false);
model.component('comp1').func('int2').set('argunit', 'm');

model.result('pg1').feature('surf1').set('expr', 'V');
model.result('pg2').feature('surf1').set('expr', 'ecs.Qsrh');

% Set the ECS thickness to our interpolated thickness data
model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');

% Run the study
model.study('std1').run;

% Save the model
mphsave(model, 'window1.mph');

% Grabbing interpolation co-ordiantes and reshaping qj_des
qj_des = mphplot(model, 'pg3');
xy_coords = qj_des{1,1}{1,1}.p(1:2,:);
xy_coords = [xy_coords; 0.0127*ones(1, length(xy_coords))];
qj_des = qj_des{1,1}{1,1}.d;
% Writing desired qj data
qj_des_data = [xy_coords' qj_des];
indices = find(isnan(qj_des_data(:,4)));
j=0;
for i=1:length(qj_des_data)
    if ismember(i, indices) == 0
        j = j + 1;
        qj_des_data_copy(j,:) = qj_des_data(i,:);
    end
end
xy_coords = qj_des_data_copy(:,1:3);
qj_des = qj_des_data_copy(:,4);
dlmwrite('qj_des_data.csv', qj_des_data_copy, 'precision', 10);

% Grabbing interpolation of intial thickness estimate
current_delta = mphplot(model, 'pg4');
current_delta = current_delta{1,1}{1,1}.d;
% Writing current delta data 

% Interpolating qj from solution using interpolation co-ordinates
current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords', 'edim', 'boundary', 'selection',1);
current_qj = current_qj';
% Writing current qj data
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
xy_coords = current_qj_data_copy(:,1:3);
current_qj = current_qj_data_copy(:,4);
current_qj_data = [xy_coords(:,1:2) current_qj];
qj_des =  qj_des_data_copy_2(:,4);
dlmwrite('current_qj.csv', current_qj_data, 'precision', 10);
current_delta_data = [xy_coords(:,1:2) 3.18E-7*ones(length(xy_coords),1)];
dlmwrite('current_delta.csv', current_delta_data, 'precision', 10);
current_delta = current_delta_data(:,3);


% figure(1);
% mphplot(model, 'pg1');
% figure(2);
% mphplot(model, 'pg2');
% figure(3);
% mphplot(model, 'pg3');

% Calculating the error
sprintf('Iteration 1')
err = sum(sum(abs(qj_des-current_qj)));
sprintf('The error is %g', err)

% Update Thickness
% Nick's updating
%updated_delta = current_delta + current_delta.*(current_qj./qj_des -1)*alpha;
updated_delta = (qj_des./current_qj).*current_delta;
updated_delta_data = [xy_coords(:,1:2) updated_delta];
dlmwrite('next_thickness.csv', updated_delta_data, 'precision', 10);
% percent_err = (qj_des-current_qj)./qj_des.*alpha;

% Hanna's Updating
% updated_delta = current_delta.*(1 - (qj_des-current_qj)./qj_des.*alpha);

%%
% Iteration counter
i = 0;
% Begin iterating
while err > tol || mean(current_qj-qj_des) >= .02 %&& i < 80
    i = i + 1;
    % Grab current delta
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
    % Interpolate solution for qj
    current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords', 'edim', 'boundary', 'selection',1);
    current_qj = current_qj';
    current_qj_data = [xy_coords(:,1:2) current_qj];
    dlmwrite('current_qj.csv', current_qj_data, 'precision', 10);
    disp(mean(current_qj-qj_des))
    disp(std(current_qj-qj_des))
    disp(max(current_qj-qj_des))
    disp(min(current_qj-qj_des))
    percent_error = abs(qj_des-current_qj)./qj_des;
    % Update delta
    updated_delta = suspect_updating(i, current_qj, qj_des, current_delta);
     
    % Pulls the x, y, thickness data from the results
    updated_delta_data = [xy_coords(:,1:2) updated_delta];
    % Writing udpated thickness data points to csv
    dlmwrite('next_thickness.csv', updated_delta_data, 'precision', 10);
%     updated_delta = current_delta.*(1 - (qj_des-current_qj)./qj_des.*alpha);
    % Calculate error

    err = sum(abs(updated_delta-current_delta)); 
    sprintf('%d Iteration', i)
    sprintf('The error is %d', err)
end

%% Commands so far
%{
To plot things:
model.study gives you the child node std1
model.study('std1').run executes the program

The results are listed under several children, called 'pg1', 'pg2',etc.
Plot these children by using the command mphplot(model,'pg1')

Command chain for getting properties of boundary condition 
mphgetproperties(model.physics('hteq#').feature('flux1'))

Set the properties of the boundary conditions
model.physics('hteq#').feature('flux1').set('q', 1, '#','g',1,'#')

% figure;
% mphgeom(model, 'geom1', 'facelabels', 'on');
% 
% figure;
% mphmesh(model);

%%
% figure;
% voltage = mphplot(model, 'pg1');

% figure;
% volumetric_heating = mphplot(model, 'pg2');
% mphmesh(model);
% [meshStats, meshData] = mphmeshstats(model);
%}
end
%% Suspect Updating
function [updated_delta] = suspect_updating(i, current_qj, qj_des, current_delta)
   
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

