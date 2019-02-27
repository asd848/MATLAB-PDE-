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
model = mphload('simple_square.mph');

% Global parameters for iterative method
alpha = 1e-3;
tol = 1e-16;
err = 1;

%% Interpolating data with COMSOL from csv file

% Interpolation of qj
%{
If you have not created an interpolation function in COMSOL yet you must
run 
%}

%model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func('int1').set('source', 'file'); % 
model.component('comp1').func('int1').set('filename', 'qj_test.csv');
model.component('comp1').func('int1').set('interp', 'linear');
model.component('comp1').func('int1').set('extrap', 'const');
model.component('comp1').func('int1').set('nargs', 2);
model.component('comp1').func('int1').set('defvars',false);
model.component('comp1').func('int1').set('argunit', 'W/m^2');

% Interpolation of initial thickness
%model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('filename', 'uniform_inital_thickness.csv');
model.component('comp1').func('int2').set('interp', 'linear');
model.component('comp1').func('int2').set('extrap', 'const');
model.component('comp1').func('int2').set('nargs', 2);
model.component('comp1').func('int2').set('defvars',false);
model.component('comp1').func('int2').set('argunit', 'm');

% model.result('pg2').feature('surf2').set('expr', 'int1(x,y)*ecs.ds/ecs.Qsrh');
% model.result('pg2').feature('surf2').set('expr', '(sqrt(ecs.Qsrh/int1(x,y)))*ecs.ds');
%model.result('pg3').feature('surf3').set('expr', 'ecs.Qsrh/int1(x,y)');
model.result('pg2').feature('surf2').set('expr', 'ecs.Qsrh');
model.result('pg3').feature('surf1').set('expr', 'ecs.ds');

% Set the ECS thickness to our interpolated thickness data
model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');

% Run the study
model.study('std1').run;

% Save the model
mphsave(model, 'simple_square.mph');

% Grabbing interpolation co-ordiantes and reshaping qj_des
qj_des = mphplot(model, 'pg4');
xy_coords = qj_des{1,1}{1,1}.p(1:2,:);
xy_coords = [xy_coords; 0.0127*ones(1, length(xy_coords))];
qj_des = qj_des{1,1}{1,1}.d;

% Grabbing interpolation of intial thickness estimate
current_delta = mphplot(model, 'pg5');
current_delta = current_delta{1,2}{1,1}.d;

% Interpolating qj from solution using interpolation co-ordinates
current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords, 'edim', 'boundary', 'selection',4);
current_qj = current_qj';

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
updated_delta = current_delta + current_delta.*(current_qj./qj_des -1)*alpha;
% percent_err = (qj_des-current_qj)./qj_des.*alpha;

% Hanna's Updating
% updated_delta = current_delta.*(1 - (qj_des-current_qj)./qj_des.*alpha);

%%
% Iteration counter
i = 1;
% Begin iterating
while err > tol
    i = i + 1;
    % Pulls the x, y, thickness data from the results
    updated_delta_data = [xy_coords(1:2,:)' updated_delta];
    % Writing udpated thickness data points to csv
    dlmwrite('next_thickness.csv', updated_delta_data, 'precision', 10);
    % Grab current delta
    current_delta = updated_delta;
    % Upload new thickness to COMSOL interpolation
    model.component('comp1').func('int2').set('filename', 'next_thickness.csv');
    % Set shell thickness to new thickness
    model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');
    % Run and svae the model
    model.study('std1').run;
    mphsave(model, 'simple_square.mph')
%     figure(1);
%     mphplot(model, 'pg1');
%     figure(2);
%     mphplot(model, 'pg2');
%     figure(3);
%     mphplot(model, 'pg3');
    % Interpolate solution for qj
    current_qj = mphinterp(model, 'ecs.Qsrh', 'coord', xy_coords, 'edim', 'boundary', 'selection',4);
    current_qj = current_qj';
    disp(current_qj(1,1))
    disp(current_qj(5000,1))
    disp(current_qj(10000,1))
    disp(mean(current_qj-qj_des))
    % Update delta
    updated_delta = current_delta + current_delta.*(current_qj./qj_des - 1)*alpha;
%     updated_delta = current_delta.*(1 - (qj_des-current_qj)./qj_des.*alpha);
    % Calculate error
    err = sum(sum(abs(qj_des-current_qj)));
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



