model = mphload('simple_square.mph');
tol = 1e-16;
err = 1;

% Interpolation of qj
%model.component('comp1').func.create('int1', 'Interpolation');
model.component('comp1').func('int1').set('source', 'file');
model.component('comp1').func('int1').set('filename', 'qj_data.csv');
model.component('comp1').func('int1').set('interp', 'linear');
model.component('comp1').func('int1').set('extrap', 'const');
model.component('comp1').func('int1').set('nargs', 2);
model.component('comp1').func('int1').set('defvars',false);
model.component('comp1').func('int1').set('argunit', 'W/m^2');


% Interpolation of initial thickness
%model.component('comp1').func.create('int2', 'Interpolation');
model.component('comp1').func('int2').set('source', 'file');
model.component('comp1').func('int2').set('filename', 'initial_thickness.csv');
model.component('comp1').func('int2').set('interp', 'linear');
model.component('comp1').func('int2').set('extrap', 'const');
model.component('comp1').func('int2').set('nargs', 2);
model.component('comp1').func('int2').set('defvars',false);
model.component('comp1').func('int2').set('argunit', 'm');

model.result('pg2').feature('surf2').set('expr', 'int1(x,y)*ecs.ds/ecs.Qsrh');
%model.result('pg2').feature('surf2').set('expr', '(ecs.Qsrh/int1(x,y))*ecs.ds');
%model.result('pg3').feature('surf3').set('expr', 'ecs.Qsrh/int1(x,y)');

model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');
model.study('std1').run;

mphsave(model, 'simple_square.mph');
new_delta = mphplot(model, 'pg2');
i = 0;

while err > tol
    i = i + 1;
    updated_delta = [new_delta{1,2}{1,1}.p(1,:)' new_delta{1,2}{1,1}.p(2,:)' new_delta{1,2}{1,1}.d];
    dlmwrite('next_thickness.csv', updated_delta, 'precision', 10);
    current_thickness = new_delta{1,2}{1,1}.d;
    model.component('comp1').func('int2').set('filename', 'next_thickness.csv');
    model.component('comp1').physics('ecs').prop('ds').set('ds', 'int2(x,y)');
    model.study('std1').run;
    mphsave(model, 'simple_square.mph')
    new_delta = mphplot(model, 'pg2');
    next_thickness = new_delta{1,2}{1,1}.d;
    err = sum(abs(current_thickness-next_thickness));
    sprintf('%d Iteration', i)
    sprintf('The error is %d', err)
    pause
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
%}



