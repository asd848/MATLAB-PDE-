model = mphload('simple_square.mph');
mphgetproperties(model.geom('geom1'))

figure;
mphgeom(model, 'geom1', 'facelabels', 'on');

figure;
mphmesh(model);

%%
figure;
voltage = mphplot(model, 'pg1');
figure;
volumetric_heating = mphplot(model, 'pg2');




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

%}
