%% Getting XY Data from a solution
voltage = 115;

%% The number of squares along one edge
sizeSquare = 10;

%% Tolerance for iteration termination
tolerance = 1e-16;

%% Initial thickness guess, (uniform)
i_thickness = 3.78e-8*ones(sizeSquare, sizeSquare);
%i_thickness = 1.0e-6*[.0199 .0441 .0789 .1231 .1717; .0166 .0367 .0657 .1026 .1431; .0133 .0294 .0526 .0821 .1145; .01 .022 .0394 .0616 .0858; .0066 .0147 .0263 .041 .0572];
%% Desired surcace heating of the global window in matrix form
load('qj_test.csv');
qjDes = qj_test;
%qjDes = 3000:-2000/9:1000;
%qjDes = [qjDes' qjDes' qjDes' qjDes' qjDes' qjDes' qjDes' qjDes' qjDes' qjDes'];

figure(1);
contour(qjDes);
set(gca, 'FontSize', 12);
xlabel('x');
ylabel('y');
title('Contour Plot of Desired Heating');
contourcbar;


figure(2);
image(i_thickness);
set(gca, 'FontSize', 12);
xlabel('x');
ylabel('y');
title('Contour Plot of Initial Thickness');


%% Initial call to jouleHeater to get a matrix of thicknesses and surface joule heating
[qj] = resisty_joule_heating(sizeSquare,voltage, i_thickness);
fprintf('The watt density is \n');
disp(qj)

out_thickness = i_thickness*(qjDes./qj);

figure(3);
contour(out_thickness);
set(gca, 'FontSize', 12);
xlabel('x');
ylabel('y');
title('Contour Plot of Iteration 1 Thickness Profile');
contourcbar;

%% Initial error calculation
err = abs(sum(sum(out_thickness - i_thickness)));
fprintf('The applied voltage is %d V \n', voltage)
fprintf('The desired watt density is \n')
disp(qjDes)

i = 0;

%% Iterations terminate if error drops below tolerance
while(err > tolerance)
   % count iterations
   i = i + 1;
   
   % update input thickness to be previous iterations output thickness
   i_thickness = out_thickness;
   
   % call to jouleHeater to solve for new thicknesses
   [qj] = resisty_joule_heating(sizeSquare,voltage, i_thickness);
   
   if mod(i,4) == 1
       out_thickness = i_thickness.*(qj./qjDes);
   elseif mod(i,4) == 2
       out_thickness = i_thickness.*(qjDes./qj).^2;
   elseif mod(i,4) == 3
       out_thickness = i_thickness.*(qj./qjDes).^2; 
   else
       out_thickness = i_thickness.*(qjDes./qj);
   end
   
   % recalculate the error between thickness iterations currently looking 
   % at the sum of all the entries in the matrix as our metric
   err = abs(sum(sum(out_thickness - i_thickness)));
   fprintf('Iteration %d \n', i);
   fprintf('This is the error %d\n', err);
   fprintf('The thickness is \n');
   disp(out_thickness);
   fprintf('The watt density is \n');
   disp(qj)
end
