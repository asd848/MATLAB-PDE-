function [conductance] = conductance_anyPoint(conductivity, voltageMap,x,y,qjDesFit)

    [gradx,grady] = evaluateGradient(voltageMap, x, y);
   	Qj = conductivity .* (gradx.^2 + grady.^2);
    conductance = conductivity .* qjDesFit(x,y)./Qj.*;
end

