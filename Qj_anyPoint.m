function [Qj] = Qj_anyPoint(conductivity, voltageMap,x,y,qjDesFit)

    [gradx,grady] = evaluateGradient(voltageMap, x, y);
   	Qj = conductivity .* (gradx.^2 + grady.^2);
    conductance = conductivity .* qjDesFit(x,y)./Qj.*sqrt((x-[x(2:end) 0]).^2+(y-[y(2:end) 0)^2)
end

