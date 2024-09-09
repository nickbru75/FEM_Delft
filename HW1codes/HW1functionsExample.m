function [u, strain, stress] = HW1functionsExample(L, Nelem, F, E, BC);
% this is where the code is starting, all functions defined below are
% called from this function 

u = zeros(1,Nelem);

strain = calcStrain(Nelem); % just an example of a function, input is also just to have an input
stress = calcStress(Nelem); % just an example of a function, input is also just to have an input

end

function strain = calcStrain(Nelem)
    strain = ones(1,Nelem);
end

function stress = calcStress(Nelem)
    stress = 2*ones(1,Nelem);
end