function [u, strain, stress, reactions] = HW2functionsExample(E,BC,F,Con,NodePos,Area)

u = zeros(1,2*size(NodePos,1));

strain = calcStrain(size(Con,2)); % just an example of a function, input is also just to have an input
stress = calcStress(size(Con,2)); % just an example of a function, input is also just to have an input

reactions = -F; % just a placeholder to have a vector of the right size

end

function strain = calcStrain(Nelem)
    strain = ones(1,Nelem);
end

function stress = calcStress(Nelem)
    stress = 2*ones(1,Nelem);
end