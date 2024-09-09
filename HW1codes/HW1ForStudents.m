clc;
clear all;
close all;

% input
L =  ; % total length of the bar (scalar value)
Nelem = ; % number of elements (scalar value)
F =  ; % applied force on the final nodes (scalar value)
E =  ; % elastic modulus (scalar value)
BC = ; %form: [node numbers; applied displacements]; each row (i.e., node numbers, applied displacements) can have multiple inputs
% example BC: [i; x] should apply a displacement of x on node i 
% example 2 BC: BC = [i j; x y] should apply a displacement of x on node i and a displacement of y on onde j

% call functions file 
[u, strain, stress] = HW1functionsExample(L, Nelem, F, E, BC);

% output
for el=1:Nelem
	fprintf(['the stress in element ', num2str(el) , ' equals ', num2str(stress(el)), ' MPa \n'])
end
	
fprintf(['The displacement of the final node is ', num2str(u(end)),' mm \n'])