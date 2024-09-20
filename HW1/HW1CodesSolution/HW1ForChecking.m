clc;
clear all;
close all;

% input
L =  250; % total length of the bar
Nelem = 4; % number of elements
F =  1000; % applied force on the final nodes
E =  70e3; % elastic modulus 
BC = [1; 0]; %form: [node numbers; applied displacements]

% functions here
[u, strain, stress] = HW1functions(L, Nelem, F, E, BC);

% output
for el=1:Nelem
	fprintf(['the stress in element ', num2str(el) , ' equals ', num2str(stress(el)), ' MPa \n'])
end
	
fprintf(['The displacement of the final node is ', num2str(u(end)),' mm \n'])