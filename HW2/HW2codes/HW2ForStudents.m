clc;
clear all;
close all;

% input for the example in the lecture notes
E = 70e3; % E-modulus
BC = [1 2 5 6; 0 0 0 0]; %boundary conditions; form: [degrees of freedom; applied displacements]
F = [0 0 0 -1000 0 0]; % force vector
Con = [1 2; 3 2];  %connectivity matrix; form: [element 1; element 2;etc.]
NodePos = [0 0; 125 100; 0 100]; % node positions; form: [node 1; node 2; etc]
Area = [10 10]; % vector with the area for each element

% functions here
[u, strain, stress, reactions] = HW2functionsExample(E,BC,F,Con,NodePos,Area);

% output: do not change!
for el= 1:size(Con,1)
	fprintf(['the stress in element ', num2str((el)) ' equals ', num2str(stress(el)), ' MPa \n'])
end

for node=1:size(NodePos,1)
	fprintf(['The displacement of node ', num2str(node) , ' in x-direction is ', num2str(u(2*(node)-1)), ' mm \n'])
	fprintf(['The displacement of node ', num2str(node) , ' in y-direction is ', num2str(u(2*(node))), ' mm \n'])
	fprintf(['The reaction force at node ', num2str(node) , ' in x-direction is ', num2str(reactions(2*(node)-1)), ' N \n'])
	fprintf(['The reaction force at node ', num2str(node) , ' in y-direction is ', num2str(reactions(2*(node))), ' N \n']) 
end