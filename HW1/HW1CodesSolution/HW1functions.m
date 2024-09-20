function [uTotal, strain, stress] = HW1functions(L, Nelem, F, E, BC)

Fvec = zeros(Nelem+1,1); % create force vector of right size
Fvec(end) = F; % put force on final node

NodePos = 0:L/Nelem:L; % calculate node positions

Con = [1:Nelem;2:Nelem+1]'; % define connectivity matrix for this case

Area = CalcArea(NodePos,L); % call function for area calculation

len = NodePos(2:end)-NodePos(1:end-1); % calculated length of each element

K = StiffMat(Con,Area,E,len); % call function to build global stiffness matrix

remov = BC(1,:); % these are the degrees of freedom (nodes) that need to be removed
keep = setdiff(1:length(NodePos),remov); % all other degrees of freedom remain 

Kred = K(keep,keep); % reduced stiffness matrix with only free degrees of freedom: the ones linked to boundary conditions have been removed

Freduced = CalFred(Fvec,keep,BC,K); % calculate the reduced force vector in a function

u = Kred\Freduced; % solve for unknown displacements
uTotal = zeros(1,length(NodePos)); % create vector of the right size
uTotal(keep) = u; % put calculated displacement on right locations
uTotal(BC(1,:)) = BC(2,:); % insert the applied BC

strain = calcStrain(uTotal,len); % calculate strain in function
stress = calcStress(E,strain); % calculate stress in function

end

function Area = CalcArea(NodePos,L)
    w1 = 50;
	w2 = 25;
	t = 3.125;
    Area = (w1 + (w2-w1)/L*(NodePos(1:end-1)+NodePos(2:end))/2)*t; % area calculated in the middle, entire vector in one go
end

function K = StiffMat(Con,Area,E,len)
    K = zeros(length(len)+1); % create global stiffness matrix
    for a = 1:length(len)
        k = E*Area(a)/len(a); % calculate equivalent stiffness of elements
        Kelem = ElemStiffMat(k);
        % add the element stiffness matrix at the correct location of the elements using the connectivity matrix 
		K(Con(a,1),Con(a,1)) = K(Con(a,1),Con(a,1)) + Kelem(1,1);
		K(Con(a,2),Con(a,2)) = K(Con(a,2),Con(a,2)) + Kelem(2,2);
		K(Con(a,1),Con(a,2)) = K(Con(a,1),Con(a,2)) + Kelem(1,2);
		K(Con(a,2),Con(a,1)) = K(Con(a,2),Con(a,1)) + Kelem(2,1);
    end
end

function Kelem = ElemStiffMat(k)
    Kelem = zeros(2,2); % define element stiffness matrix of correct size
    % on following lines the stiffness elements are placed on the correct
    % location.
    Kelem(1,1) = k;
    Kelem(2,2) = k;
    Kelem(1,2) = -k;
    Kelem(2,1) = -k;
end

function Fred = CalFred(F,keep,BC,K)

    Fred = F(keep); % initially the reduced force vector is equal to the forces applied at these nodes
    nodes = BC(1,:); % the nodes that are fixed have to be added (if displacement is non-zero)
    for i = 1:length(nodes)
        for i2 = 1:length(Fred)
            Fred(i2) = Fred(i2) - BC(2,i) * K(keep(i2),nodes(i)); % subtract what was removed from Left-hand side to the force vector
        end
    end
end

function strain = calcStrain(u,len)
    strain = (u(2:end)-u(1:end-1))./len;
end

function stress = calcStress(E,strain)
    stress = E*strain;
end