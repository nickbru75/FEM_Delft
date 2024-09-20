
import numpy as np

def StiffMat(Con,Area,E,length):    # define stiffness matrix
	K = np.zeros([len(length)+1,len(length)+1]) # global stiffness matrix of correct size
	a = 0
	while a < len(length): # loop over the elements 
		k = E*Area[a]/length[a] # calculate equivalent stiffness 
		Kelem = ElemStiffMat(k)
		# on following lines the parts of the element stiffness matrix are placed in the global stiffness matrix
		K[Con[a,0]-1,Con[a,0]-1] = K[Con[a,0]-1,Con[a,0]-1] + Kelem[0,0]
		K[Con[a,1]-1,Con[a,1]-1] = K[Con[a,1]-1,Con[a,1]-1] + Kelem[1,1]
		K[Con[a,0]-1,Con[a,1]-1] = K[Con[a,0]-1,Con[a,1]-1] + Kelem[0,1]
		K[Con[a,1]-1,Con[a,0]-1] = K[Con[a,1]-1,Con[a,0]-1] + Kelem[1,0]
		a = a+1
	return K

def ElemStiffMat(k):
	Kelem = np.zeros([2,2]) # define element stiffness matrix of correct size
    # on following lines the stiffness elements are placed on the correct location
	Kelem[0,0] = k
	Kelem[1,1] = k
	Kelem[0,1] = -k
	Kelem[1,0] = -k
	return Kelem
	
def ApplyBC(F,BC):
	keep = np.arange(len(F)) # initially all nodes are to be kept 
	nodes = BC[0,:] # fixed nodes are defined in first part of BC
	ind =np.zeros(len(nodes))
	a = 0
	while a < len(nodes):
		ind[a] = np.where(keep == nodes[a]-1)[0][0] # check whether node has a boundary condition applied or not 
		a = a + 1
	keep = (set(keep) - set(ind))
	keep = list(keep)
	
	return keep 
	
def CalFred(F,keep,BC,K):
	Freduced = F[keep] # initially the reduced force vector is equal to the forces applied at these nodes
	nodes = BC[0,:] # the nodes that are fixed have to be added (if displacement is non-zero)
	nodes = nodes.astype(int) # make sure it are all integers
	nodebc = 0
	while nodebc < len(nodes):
		a = 0
		while a < len(Freduced):
			Freduced[a] = Freduced[a] - BC[1,nodebc] * K[keep[a],nodes[nodebc]-1] # subtract what was removed from Left-hand side to the force vector
			
			a = a + 1
		nodebc = nodebc + 1
	
	return Freduced
	
def calcStrain(u,length): # calculate strain for each element 
	strain = np.zeros(len(length)) # calculate vector of right size 
	a = 0
	while a < len(strain):
		strain[a] = (u[a+1]-u[a])/length[a] # calculate strain of an element 
		a = a + 1
	return strain
	
def calcStress(strain, E): # calculate stress for each element 
	stress = strain * E
	return stress 

def DefineConn(Nelem): # define connectivity matrix 
	Con = np.zeros([Nelem,2]) # create vector of right size
	a = 0
	while a < Nelem: # fill up the connectivity matrix for this specific case 
		Con[a,0] = a+1 
		Con[a,1] = a+2
		a = a + 1
	Con = Con.astype(int) # make sure it are all integers
	return Con

def DefForceVec(Nelem,F): # define force vector (force on last node)
	Fvec = np.zeros(Nelem+1) # vector of right size, all zeros 
	Fvec[-1]=F # put the force on the final node 
	return Fvec

def CalcNodePos(Nelem,L):
	NodePos = np.zeros(Nelem+1) # create vector of correct size 
	a = 0
	while a < Nelem:
		NodePos[a+1] = NodePos[a] + L/Nelem # fill up the node positions
		a = a + 1
	return NodePos 

def CalcArea(NodePos,L): # function to calculate the area 
	Nelem = len(NodePos) - 1
	Area = np.zeros(Nelem) # create vector of right size
	w1 = 50
	w2 = 25
	t = 3.125
	a = 0
	while a < Nelem:
		Area[a] = (w1 + (w2-w1)/L*(NodePos[a]+NodePos[a+1])/2)*t #area calculated in centre of element 
		a = a + 1
	return Area
	
def CalcLength(L,Nelem):
	Length = np.zeros(Nelem) # create vector of right size 
	a = 0
	while a < Nelem:
		Length[a] = L/Nelem # put length of element at right position 
		a = a + 1
	return Length

def calcDisp(L, Nelem, F, E, BC): # main function called from other program 
	F = DefForceVec(Nelem,F) # calculate force vector in function
	Con = DefineConn(Nelem) # calculate connectivity matrix in function
	NodePos = CalcNodePos(Nelem,L) # calculate Node position vector in function
	Area = CalcArea(NodePos,L) # calculate area in function
	length = CalcLength(L,Nelem) # calculate length in function
	
	K = StiffMat(Con,Area,E,length) # calculate global stiffness matrix  in function
	
	keep = ApplyBC(F,BC) #determine which indices to keep 
	
	KredInt = K[keep,:] # only keep the required rows
	Kreduced = KredInt[:,keep] # only keep required columns 
	
	Freduced = CalFred(F,keep,BC,K) # calculate reduced force vector 
	
	u = np.linalg.solve(Kreduced, Freduced) # solve for unknown displacements
	
	utotal = np.zeros(len(F)) # vector of right size for displacement 
	utotal[keep] = u # put calculated displacement on right locations
	utotal[BC[0,:]-1] = BC[1,:] # insert the applied BC
	
	strain = calcStrain(utotal,length) # calculate strain in function
	stress = calcStress(strain, E) # calculate stress in function
	
	return utotal, strain, stress 