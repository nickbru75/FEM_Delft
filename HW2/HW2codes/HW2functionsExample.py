import numpy as np

def calcStrain(Nelem):
	strain = np.ones(Nelem)
	return strain 
	
def calcStress(Nelem):
	stress = 2*np.ones(Nelem)
	return stress

def calcUSSR(E,BC,F,Con,NodePos,Area):
	
	u = np.zeros(2*len(NodePos))
	
	strain = calcStrain(len(Con)) # just an example of a function, input is also just to have an input
	stress = calcStress(len(Con)) # just an example of a function, input is also just to have an input
	
	reactions = -F # just a placeholder to have a vector of the right size
	
	return u, strain, stress, reactions

