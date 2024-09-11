import numpy as np

def calcStrain(Nelem):
	strain = np.ones(Nelem)
	return strain 
	
def calcStress(Nelem):
	stress = 2*np.ones(Nelem)
	return stress

def calcDisp(L, Nelem, F, E, BC):
	# this is where the input is tranferred to; the functions defined above are called from line 15 and 16 in this example
	u = np.zeros(Nelem)
	
	strain = calcStrain(Nelem) # just an example of a function, input is also just to have an input
	stress = calcStress(Nelem) # just an example of a function, input is also just to have an input
	
	return u, strain, stress 

