import HW1functions as HW1f
import numpy as np

# input
L =  250 # total length of the bar
Nelem = 4 # number of elements
F =  1000# applied force on the final nodes
E =  70e3 # elastic modulus 
BC = np.array([[1],[0]]) #form: [node numbers], [applied displacements]

# functions here
u, strain, stress = HW1f.calcDisp(L, Nelem, F, E, BC)

# output
for el in range(Nelem):
	print ("the stress in element ", (el+1) , " equals ", "%.2f" % stress[el], "MPa")
	
print("The displacement of the final node is ", "%.4f" % u[-1], "mm")