import HW1functionsExample as HW1f
import numpy as np

# input
L =  1 # total length of the bar (scalar value)
Nelem = 1 # number of elements (scalar value)
F = 1  # applied force on the final nodes (scalar value)
E = 1 # elastic modulus (scalar value)
BC = np.array([[1 ],[1 ]]) #form: [node numbers], [applied displacements]; each vector can have multiple inputs
# example BC: np.array([[i],[x]]) should apply a displacement of x on node i
# example 2 BC: np.array([[i, j],[x, y]]) should apply a displacement of x on node i and a displacement of y on onde j
# additional note: node numbering starts at 1, not at 0 for these inputs.

# call functions file 
u, strain, stress = HW1f.calcDisp(L, Nelem, F, E, BC)

# output
for el in range(Nelem):
	print ("the stress in element ", (el+1) , " equals ", "%.2f" % stress[el], "MPa")
	
print("The displacement of the final node is ", "%.4f" % u[-1], "mm")