import HW1functionsExample as HW1f
import numpy as np

# input
L = 250 							# total length of the bar (scalar value)
Nelem = 4	    				# number of elements (scalar value)
F = 1e3  						# applied force on the final nodes (scalar value)
E = 70e9	 					# elastic modulus (scalar value)
BC = np.array([[1], [0]]) 		# form: [node numbers], [applied displacements]; each vector can have multiple inputs
# example BC: np.array([[i],[x]]) should apply a displacement of x on node i
# example 2 BC: np.array([[i, j],[x, y]]) should apply a displacement of x on node i and a displacement of y on onde j
# additional note: node numbering starts at 1, not at 0 for these inputs.

# call functions file 
u, strain, stress = HW1f.calcDisp(L, Nelem, F, E, BC)

# output
for el in range(Nelem):
	print("the stress in element ", (el+1), " equals ", "%.5f" % stress[el], "MPa")
	print("the displacement in element ", (el+1), " equals ", "%.5f" % u[el], "mm")
print("The displacement of the final node is ", "%.5f" % u[-1], "mm")

test = [3.828125e+10, 3.281250e+10, 2.734375e+10, 2.187500e+10]
coeff = test[-1]/8.75
for i in range(len(test)):
	test[i] = test[i]/coeff
print(test)