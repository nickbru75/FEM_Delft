import HW2functions_5943183 as HW2f
import numpy as np

# input

E = 70e3 # E-modulus
BC = np.array([[1, 2, 4, 5, 6, 8], [0, 0, 0, 0, 0, 0]]) #boundary conditions; form: [degrees of freedom], [applied displacements]
F = np.array([0., 0., 1500., 0., 0., 0., -500., 0., 0., 0.]) # force vector
Con = np.array([[1, 2], [5, 3], [4, 5], [5, 1]]) # connectivity matrix; form: [[element 1], [element 2],etc]
NodePos = np.array([[1250, 750], [2000, 0], [750, 0], [0, 0], [500, 500]]) # node position; form: [[element 1],[element 2],etc]
Area = np.array([50, 30, 30, 50]) # vector with the area for each element

# functions here
u, strain, stress, reactions = HW2f.calcUSSR(E,BC,F,Con,NodePos,Area)

# output: do not change!
for el in range(len(Con)):
	print ("the stress in element ", (el+1) , " equals ", "%.2f" % stress[el], "MPa")

for node in range(len(NodePos)):
	print("The displacement of node", node+1 , " in x-direction is ", "%.3f" % u[2*(node+1)-2], "mm")
	print("The displacement of node", node+1 , " in y-direction is ", "%.3f" % u[2*(node+1)-1], "mm")
	print("The reaction force at node", node+1 , " in x-direction is ", "%.0f" % reactions[2*(node+1)-2], "N")
	print("The reaction force at node", node+1 , " in y-direction is ", "%.0f" % reactions[2*(node+1)-1], "N")
	
pass