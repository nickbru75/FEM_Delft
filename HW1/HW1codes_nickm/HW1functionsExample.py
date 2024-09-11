import numpy as np


class Shape:
	def __init__(self, length: float, n: int, w_1: float, w_end: float, t: float):
		self.length = length
		self.n = n
		self.w_1 = w_1
		self.w_end = w_end
		self.t = t
		self.sections = self.compute_sections()

	def compute_sections(self):
		l = self.length / self.n
		h = np.array([self.w_1 - (self.w_1 - self.w_end)/self.n * (i+1) for i in range(self.n)])
		areas = h * self.t
		return areas


def calc_equivalent_stifness(sections: np.array, E: float, L: float):
	return sections * E / L


def calcStrain(E: float, stress: np.array):
	return stress / E


def calcStress(F: float, A: np.array):
	return F / A


def apply_boundary_conditions(K, F, BC):
	# Convert node numbers to zero-based indices
	node_indices = BC[0] - 1  # Convert to zero-based indexing
	displacements = BC[1]  # Corresponding displacements

	# Make copies of the original stiffness matrix and load vector
	K_mod = np.copy(K)
	F_mod = np.copy(F)

	# Process each boundary condition
	for idx, displacement in sorted(zip(node_indices, displacements), key=lambda x: x[0], reverse=True):
		if displacement != 0:  # Non-zero displacement boundary condition
			# Subtract the column times the applied displacement from the load vector
			F_mod -= K_mod[:, idx] * displacement

		# Remove the row from the stiffness matrix and the corresponding entry from the load vector
		K_mod = np.delete(K_mod, idx, axis=0)
		F_mod = np.delete(F_mod, idx)

		# Remove the column from the stiffness matrix
		K_mod = np.delete(K_mod, idx, axis=1)

	return K_mod, F_mod


def calcDisp(L, Nelem, F, E, BC):
	bar = Shape(L, Nelem, 50, 25, 3.125)
	k = calc_equivalent_stifness(bar.sections, E, L)
	print(bar.sections, k)
	stress = calcStress(F, bar.sections)
	strain = calcStrain(E, stress)
	# Stiffness matrix
	K = np.zeros((Nelem+1, Nelem+1))
	for i in range(Nelem):
		K[i, i] += k[i]
		K[i, i+1] -= k[i]
		K[i+1, i] -= k[i]
		K[i+1, i+1] += k[i]
	# Load vector
	load = np.zeros(Nelem+1)
	load[0] = -F
	load[-1] = F
	# Boundary conditions
	nodes_modified = [BC[0][i] for i in range(len(BC[0]))]
	K, load = apply_boundary_conditions(K, load, BC)
	u_partial = np.linalg.solve(K, load)
	u = np.zeros(Nelem)
	c_mod = 0
	c_not_mod = 0
	for i in range(Nelem):
		if i+1 in nodes_modified:
			u[i] = BC[1][c_mod]
			c_mod += 1
		else:
			u[i] = u_partial[c_not_mod]
			c_not_mod += 1

	return u*1e6, strain, stress

