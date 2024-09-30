import numpy as np


def calcStress(strain, E):  # calculate stress for each element
	stress = strain * E
	return stress


def calcStrain(con, nodepos, numel, U, area, e):
	strain = np.zeros(numel)
	F_dummy = np.zeros(numel)
	for i in range(numel):
		node1, node2 = con[i] - 1  # Adjust for zero-indexing
		x1, y1 = nodepos[node1]
		x2, y2 = nodepos[node2]

		L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
		cos_theta = (x2 - x1) / L
		sin_theta = (y2 - y1) / L

		u = np.array([U[2 * node1], U[2 * node1 + 1], U[2 * node2], U[2 * node2 + 1]])

		# Internal force in the element
		F_dummy[i] = (e * area[i] / L) * np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]) @ u

		# Strain in the element
		strain[i] = (1 / L) * np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]) @ u

	return strain


def calcUSSR(E, BC, F, Con, NodePos, Area):
	F = F.astype('float64')
	num_nodes = NodePos.shape[0]
	num_elements = Con.shape[0]
	dof_per_node = 2
	u = np.zeros(num_nodes)
	total_dof = num_nodes * dof_per_node
	K_global = np.zeros((total_dof, total_dof))

	for i in range(num_elements):
		node1, node2 = Con[i] - 1
		x1, y1 = NodePos[node1]
		x2, y2 = NodePos[node2]

		L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
		cos_theta = (x2 - x1) / L
		sin_theta = (y2 - y1) / L

		k = (E * Area[i] / L) * np.array([[1, -1], [-1, 1]])
		T = np.array([[cos_theta, sin_theta, 0, 0], [0, 0, cos_theta, sin_theta]])
		k_global = T.T @ k @ T

		dof_indices = np.array([2 * node1, 2 * node1 + 1, 2 * node2, 2 * node2 + 1])

		for a in range(4):
			for b in range(4):
				K_global[dof_indices[a], dof_indices[b]] += k_global[a, b]
	K_global_nobc = np.copy(K_global)
	for i, dof in enumerate(BC[0]-1):
		displacement = BC[1][i]
		F -= K_global[:, dof] * displacement
		F[dof] = displacement
		K_global[dof, :] = 0
		K_global[:, dof] = 0
		K_global[dof, dof] = 1

	U = np.linalg.solve(K_global, F)

	F_int = np.zeros(num_elements)
	for i in range(num_elements):
		node1, node2 = Con[i] - 1
		x1, y1 = NodePos[node1]
		x2, y2 = NodePos[node2]

		L = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
		cos_theta = (x2 - x1) / L
		sin_theta = (y2 - y1) / L

		u = np.array([U[2 * node1], U[2 * node1 + 1], U[2 * node2], U[2 * node2 + 1]])
		F_int[i] = (E * Area[i] / L) * np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]) @ u

	strain = calcStrain(Con, NodePos, num_elements, U, Area, E)
	stress = calcStress(strain, E)

	reactions = K_global_nobc @ U - F

	return U, strain, stress, reactions
