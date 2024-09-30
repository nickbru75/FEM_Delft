import numpy as np


def calcStress(strain, E):  # Calculate stress for each element
    return np.multiply(strain, E)  # Alternative to using strain * E


def calcStrain(connectivity, node_positions, num_elements, displacements, areas, E_modulus):
    strain = np.zeros(num_elements)
    internal_forces = np.zeros(num_elements)

    for i in range(num_elements):
        # Get node indices for the current element
        node1, node2 = connectivity[i] - 1
        # Get coordinates of the nodes
        x1, y1 = node_positions[node1]
        x2, y2 = node_positions[node2]

        # Calculate the length of the element
        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        # Calculate direction cosines
        cos_theta = (x2 - x1) / length
        sin_theta = (y2 - y1) / length

        # Displacement vector for the current element
        u = np.array([displacements[2 * node1], displacements[2 * node1 + 1],
                      displacements[2 * node2], displacements[2 * node2 + 1]])

        # Calculate internal axial force for the element
        internal_forces[i] = (E_modulus * areas[i] / length) * np.matmul(np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]), u)
        # Compute strain for the element
        strain[i] = np.matmul(np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]), u) / length

    return strain


def calcUSSR(E_modulus, boundary_conditions, external_forces, connectivity, node_positions, areas):
    # Convert force array to float type
    external_forces = external_forces.astype('float64')
    num_nodes = node_positions.shape[0]  # Total number of nodes
    num_elements = connectivity.shape[0]  # Total number of elements
    dof_per_node = 2  # Degrees of freedom per node (2D problem)
    total_dof = num_nodes * dof_per_node  # Total degrees of freedom in the system

    # Initialize global stiffness matrix
    global_stiffness_matrix = np.zeros((total_dof, total_dof))

    # Assemble the global stiffness matrix
    for i in range(num_elements):
        # Get node indices for the current element
        node1, node2 = connectivity[i] - 1
        # Get coordinates of the nodes
        x1, y1 = node_positions[node1]
        x2, y2 = node_positions[node2]

        # Calculate the length of the element
        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        # Calculate direction cosines
        cos_theta = (x2 - x1) / length
        sin_theta = (y2 - y1) / length

        # Local stiffness matrix for the element
        local_stiffness = (E_modulus * areas[i] / length) * np.array([[1, -1], [-1, 1]])
        # Transformation matrix for the element
        transformation_matrix = np.array([[cos_theta, sin_theta, 0, 0], [0, 0, cos_theta, sin_theta]])
        # Global stiffness matrix for the element
        element_global_stiffness = np.matmul(np.matmul(transformation_matrix.T, local_stiffness), transformation_matrix)

        # Degrees of freedom indices for the current element
        dof_indices = np.array([2 * node1, 2 * node1 + 1, 2 * node2, 2 * node2 + 1])

        # Assemble global stiffness matrix
        for j in range(4):
            for k in range(4):
                global_stiffness_matrix[dof_indices[j], dof_indices[k]] += element_global_stiffness[j, k]

    # Make a copy of the global stiffness matrix before applying boundary conditions
    global_stiffness_matrix_nobc = np.copy(global_stiffness_matrix)

    # Apply boundary conditions
    for i, dof in enumerate(boundary_conditions[0] - 1):
        displacement = boundary_conditions[1][i]
        # Adjust the external force vector
        external_forces -= global_stiffness_matrix[:, dof] * displacement
        # Set the force at the constrained degree of freedom to the displacement value
        external_forces[dof] = displacement
        # Modify the global stiffness matrix to account for the boundary condition
        global_stiffness_matrix[dof, :] = 0
        global_stiffness_matrix[:, dof] = 0
        global_stiffness_matrix[dof, dof] = 1

    # Solve for displacements
    displacements = np.linalg.solve(global_stiffness_matrix, external_forces)

    # Initialize internal force array
    internal_forces = np.zeros(num_elements)

    # Calculate internal forces for each element
    for i in range(num_elements):
        # Get node indices for the current element
        node1, node2 = connectivity[i] - 1
        # Get coordinates of the nodes
        x1, y1 = node_positions[node1]
        x2, y2 = node_positions[node2]

        # Calculate the length of the element
        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        # Calculate direction cosines
        cos_theta = (x2 - x1) / length
        sin_theta = (y2 - y1) / length

        # Displacement vector for the current element
        u = np.array([displacements[2 * node1], displacements[2 * node1 + 1],
                      displacements[2 * node2], displacements[2 * node2 + 1]])

        # Calculate internal axial force for the element
        internal_forces[i] = (E_modulus * areas[i] / length) * np.matmul(np.array([-cos_theta, -sin_theta, cos_theta, sin_theta]), u)

    # Calculate strain in each element
    strain = calcStrain(connectivity, node_positions, num_elements, displacements, areas, E_modulus)
    # Calculate stress in each element
    stress = calcStress(strain, E_modulus)
    # Calculate reaction forces
    reactions = np.matmul(global_stiffness_matrix_nobc, displacements) - external_forces

    return displacements, strain, stress, reactions
