import numpy as np


# Function to calculate strain for each element
def calcStrain(u, length):
    # Using a for loop instead of while and using NumPy for vectorized calculation
    strain = np.zeros(len(length))  # Initialize strain array with zeros
    for idx in range(len(strain)):
        strain[idx] = (u[idx + 1] - u[idx]) / length[idx]  # Vectorized calculation
    return strain


# Function to calculate stress for each element using strain and Young's Modulus
def calcStress(strain, E):
    return strain * E  # Direct element-wise multiplication


# Function to calculate strain based on nodal displacements and positions
def calcStrain(connections, node_positions, num_elements, displacements, cross_sectional_area, modulus_of_elasticity):
    strain = np.zeros(num_elements)  # Initialize the strain array
    internal_force_component = np.zeros(num_elements)  # Array for internal force components

    # Precompute node positions and displacement indices
    node_positions_diff = node_positions[connections[:, 1] - 1] - node_positions[connections[:, 0] - 1]
    lengths = np.linalg.norm(node_positions_diff, axis=1)
    cosine_angle = node_positions_diff[:, 0] / lengths
    sine_angle = node_positions_diff[:, 1] / lengths

    # Loop over elements and calculate strain and internal force
    for element in range(num_elements):
        node1, node2 = connections[element] - 1

        # Fetch global displacements for these nodes
        local_displacements = np.array([
            displacements[2 * node1], displacements[2 * node1 + 1],
            displacements[2 * node2], displacements[2 * node2 + 1]
        ])

        # Create the transformation vector
        transformation_vector = np.array([-cosine_angle[element], -sine_angle[element],
                                          cosine_angle[element], sine_angle[element]])

        # Calculate internal force and strain
        internal_force_component[element] = (modulus_of_elasticity * cross_sectional_area[element] / lengths[element]) * \
                                            (transformation_vector @ local_displacements)

        strain[element] = (1 / lengths[element]) * (transformation_vector @ local_displacements)

    return strain


# Function to calculate displacement, strain, stress, and reactions
def calcUSSR(E, BC, F, Con, NodePos, Area):
    num_nodes = NodePos.shape[0]  # Total number of nodes
    num_elements = Con.shape[0]  # Total number of elements
    degrees_of_freedom_per_node = 2  # Typically, 2 DOFs per node (x and y directions)
    total_dof = num_nodes * degrees_of_freedom_per_node  # Total degrees of freedom
    displacements = np.zeros(total_dof)  # Adjust this to store total DOFs

    global_stiffness = np.zeros((total_dof, total_dof))  # Global stiffness matrix

    # Construct the global stiffness matrix
    for element in range(num_elements):
        node1, node2 = Con[element] - 1
        x1, y1 = NodePos[node1]
        x2, y2 = NodePos[node2]

        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)  # Length of the element
        cosine_angle = (x2 - x1) / length
        sine_angle = (y2 - y1) / length

        # Element stiffness matrix in local coordinates
        local_stiffness = (E * Area[element] / length) * np.array([[1, -1], [-1, 1]])

        # Use matrix multiplication with transformation matrix directly
        transformation_matrix = np.array([
            [cosine_angle, sine_angle, 0, 0],
            [0, 0, cosine_angle, sine_angle]
        ])

        global_stiffness_contrib = transformation_matrix.T @ local_stiffness @ transformation_matrix

        dof_indices = np.array([2 * node1, 2 * node1 + 1, 2 * node2, 2 * node2 + 1])

        # Directly update the global stiffness matrix
        np.add.at(global_stiffness, (dof_indices[:, None], dof_indices), global_stiffness_contrib)

    # Apply boundary conditions
    stiffness_with_bc = np.copy(global_stiffness)

    # Use np.ix_ for efficient matrix indexing and setting boundary conditions
    constrained_dof = BC[0] - 1
    constrained_displacements = BC[1]
    free_dof = np.setdiff1d(np.arange(total_dof), constrained_dof)

    F -= global_stiffness[:, constrained_dof] @ constrained_displacements
    F[constrained_dof] = constrained_displacements
    global_stiffness[constrained_dof, :] = 0
    global_stiffness[:, constrained_dof] = 0
    global_stiffness[constrained_dof, constrained_dof] = 1

    # Solve the linear system only for free degrees of freedom
    final_displacements = np.zeros_like(F)
    final_displacements[free_dof] = np.linalg.solve(global_stiffness[np.ix_(free_dof, free_dof)], F[free_dof])

    # Calculate internal forces for each element using a vectorized approach
    internal_forces = np.zeros(num_elements)
    for element in range(num_elements):
        node1, node2 = Con[element] - 1
        x1, y1 = NodePos[node1]
        x2, y2 = NodePos[node2]

        length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        cosine_angle = (x2 - x1) / length
        sine_angle = (y2 - y1) / length

        element_displacements = np.array([
            final_displacements[2 * node1], final_displacements[2 * node1 + 1],
            final_displacements[2 * node2], final_displacements[2 * node2 + 1]
        ])

        transformation_vector = np.array([-cosine_angle, -sine_angle, cosine_angle, sine_angle])
        internal_forces[element] = (E * Area[element] / length) * (transformation_vector @ element_displacements)

    # Calculate strain and stress for the system
    element_strain = calcStrain(Con, NodePos, num_elements, final_displacements, Area, E)
    element_stress = calcStress(element_strain, E)

    # Calculate reactions at the nodes
    reactions = stiffness_with_bc @ final_displacements - F

    return final_displacements, element_strain, element_stress, reactions