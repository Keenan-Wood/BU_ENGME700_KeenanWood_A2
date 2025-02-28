import numpy as np

# Frame geometry definition
nodes = np.array([[0,0,10,0,0,0], [15,0,10,0,0,0], [15,0,0,0,0,0]])
zvec = np.array([[0,0,1], [1,0,0]])
elements = [[0,1,0,zvec[0,:]], [1,2,0,zvec[1,:]]]

# Cross section list
E = 1000
A = b*h
I_y = h*b**3/12
I_z = b*h**3/12
J = b*h*(b**2+h**2)/12
v = .3
xsection = [[E, A, I_y, I_z, J, v]]

# Constraint list (node_id, type) - 1 for pinned, 2 for fixed
constraints = [[0,2], [2,1]]

# Force list (node_id, forces on each DOF)
forces = [[1,-0.05,0.075,0.1,-0.05,0.1,-0.25]]


def load_frame(nodes: np.array, elements: list, xsections_list: list, constraint_list: list, force_list: list):
    N_nodes = len(nodes)
    N_elements = len(elements)

    # Build element array
    #elements = zeros((N_nodes, N_nodes))
    #for node_pair in element_list:
    #    elements[node_pair[0], node_pair[1]] = node_pair[3]

    # Build constraint array
    constrained = np.zeros((N_nodes, 6))
    for constraint in constraint_list:
        if constraint[1] == 1: constrained[constraint[0], range(3,6)] = 1
        elif constraint[1] == 2: constrained[constraint[0], range(0,6)] = 1

    # Build force array
    forces = np.zeros((N_nodes, 6))
    for force in force_list:
        forces[force[0]] = force[1:6]
    
    # Fill z_vec_list if not provided
    for i_el in range(N_elements):
        if not elements[i_el][3]:
            node_pair = np.vstack(nodes[elements[i_el][0]], nodes[elements[i_el][1]])
            elements[i_el][3] = get_assumed_z_vec(node_pair)

####################
    # Assemble global stiffness matrix
    frame_K_e = np.zeros(len(self.nodes), len(self.nodes))
    for elem in elements:
        frame_K_e = overlay_array(frame_K_e, elem_K_e, elem.node_a.id, elem.node_b.id)
##################################


def get_assumed_z_vec(node_pair: np.array):
    x_vec = node_pair[1, :] - node_pair[0, :]
    z_vec = np.cross(x_vec, np.array([1, 1, 0]))
    if np.linalg.norm(z_vec) / np.linalg.norm(x_vec) < 0.01:
        z_vec = np.cross(x_vec, np.array([1, -1, 0]))
    z_vec = z_vec / np.linalg.norm(z_vec)
    return z_vec

def calc_coord_transform(node_pair: np.array, z_vec: np.array):
    x_vec = node_pair[1, :] - node_pair[0, :]
    y_vec = np.cross(z_vec, x_vec)
    gam_small = np.hstack(x_vec, y_vec, z_vec)
    gam_full = np.zeros(12, 12)
    for i in range(0,4):
        gam_full[3*i:3*i+2, 3*i:3*i+2] = gam_small
    return gam_full

def calc_local_stiffness(node_pair: np.array, xsec: list):
    L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
    (E, A, Iy, Iz, J, v) = tuple(xsec)
    k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz]) / A
    k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz]) / A
    k_cross_1 = 6/L**2 * np.array([Iz, -Iy, 0, -Iy, Iz]) / A
    k_cross_2 = 6/L**2 * np.array([Iz, -Iy, 0, Iy, -Iz]) / A
    k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
    k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
    k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_1, -1))
    k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_2, -1))
    K_local = E*A/L * np.vstack(np.hstack(k_a, k_b), np.hstack(k_c, k_d))
    return K_local

def overlay_array(ar1: np.array, ar2: np.array, a: int, b: int):
    (a, b) = handle_inputs_overlay_array(ar1, ar2, a, b)
    w2 = len(ar2)/2
    a_rng = range(a*w2, (a+1)*w2)
    b_rng = range(b*w2, (b+1)*w2)
    ar1[a_rng, a_rng] = ar1[a_rng, a_rng] + ar2[0:w2, 0:w2]
    ar1[a_rng, b_rng] = ar1[a_rng, b_rng] + ar2[0:w2, w2:-1]
    ar1[b_rng, a_rng] = ar1[b_rng, a_rng] + ar2[w2:-1, 0:w2]
    ar1[b_rng, b_rng] = ar1[b_rng, b_rng] + ar2[w2:-1, w2:-1]
    return ar1
    
