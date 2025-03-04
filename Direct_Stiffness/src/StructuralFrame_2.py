import numpy as np
import scipy as sp

def load_frame(nodes: np.array, elements: list, xsections_list: list, constraint_list: list = [[]], applied_forces: list = [[]]):
    N_nodes = len(nodes)
    N_elements = len(elements)

    # Build constraint array
    constrained = np.zeros((N_nodes, 6))
    for constraint in constraint_list: constrained[constraint[0]] = constraint[1:7]

    # Build force array
    forces = np.zeros((N_nodes, 6))
    for force in applied_forces: forces[force[0], :] = force[1:7]
    
    # Fill z_vec_list if not provided
    for i, element in enumerate(elements):
        if len(element[3]) == 0:
            node_pair = np.vstack((nodes[element[0], :], nodes[element[1], :]))
            elements[i][3] = get_assumed_z_vec(node_pair)

    # Assemble frame stiffness matrices
    Ke = np.zeros((6*N_nodes, 6*N_nodes))
    Kg = np.zeros((6*N_nodes, 6*N_nodes))
    for element in elements:
        node_pair = np.vstack((nodes[element[0], :], nodes[element[1], :]))
        local_Ke = calc_local_stiffness(node_pair, xsections_list[element[2]])
        
        gam_full = calc_coord_transform(node_pair, element[3])
        global_Ke = np.matmul(gam_full.T, np.matmul(local_Ke, gam_full))
        a_rng = np.arange(6*element[0], 6*(element[0]+1))
        b_rng = np.arange(6*element[1], 6*(element[1]+1))
        Ke[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += global_Ke

        ### Change!! - Input forces from output of elastic stiffness analysis (transformed into local coords)
        ### - Then solve for lambda with Kg_ff and Ke_ff (scipy.linalg.eig(Ke_ff, -Kg_ff)) - lambda_crit = smallest positive eigenvalue
        ## Part 2: Elastic critical load (eigenvalue); Buckled shape (eigenvector)
        ## Post processing - internal (interpolated) forces, displacements
        force_vec = np.vstack((forces[element[0], :], forces[element[1], :]))
        local_Kg = calc_geometric_local_stiffness(node_pair, xsections_list[element[2]], force_vec)
        global_Kg = np.matmul(gam_full.T, np.matmul(local_Kg, gam_full))
        Kg[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += global_Kg

    # Partition frame stiffness matrix
    free_ind = np.flatnonzero(constrained == 0)
    fixed_ind = np.flatnonzero(constrained)
    Ke_ff = Ke[np.ix_(free_ind, free_ind)]
    Ke_sf = Ke[np.ix_(fixed_ind, free_ind)]
    Kg_ff = Kg[np.ix_(free_ind, free_ind)]
    Kg_sf = Kg[np.ix_(fixed_ind, free_ind)]

    # Calculate displacements and forces
    disps = np.zeros((N_nodes, 6))
    disps_vec = disps.reshape(-1)
    forces_vec = forces.reshape(-1)
    free_disps = np.matmul(np.linalg.inv(Ke_ff), forces_vec[free_ind])
    support_forces = np.matmul(Ke_sf, free_disps)
    disps_vec[free_ind] = free_disps
    forces_vec[fixed_ind] = support_forces

    return (disps, forces)

def get_assumed_z_vec(node_pair: np.array):
    x_vec = node_pair[1, 0:3] - node_pair[0, 0:3]
    z_vec = np.cross(x_vec, np.array([1, 1, 0]))
    if np.linalg.norm(z_vec) / np.linalg.norm(x_vec) < 0.01:
        z_vec = np.cross(x_vec, np.array([1, -1, 0]))
    return z_vec

def calc_local_stiffness(node_pair: np.array, xsec: list):
    L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
    (E, A, Iy, Iz, Ip, J, v) = tuple(xsec)
    k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz])
    k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz])
    k_cross_1 = 6/L * np.array([Iz, -Iy, 0, -Iy, Iz])
    k_cross_2 = 6/L * np.array([Iz, -Iy, 0, Iy, -Iz])
    k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
    k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
    k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_2, -1))
    k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_1, -1))
    K_local = E/L * np.vstack((np.hstack((k_a, k_b)), np.hstack((k_c, k_d))))
    return K_local

def calc_geometric_local_stiffness(node_pair: np.array, xsec, forces: np.array):
    L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
    (_, A, _, _, Ip, _, _) = tuple(xsec)
    (_, _, _, _, My1, Mz1) = tuple(forces[0,:])
    (Fx2, _, _, Mx2, My2, Mz2) = tuple(forces[1,:])
    k_diag_1 = Fx2/L * np.array([1, 6/5, 6/5, Ip/A, 2/15*L**2, 2/15*L**2])
    k_diag_2 = -k_diag_1 + Fx2 * np.array([0, 0, 0, 0, L/10, L/10])

    (k_a, k_d) = (np.zeros((6, 6)), np.zeros((6,6)))
    k_a[np.ix_([1,2,3], [3,4,5])] = 1/L * np.array([[My1, Mx2, Fx2*L/10], [Mz1, -Fx2*L/10, Mx2], [0, -L*(2*Mz1-Mz2)/6, L*(2*My1-My2)/6]])
    k_d[np.ix_([1,2,3], [3,4,5])] = 1/L * np.array([[-My2, Mx2, -Fx2*L/10], [-Mz2, Fx2*L/10, Mx2], [0, -L*(2*Mz2-Mz1)/6, L*(2*My2-My1)/6]])

    k_b = -k_d - np.transpose(k_a) + np.diag(k_diag_2)
    k_b[3, 4:6] = np.array([-(Mz1+Mz2)/6, (My1+My2)/6])
    k_b[4:6, 3] = k_b[3, 4:6]
    k_b[np.array([4,5]), np.array([5,4])] = np.array([Mx2/2, -Mx2/2])

    Kg = np.vstack((np.hstack((k_a, k_b)), np.hstack((np.zeros((6,6)), k_d))))
    Kg += np.transpose(Kg) + np.diag(np.hstack((k_diag_1, k_diag_1)))
    return Kg

def calc_coord_transform(node_pair: np.array, z_vec: np.array):
    x_vec = node_pair[1, 0:3] - node_pair[0, 0:3]
    gam = np.vstack((x_vec, np.cross(z_vec, x_vec), z_vec))
    row_norms = np.linalg.norm(gam, axis=1)
    gam = gam / row_norms[:, np.newaxis]
    gam_full = sp.sparse.block_diag((gam, gam, gam, gam)).toarray()
    return gam_full