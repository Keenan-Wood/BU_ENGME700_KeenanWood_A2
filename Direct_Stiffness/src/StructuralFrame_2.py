import numpy as np
import math_utils as mu

def load_frame(nodes: np.array, elements: list, xsections_list: list, constraint_list: list = [[]], applied_forces: list = [[]]):
    N_nodes = len(nodes)
    N_elements = len(elements)

    # Build constraint array
    constrained = np.zeros((N_nodes, 6))
    for constraint in constraint_list:
        constrained[constraint[0]] = constraint[1:7]

    # Build force array
    forces = np.zeros((N_nodes, 6))
    for force in applied_forces:
        forces[force[0], :] = force[1:7]
    
    # Fill z_vec_list if not provided
    for i, element in enumerate(elements):
        if len(element[3]) == 0:
            node_pair = np.vstack((nodes[element[0], :], nodes[element[1], :]))
            elements[i][3] = get_assumed_z_vec(node_pair)

    # Assemble frame stiffness matrices
    Ke = np.zeros((6*N_nodes, 6*N_nodes))
    for element in elements:
        node_pair = np.vstack((nodes[element[0], :], nodes[element[1], :]))
        local_Ke = calc_local_stiffness(node_pair, xsections_list[element[2]])

        # Compare with provided utility function
        (E, A, Iy, Iz, Ip, J, v) = tuple(xsections_list[element[2]])
        L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
        check_Ke = mu.local_elastic_stiffness_matrix_3D_beam(E, v, A, L, Iy, Iz, J)
        Ke_dif = local_Ke - check_Ke
        if np.max(abs(Ke_dif)) > 10**-10: raise Exception("Incorrect Local Stiffness Matrix")
        
        gam_full = calc_coord_transform(node_pair, element[3])
        global_Ke = np.matmul(np.transpose(gam_full), np.matmul(local_Ke, gam_full))
        a_rng = np.arange(6*element[0], 6*(element[0]+1))
        b_rng = np.arange(6*element[1], 6*(element[1]+1))
        Ke[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += global_Ke

        force_vec = np.vstack((forces[element[0], :], forces[element[1], :]))
        Kg = calc_geometric_local_stiffness(node_pair, xsections_list[element[2]], force_vec)

        # Compare with provided utility function
        L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
        (_, A, _, _, Ip, _, _) = tuple(xsections_list[element[2]])
        (_, _, _, _, My1, Mz1) = tuple(force_vec[0,:])
        (Fx2, _, _, Mx2, My2, Mz2) = tuple(force_vec[1,:])
        check_Kg = mu.local_geometric_stiffness_matrix_3D_beam(L, A, Ip, Fx2, Mx2, My1, Mz1, My2, Mz2)
        Kg_dif = Kg - check_Kg
        if np.max(abs(Kg_dif)) > 10**-10: raise Exception("Incorrect Local Geometric Stiffness Matrix")

    # Partition frame stiffness matrix
    free_ind = np.flatnonzero(constrained == 0)
    fixed_ind = np.flatnonzero(constrained)
    Ke_ff = Ke[np.ix_(free_ind, free_ind)]
    Ke_sf = Ke[np.ix_(fixed_ind, free_ind)]

    # Calculate displacements and forces
    all_disps = np.zeros(6*N_nodes)
    all_forces = forces.reshape(-1)
    free_disps = np.matmul(np.linalg.inv(Ke_ff), all_forces[free_ind])
    support_forces = np.matmul(Ke_sf, free_disps)
    all_disps[free_ind] = free_disps
    all_forces[fixed_ind] = support_forces
    all_disps = np.reshape(all_disps, (N_nodes, 6))
    all_forces = np.reshape(all_forces, (N_nodes, 6))

    return (all_disps, all_forces)

def get_assumed_z_vec(node_pair: np.array):
    x_vec = node_pair[1, 0:3] - node_pair[0, 0:3]
    z_vec = np.cross(x_vec, np.array([1, 1, 0]))
    if np.linalg.norm(z_vec) / np.linalg.norm(x_vec) < 0.01:
        z_vec = np.cross(x_vec, np.array([1, -1, 0]))
    z_vec = z_vec / np.linalg.norm(z_vec)
    return z_vec

def calc_local_stiffness(node_pair: np.array, xsec: list):
    L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
    (E, A, Iy, Iz, Ip, J, v) = tuple(xsec)
    k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz]) / A
    k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz]) / A
    k_cross_1 = 6/L * np.array([Iz, -Iy, 0, -Iy, Iz]) / A
    k_cross_2 = 6/L * np.array([Iz, -Iy, 0, Iy, -Iz]) / A
    k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
    k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
    k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_2, -1))
    k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_1, -1))
    K_local = E*A/L * np.vstack((np.hstack((k_a, k_b)), np.hstack((k_c, k_d))))
    return K_local

def calc_geometric_local_stiffness(node_pair: np.array, xsec, forces: np.array):
    L = np.linalg.norm(node_pair[1, :] - node_pair[0, :])
    (_, A, _, _, Ip, _, _) = tuple(xsec)
    (_, _, _, _, My1, Mz1) = tuple(forces[0,:])
    (Fx2, _, _, Mx2, My2, Mz2) = tuple(forces[1,:])
    k_diag_1 = Fx2/L * np.array([1, 6/5, 6/5, Ip/A, 2/15*L**2, 2/15*L**2])
    k_diag_2 = -k_diag_1 + Fx2 * np.array([0, 0, 0, 0, L/10, L/10])

    k_a = np.zeros((6, 6))
    k_a[1, 3:6] = 1/L * np.array([My1, Mx2, Fx2*L/10])
    k_a[2, 3:6] = 1/L * np.array([Mz1, -Fx2*L/10, Mx2])
    k_a[3, 4:6] = np.array([-(2*Mz1-Mz2)/6, (2*My1-My2)/6])

    k_c = np.zeros((6, 6))
    k_c[1, 3:6] = 1/L * np.array([-My2, Mx2, -Fx2*L/10])
    k_c[2, 3:6] = 1/L * np.array([-Mz2, Fx2*L/10, Mx2])
    k_c[3, 4:6] = np.array([-(2*Mz2-Mz1)/6, (2*My2-My1)/6])

    k_b = -k_c - np.transpose(k_a) + np.diag(k_diag_2)
    k_b[3, 4:6] = np.array([-(Mz1+Mz2)/6, (My1+My2)/6])
    k_b[4, 5] = Mx2/2
    k_b[4:6, 3] = k_b[3, 4:6]
    k_b[5, 4] = -Mx2/2

    Kg = np.zeros((12, 12))
    Kg[0:6, 0:6] += k_a
    Kg[0:6, 6:12] += k_b
    Kg[6:12, 6:12] += k_c
    Kg += np.transpose(Kg) + np.diag(np.hstack((k_diag_1, k_diag_1)))
    return Kg

def calc_coord_transform(node_pair: np.array, z_vec: np.array):
    x_vec = node_pair[1, 0:3] - node_pair[0, 0:3]
    y_vec = np.cross(z_vec, x_vec)
    x_vec = x_vec / np.linalg.norm(x_vec)
    y_vec = y_vec / np.linalg.norm(y_vec)
    z_vec = z_vec / np.linalg.norm(z_vec)
    gam_small = np.vstack((x_vec, y_vec, z_vec))

    # Check gam_small with provided utility function
    (x1, y1, z1) = tuple(node_pair[0, 0:3])
    (x2, y2, z2) = tuple(node_pair[1, 0:3])
    v_temp = z_vec
    check_gam = mu.rotation_matrix_3D(x1, y1, z1, x2, y2, z2, v_temp)
    gam_dif = gam_small - check_gam
    if np.max(abs(gam_dif)) > 10**-10: raise Exception("Incorrect Transform Matrix")

    gam_full = np.zeros((12, 12))
    for i in range(0,4):
        gam_full[3*i:3*i+3, 3*i:3*i+3] = gam_small

    # Check gam_full with provided utility function
    gam_full_dif = gam_full - mu.transformation_matrix_3D(check_gam)
    if np.max(abs(gam_full_dif)) > 10**-10: raise Exception("Incorrect Full Transform Matrix")

    return gam_full