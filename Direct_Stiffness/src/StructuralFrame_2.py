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

    # Check assembled stiffness matrix (for example 2)
    Ke_true = true_Ke_ex2()
    Ke_full_dif = Ke - Ke_true
    if np.max(abs(Ke_full_dif)) > 10**-10: raise Exception("Incorrect Assembled Stiffness Matrix")

    # Partition frame stiffness matrix
    free_ind = np.flatnonzero(constrained == 0)
    fixed_ind = np.flatnonzero(constrained)
    Ke_ff = Ke[np.ix_(free_ind, free_ind)]
    Ke_sf = Ke[np.ix_(fixed_ind, free_ind)]
    Kg_ff = Kg[np.ix_(free_ind, free_ind)]
    Kg_sf = Kg[np.ix_(fixed_ind, free_ind)]

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

    # Compare with provided utility function
    check_Ke = mu.local_elastic_stiffness_matrix_3D_beam(E, v, A, L, Iy, Iz, J)
    Ke_dif = K_local - check_Ke
    if np.max(abs(Ke_dif)) > 10**-10: raise Exception("Incorrect Local Stiffness Matrix")

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

    # Compare with provided utility function
    check_Kg = mu.local_geometric_stiffness_matrix_3D_beam(L, A, Ip, Fx2, Mx2, My1, Mz1, My2, Mz2)
    Kg_dif = Kg - check_Kg
    if np.max(abs(Kg_dif)) > 10**-10: raise Exception("Incorrect Local Geometric Stiffness Matrix")

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

def true_Ke_ex2():
    Ke = np.array(
    [[3.043617640661192070e+01, -5.420865476416545370e+00, -5.420865476416545192e+01, 0.000000000000000000e+00, 1.665924512264597013e+01, -1.665924512264596924e+00, -3.043617640661192070e+01, 5.420865476416545370e+00, 5.420865476416545192e+01, 0.000000000000000000e+00, 1.665924512264597013e+01, -1.665924512264596924e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [-5.420865476416545370e+00, 4.416022119812502922e+00, 1.084173095283308896e+01, -1.665924512264597013e+01, 0.000000000000000000e+00, -8.329622561322985064e+00, 5.420865476416545370e+00, -4.416022119812502922e+00, -1.084173095283308896e+01, -1.665924512264597013e+01, 0.000000000000000000e+00, -8.329622561322985064e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [-5.420865476416545192e+01, 1.084173095283308896e+01, 1.117491585528600950e+02, 1.665924512264596924e+00, 8.329622561322985064e+00, 0.000000000000000000e+00, 5.420865476416545192e+01, -1.084173095283308896e+01, -1.117491585528600950e+02, 1.665924512264596924e+00, 8.329622561322985064e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, -1.665924512264597013e+01, 1.665924512264596924e+00, 1.175117521343565841e+02, 4.485181379173909022e+00, 4.485181379173916127e+01, 0.000000000000000000e+00, 1.665924512264597013e+01, -1.665924512264596924e+00, 5.074662360436773412e+01, 3.844441182149067160e+00, 3.844441182149070357e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [1.665924512264597013e+01, 0.000000000000000000e+00, 8.329622561322985064e+00, 4.485181379173912575e+00, 1.390406227543913644e+02, -8.970362758347830479e+00, -1.665924512264597013e+01, 0.000000000000000000e+00, -8.329622561322985064e+00, 3.844441182149068936e+00, 6.919994127868325506e+01, -7.688882364298139649e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [-1.665924512264596924e+00, -8.329622561322985064e+00, 0.000000000000000000e+00, 4.485181379173916127e+01, -8.970362758347830479e+00, 5.023403144674784926e+01, 1.665924512264596924e+00, 8.329622561322985064e+00, 0.000000000000000000e+00, 3.844441182149070357e+01, -7.688882364298139649e+00, -6.919994127868324796e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [-3.043617640661192070e+01, 5.420865476416545370e+00, 5.420865476416545192e+01, 0.000000000000000000e+00, -1.665924512264597013e+01, 1.665924512264596924e+00, 1.371148512302059430e+02, 8.330777240695574903e+01, 1.233782364836377354e+01, 0.000000000000000000e+00, 1.026581028768658399e+01, -3.423414936817881227e+01, -1.066786748235940081e+02, -8.872863788337230062e+01, -6.654647841252922547e+01, 0.000000000000000000e+00, 2.692505541033255412e+01, -3.590007388044340786e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [5.420865476416545370e+00, -4.416022119812502922e+00, -1.084173095283308896e+01, 1.665924512264597013e+01, 0.000000000000000000e+00, 8.329622561322985064e+00, 8.330777240695574903e+01, 1.110946969434065181e+02, 7.738820936536231443e+01, -1.026581028768658399e+01, 0.000000000000000000e+00, 4.422969644176639292e+01, -8.872863788337230062e+01, -1.066786748235940081e+02, -6.654647841252922547e+01, -2.692505541033255412e+01, 0.000000000000000000e+00, 3.590007388044340786e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [5.420865476416545192e+01, -1.084173095283308896e+01, -1.117491585528600950e+02, -1.665924512264596924e+00, -8.329622561322985064e+00, 0.000000000000000000e+00, 1.233782364836375933e+01, 7.738820936536230022e+01, 1.796090543024787110e+02, 3.423414936817881227e+01, -4.422969644176639292e+01, 0.000000000000000000e+00, -6.654647841252921125e+01, -6.654647841252921125e+01, -6.785989574961861592e+01, 3.590007388044340786e+01, -3.590007388044340786e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, -1.665924512264597013e+01, 1.665924512264596924e+00, 5.074662360436773412e+01, 3.844441182149067160e+00, 3.844441182149070357e+01, 0.000000000000000000e+00, -1.026581028768658399e+01, 3.423414936817881227e+01, 2.855056876005340314e+02, -7.283805467101187503e+01, -1.314061324590019808e+01, 0.000000000000000000e+00, 2.692505541033255412e+01, -3.590007388044340786e+01, 5.638152628659381094e+01, -6.627705947158780475e+01, -4.970779460369087133e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [1.665924512264597013e+01, 0.000000000000000000e+00, 8.329622561322985064e+00, 3.844441182149068936e+00, 6.919994127868325506e+01, -7.688882364298139649e+00, 1.026581028768658399e+01, 0.000000000000000000e+00, -4.422969644176639292e+01, -7.283805467101187503e+01, 3.070345582205687833e+02, -6.696278979598719161e+01, -2.692505541033255412e+01, 0.000000000000000000e+00, 3.590007388044340786e+01, -6.627705947158780475e+01, 5.638152628659381094e+01, -4.970779460369087133e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [-1.665924512264596924e+00, -8.329622561322985064e+00, 0.000000000000000000e+00, 3.844441182149070357e+01, -7.688882364298139649e+00, -6.919994127868324796e+00, -3.423414936817881227e+01, 4.422969644176639292e+01, 0.000000000000000000e+00, -1.314061324590019808e+01, -6.696278979598719161e+01, 2.520568826848816286e+02, 3.590007388044340786e+01, -3.590007388044340786e+01, 0.000000000000000000e+00, -4.970779460369087133e+01, -4.970779460369087133e+01, 8.537773980541349772e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.066786748235940081e+02, -8.872863788337230062e+01, -6.654647841252922547e+01, 0.000000000000000000e+00, -2.692505541033255412e+01, 3.590007388044340786e+01, 3.879311019587921123e+02, 4.686867041944811518e+00, 1.212675598307455402e+02, 0.000000000000000000e+00, -1.531527161361761102e+02, -8.389498661479770192e+01, -2.267249205292773127e+02, 1.133624602646386705e+02, -1.133624602646386563e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, -1.133624602646386421e+02, -5.452750660592079157e+01, -2.932068942321118143e+01, 5.864137884642235576e+01, 0.000000000000000000e+00, -1.286520046120490868e+01, -6.432600230602454339e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -8.872863788337230062e+01, -1.066786748235940081e+02, -6.654647841252922547e+01, 2.692505541033255412e+01, 0.000000000000000000e+00, -3.590007388044340786e+01, 4.686867041944811518e+00, 3.533745751385789617e+02, -8.032534119292219543e+01, 1.531527161361761102e+02, -8.881784197001252323e-16, -1.380054837415277404e+02, 1.133624602646386705e+02, -2.267249205292773127e+02, 1.133624602646386563e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, -2.932068942321118143e+01, -1.997097978570761612e+01, 3.350935934081277168e+01, 1.286520046120491045e+01, -8.881784197001252323e-16, 1.125705040355429709e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -6.654647841252921125e+01, -6.654647841252921125e+01, -6.785989574961861592e+01, -3.590007388044340786e+01, 3.590007388044340786e+01, 0.000000000000000000e+00, 1.212675598307455118e+02, -8.032534119292222385e+01, 3.648198350758227093e+02, 8.389498661479770192e+01, 1.380054837415277404e+02, 0.000000000000000000e+00, -1.133624602646386563e+02, 1.133624602646386563e+02, -2.267249205292773127e+02, 1.133624602646386421e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, 5.864137884642235576e+01, 3.350935934081277168e+01, -7.023501879692678074e+01, 6.432600230602454339e+00, -1.125705040355429531e+01, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -2.692505541033255412e+01, 3.590007388044340786e+01, 5.638152628659381094e+01, -6.627705947158780475e+01, -4.970779460369087133e+01, 0.000000000000000000e+00, 1.531527161361761102e+02, 8.389498661479770192e+01, 5.952316112150615481e+02, 2.051345875023111986e+01, -1.315831671227086304e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, -1.133624602646386421e+02, 1.220826495157646718e+02, 1.046422710135125556e+02, -1.046422710135126124e+02, 0.000000000000000000e+00, -1.286520046120490868e+01, -6.432600230602454339e+00, 3.278152040595481509e+01, -2.078224689886946663e+01, 4.156449379773893327e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 2.692505541033255412e+01, 0.000000000000000000e+00, -3.590007388044340786e+01, -6.627705947158780475e+01, 5.638152628659381094e+01, -4.970779460369087133e+01, -1.531527161361761102e+02, 8.881784197001252323e-16, 1.380054837415277404e+02, 2.051345875023112342e+01, 6.238072007010071047e+02, 9.179988500995129641e+01, 1.133624602646386421e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, 1.046422710135125556e+02, 1.220826495157646718e+02, 1.046422710135126124e+02, 1.286520046120491045e+01, -8.881784197001252323e-16, 1.125705040355429709e+01, -2.078224689886946663e+01, 5.727488282247954032e+01, 2.375113931299367209e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -3.590007388044340786e+01, 3.590007388044340786e+01, 0.000000000000000000e+00, -4.970779460369087133e+01, -4.970779460369087133e+01, 8.537773980541349772e+01, -8.389498661479770192e+01, -1.380054837415277404e+02, 0.000000000000000000e+00, -1.315831671227086588e+02, 9.179988500995129641e+01, 6.160716226752244893e+02, 1.133624602646386421e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, -1.046422710135126124e+02, 1.046422710135126124e+02, 1.220826495157647003e+02, 6.432600230602454339e+00, -1.125705040355429531e+01, 0.000000000000000000e+00, 4.156449379773892616e+01, 2.375113931299367209e+01, 2.164817385298902863e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -2.267249205292773127e+02, 1.133624602646386705e+02, -1.133624602646386563e+02, 0.000000000000000000e+00, 1.133624602646386421e+02, 1.133624602646386421e+02, 2.267249205292773127e+02, -1.133624602646386705e+02, 1.133624602646386563e+02, 0.000000000000000000e+00, 1.133624602646386421e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 1.133624602646386705e+02, -2.267249205292773127e+02, 1.133624602646386563e+02, -1.133624602646386421e+02, 0.000000000000000000e+00, 1.133624602646386421e+02, -1.133624602646386705e+02, 2.267249205292773127e+02, -1.133624602646386563e+02, -1.133624602646386421e+02, 0.000000000000000000e+00, 1.133624602646386421e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.133624602646386563e+02, 1.133624602646386563e+02, -2.267249205292773127e+02, -1.133624602646386421e+02, -1.133624602646386421e+02, 0.000000000000000000e+00, 1.133624602646386563e+02, -1.133624602646386563e+02, 2.267249205292773127e+02, -1.133624602646386421e+02, -1.133624602646386421e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 1.133624602646386421e+02, 1.133624602646386421e+02, 1.220826495157646718e+02, 1.046422710135125556e+02, -1.046422710135126124e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, -1.133624602646386421e+02, 3.313671915427898398e+02, 1.220826495157646150e+02, -1.220826495157647003e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.133624602646386421e+02, 0.000000000000000000e+00, 1.133624602646386421e+02, 1.046422710135125556e+02, 1.220826495157646718e+02, 1.046422710135126124e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, -1.133624602646386421e+02, 1.220826495157646150e+02, 3.313671915427898398e+02, 1.220826495157647003e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.133624602646386421e+02, -1.133624602646386421e+02, 0.000000000000000000e+00, -1.046422710135126124e+02, 1.046422710135126124e+02, 1.220826495157647003e+02, 1.133624602646386421e+02, 1.133624602646386421e+02, 0.000000000000000000e+00, -1.220826495157647003e+02, 1.220826495157647003e+02, 3.313671915427898966e+02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -5.452750660592079157e+01, -2.932068942321118143e+01, 5.864137884642235576e+01, 0.000000000000000000e+00, 1.286520046120490868e+01, 6.432600230602454339e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 5.452750660592079157e+01, 2.932068942321118143e+01, -5.864137884642235576e+01, 0.000000000000000000e+00, 1.286520046120490868e+01, 6.432600230602454339e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -2.932068942321118143e+01, -1.997097978570761612e+01, 3.350935934081277168e+01, -1.286520046120491045e+01, 8.881784197001252323e-16, -1.125705040355429709e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 2.932068942321118143e+01, 1.997097978570761612e+01, -3.350935934081277168e+01, -1.286520046120491045e+01, 8.881784197001252323e-16, -1.125705040355429709e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 5.864137884642235576e+01, 3.350935934081277168e+01, -7.023501879692678074e+01, -6.432600230602454339e+00, 1.125705040355429531e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -5.864137884642235576e+01, -3.350935934081277168e+01, 7.023501879692678074e+01, -6.432600230602454339e+00, 1.125705040355429531e+01, 0.000000000000000000e+00],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 1.286520046120490868e+01, 6.432600230602454339e+00, 3.278152040595481509e+01, -2.078224689886946663e+01, 4.156449379773893327e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.286520046120490868e+01, -6.432600230602454339e+00, 9.587048420609426103e+01, -2.424595471534771107e+01, 4.849190943069542925e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -1.286520046120491045e+01, 8.881784197001252323e-16, -1.125705040355429709e+01, -2.078224689886946663e+01, 5.727488282247954032e+01, 2.375113931299367209e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 1.286520046120491045e+01, -8.881784197001252323e-16, 1.125705040355429709e+01, -2.424595471534770752e+01, 1.244460736920397608e+02, 2.770966253182595551e+01],
    [0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, -6.432600230602454339e+00, 1.125705040355429531e+01, 0.000000000000000000e+00, 4.156449379773892616e+01, 2.375113931299367209e+01, 2.164817385298902863e+01, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00, 6.432600230602454339e+00, -1.125705040355429531e+01, 0.000000000000000000e+00, 4.849190943069542215e+01, 2.770966253182595551e+01, 8.288157989430084172e+01]]
    )
    return Ke