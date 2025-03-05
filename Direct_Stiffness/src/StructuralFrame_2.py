import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def load_frame(nodes: np.array, elements: list, xsections_list: list, constraint_list: list = [[]], applied_forces: list = [[]]):
    (N_nodes, N_elements) = (len(nodes), len(elements))

    # Build constraint and force arrays
    (constrained_dof, forces) = (np.zeros((N_nodes, 6)), np.zeros((N_nodes, 6)))
    for constraint in constraint_list: constrained_dof[constraint[0]] = constraint[1:7]
    for force in applied_forces: forces[force[0], :] = force[1:7]
    
    # Append element lengths; Compute a z-vector for each element where not provided
    for i, element in enumerate(elements):
        elements[i].append(np.linalg.norm(nodes[element[1], 0:3] - nodes[element[0], 0:3]))
        if len(element[3]) == 0:
            x_vec = nodes[element[1], 0:3] - nodes[element[0], 0:3]
            elements[i][3] = get_assumed_z_vec(x_vec)

    # Assemble and partition elastic frame stiffness matrix
    Ke = assemble_stiffness_matrix(nodes, elements, xsections_list, np.zeros((N_elements, 12)), calc_local_stiffness)
    free_ind = np.flatnonzero(constrained_dof == 0)
    fixed_ind = np.flatnonzero(constrained_dof)
    ff_ind = np.ix_(free_ind, free_ind)
    sf_ind = np.ix_(fixed_ind, free_ind)
    (Ke_ff, Ke_sf) = (Ke[ff_ind], Ke[sf_ind])

    # Calculate free displacements and support forces; Combine with known displacements and forces
    disps = np.zeros((N_nodes, 6))
    (disps_vec, forces_vec) = (disps.reshape(-1), forces.reshape(-1))
    free_disps = np.matmul(np.linalg.inv(Ke_ff), forces_vec[free_ind])
    support_forces = np.matmul(Ke_sf, free_disps)
    disps_vec[free_ind] = free_disps
    forces_vec[fixed_ind] = support_forces

    # Calculate element forces
    ax = plt.figure().add_subplot(projection='3d')
    (el_disps, el_forces) = (np.zeros((N_elements, 12)), np.zeros((N_elements, 12)))
    for i, element in enumerate(elements):
        x_vec = nodes[element[1], 0:3] - nodes[element[0], 0:3]
        gam_full = calc_transformation_matrix(x_vec, element[3])
        el_disp = np.hstack((disps[element[0],:], disps[element[1],:]))
        el_disps[i,:] = np.matmul(gam_full, el_disp)
        el_Ke = calc_local_stiffness(xsections_list[element[2]], element[4], [])
        el_forces[i,:] = np.matmul(el_Ke, el_disps[i,:])

        # Plot interpolated element displacement using shape functions
        f_N = lambda x, L: np.array([1 - (x/L)**2*(3 - 2*x/L), x*(1 - x/L)**2, (x/L)**2*(3 - 2*x/L), x**2/L*(x/L - 1)])
        f_v = lambda x, L, vec: np.matmul(f_N(x, L), vec)
        el_x = np.arange(0, element[4], 100)
        el_y = np.array([f_v(x, element[4], el_disps[i,[1,5,6,11]]) for x in el_x])
        el_z = np.array([f_v(x, element[4], el_disps[i,[2,4,7,10]]) for x in el_x])
        y_vec = np.cross(x_vec, element[3])
        plt_x = el_x*x_vec[0] + el_y*y_vec[0] + el_z*element[3][0]
        plt_y = el_x*x_vec[1] + el_y*y_vec[1] + el_z*element[3][1]
        plt_z = el_x*x_vec[2] + el_y*y_vec[2] + el_z*element[3][2]
        ax.plot(el_x, el_y, el_z)

    plt.show()

    # Assemble and Partition geometric stiffness matrix; Calculate critical load factor and vector
    Kg = assemble_stiffness_matrix(nodes, elements, xsections_list, el_forces, calc_geometric_local_stiffness)
    Kg_ff = Kg[ff_ind]
    (eig_vals, eig_vects) = sp.linalg.eig(Ke_ff, -Kg_ff)
    crit_ind = np.where(np.real(eig_vals) > 0, np.real(eig_vals), np.inf).argmin()
    crit_load_factor = np.real(eig_vals[crit_ind])
    crit_load_vec = np.zeros(6*N_nodes)
    crit_load_vec[free_ind] = eig_vects[crit_ind, :]
    crit_load_vec = crit_load_vec.reshape((N_nodes, 6))

    return (disps, forces, crit_load_factor, crit_load_vec)

def get_assumed_z_vec(x_vec: np.array):
    # Return a vector orthogonal to x_vec using the cross product; Cross x_vec with [1,-1,0] if x_vec is nearly parallel to [1,1,0]
    z_vec = np.cross(x_vec, np.array([1, 1, 0]))
    if np.linalg.norm(z_vec) / np.linalg.norm(x_vec) < 0.01:
        z_vec = np.cross(x_vec, np.array([1, -1, 0]))
    return z_vec

def assemble_stiffness_matrix(nodes, elements, xsections_list, forces, fun_local):
    K_full = np.zeros((6*len(nodes), 6*len(nodes)))
    for i, element in enumerate(elements):
        # Calculate stiffness matrix in local coordinates
        local_K = fun_local(xsections_list[element[2]], element[4], forces[i,:])
        
        # Transform stiffness matrix to global coordinates with full (12x12) rotation matrix (gam_full)
        x_vec = nodes[element[1], 0:3] - nodes[element[0], 0:3]
        gam_full = calc_transformation_matrix(x_vec, element[3])
        global_K = np.matmul(gam_full.T, np.matmul(local_K, gam_full))

        # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
        a_rng = np.arange(6*element[0], 6*(element[0]+1))
        b_rng = np.arange(6*element[1], 6*(element[1]+1))
        K_full[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += global_K
    return K_full

def calc_transformation_matrix(x_vec: np.array, z_vec: np.array):
    gam = np.vstack((x_vec, np.cross(z_vec, x_vec), z_vec))
    row_norms = np.linalg.norm(gam, axis=1)
    gam = gam / row_norms[:, np.newaxis]
    gam_full = sp.sparse.block_diag((gam, gam, gam, gam)).toarray()
    return gam_full

def calc_local_stiffness(xsec: list, L: float, _):
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

def calc_geometric_local_stiffness(xsec: list, L: float, forces: np.array):
    (_, A, _, _, Ip, _, _) = tuple(xsec)
    (_, _, _, _, My1, Mz1) = tuple(forces[0:6])
    (Fx2, _, _, Mx2, My2, Mz2) = tuple(forces[6:12])
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