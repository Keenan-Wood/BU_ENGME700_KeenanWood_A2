import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from texttable import Texttable
import copy as cp

def load_frame(nodes: np.array, elements_in: list, xsections: list, constraint_list: list = [[]], applied_forces: list = [[]], N_pts: int = 10):
    elements = cp.deepcopy(elements_in)
    (N_nodes, N_elements) = (len(nodes), len(elements))
    if np.size(nodes, 1) < 6: nodes = np.hstack((nodes, np.zeros((N_nodes, 6 - np.size(nodes, 1)))))

    # Build constraint and force arrays
    (constrained_dof, forces) = (np.zeros((N_nodes, 6)), np.zeros((N_nodes, 6)))
    for constraint in constraint_list: constrained_dof[constraint[0]] = constraint[1:7]
    for force in applied_forces: forces[force[0], :] = force[1:7]
    
    # Append element lengths; Compute a z-vector (perpindicular to x_vec and [1,1,0] or [1,-1,0]) for each element where not provided
    for i, element in enumerate(elements):
        elements[i].append(np.linalg.norm(nodes[element[1], 0:3] - nodes[element[0], 0:3]))
        x_vec = (nodes[element[1], 0:3] - nodes[element[0], 0:3]) / elements[i][4]
        z_vec = element[3]
        if len(z_vec) == 0:
            if x_vec[0] < 0.95: z_vec = np.cross(x_vec, np.array([1, 0, 0]))
            else: z_vec = np.cross(x_vec, np.array([0, 1, 0]))
        z_vec = z_vec / np.linalg.norm(z_vec)
        # Insert corresponding coord transform matrix (gamma) into each element
        elements[i][3] = np.vstack((x_vec, np.cross(z_vec, x_vec), z_vec))

    # Assemble and partition elastic frame stiffness matrix
    Ke = assemble_stiffness_matrix(nodes, elements, xsections, local_stiffness)
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

    # Calculate element displacements and forces
    el_disps = np.array([np.hstack((
        el[3] @ disps[el[0], 0:3], el[3] @ disps[el[0], 3:6],
        el[3] @ disps[el[1], 0:3], el[3] @ disps[el[1], 3:6]
        )) for el in elements])
    el_forces = np.array([local_stiffness(xsections[el[2]], el[4], []) @ el_disps[i, :] for i, el in enumerate(elements)])

    # Calculate local interpolated displacements with bending shape functions (and linear axial displacement)
    el_lengths = np.array([np.linalg.norm(nodes[el[1], 0:3] + disps[el[0], 0:3] - nodes[el[0], 0:3] - disps[el[1], 0:3]) for el in elements])
    shapefuncs = bending_shape_functions(el_lengths, N_pts)
    el_x = np.linspace(el_disps[:,0], el_lengths + el_disps[:,6], N_pts)
    el_y = np.sum(shapefuncs * el_disps[np.newaxis, :, [1,5,7,11]], axis=2)
    el_z = np.sum(shapefuncs * el_disps[np.newaxis, :, [2,4,8,10]], axis=2)
    local_coords = np.stack((el_x, el_y, el_z), 2)
    # Transform to global coordinates
    el_gam = np.array([el[3].T for el in elements])
    start_coords = np.array([nodes[el[0], 0:3] for el in elements])
    inter_coords = start_coords[np.newaxis, ...] + np.squeeze(el_gam[np.newaxis, ...] @ local_coords[..., np.newaxis])

    # Assemble and Partition geometric stiffness matrix; Calculate critical load factor and vector
    Kg = assemble_stiffness_matrix(nodes, elements, xsections, geo_local_stiffness, el_forces)
    Kg_ff = Kg[ff_ind]
    (eig_vals, eig_vects) = sp.linalg.eig(Ke_ff, -Kg_ff)
    crit_ind = np.where(np.real(eig_vals) > 0, np.real(eig_vals), np.inf).argmin()
    crit_load_factor = np.real(eig_vals[crit_ind])
    crit_load_vec = np.zeros(6*N_nodes)
    crit_load_vec[free_ind] = eig_vects[crit_ind, :]
    crit_load_vec = crit_load_vec.reshape((N_nodes, 6))

    return (disps, forces, el_disps, el_forces, inter_coords, crit_load_factor, crit_load_vec)

def bending_shape_functions(L, N_pts):
    # Calculate Hermite shape functions for bending displacement interpolation
    x = np.outer(np.linspace(0,1,N_pts), L)
    L = np.outer(np.ones(N_pts), L)
    shape_mat = np.array([1 - (x/L)**2*(3 - 2*x/L), x*(1 - x/L)**2, (x/L)**2*(3 - 2*x/L), x**2/L*(x/L - 1)])
    shape_mat = np.permute_dims(shape_mat, (1, 2, 0))
    return shape_mat

def assemble_stiffness_matrix(nodes, elements, xsections, fun_local, forces = None):
    # Assemble element stiffness matrices to get frame stiffness matrix
    if forces is None: forces = np.zeros((len(elements), 12))
    K_full = np.zeros((6*len(nodes), 6*len(nodes)))
    for i, el in enumerate(elements):
        # Calculate stiffness matrix in local coordinates; Transform to global coordinates
        local_K = fun_local(xsections[el[2]], el[4], forces[i,:])
        gam_full = sp.sparse.block_diag((el[3], el[3], el[3], el[3])).toarray()
        global_K = np.matmul(gam_full.T, np.matmul(local_K, gam_full))
        # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
        a_rng = np.arange(6*el[0], 6*(el[0]+1))
        b_rng = np.arange(6*el[1], 6*(el[1]+1))
        K_full[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += global_K
    return K_full

def local_stiffness(xsec: list, L: float, _):
    # Calculate elastic stiffness matrix (in local coordinates)
    (E, A, Iy, Iz, Ip, J, v) = tuple(xsec)
    k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz])
    k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz])
    k_cross_1 = 6/L * np.array([Iz, -Iy, 0, -Iy, Iz])
    k_cross_2 = 6/L * np.array([Iz, -Iy, 0, Iy, -Iz])
    # Build the quadrants of the full matrix, then combine
    k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
    k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
    k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_2, -1))
    k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_1, -1))
    K_local = E/L * np.vstack((np.hstack((k_a, k_b)), np.hstack((k_c, k_d))))
    return K_local

def geo_local_stiffness(xsec: list, L: float, forces: np.array):
    # Calculate geometric stiffness matrix (in local coordinates)
    (_, A, _, _, Ip, _, _) = tuple(xsec)
    (_, _, _, _, My1, Mz1) = tuple(forces[0:6])
    (Fx2, _, _, Mx2, My2, Mz2) = tuple(forces[6:12])
    k_diag_1 = Fx2/L * np.array([1, 6/5, 6/5, Ip/A, 2/15*L**2, 2/15*L**2])
    k_diag_2 = -k_diag_1 + Fx2 * np.array([0, 0, 0, 0, L/10, L/10])
    # Build the quadrants of the full matrix (neglecting lower triangle), then combine
    (k_a, k_d) = (np.zeros((6, 6)), np.zeros((6,6)))
    k_a[np.ix_([1,2,3], [3,4,5])] = 1/L * np.array([[My1, Mx2, Fx2*L/10], [Mz1, -Fx2*L/10, Mx2], [0, -L*(2*Mz1-Mz2)/6, L*(2*My1-My2)/6]])
    k_d[np.ix_([1,2,3], [3,4,5])] = 1/L * np.array([[-My2, Mx2, -Fx2*L/10], [-Mz2, Fx2*L/10, Mx2], [0, -L*(2*Mz2-Mz1)/6, L*(2*My2-My1)/6]])
    k_b = -k_d - np.transpose(k_a) + np.diag(k_diag_2)
    k_b[3, 4:6] = np.array([-(Mz1+Mz2)/6, (My1+My2)/6])
    k_b[4:6, 3] = k_b[3, 4:6]
    k_b[np.array([4,5]), np.array([5,4])] = np.array([Mx2/2, -Mx2/2])
    Kg = np.vstack((np.hstack((k_a, k_b)), np.hstack((np.zeros((6,6)), k_d))))
    # Add the transpose to fill the lower triangle (and make symmetric)
    Kg += np.transpose(Kg) + np.diag(np.hstack((k_diag_1, k_diag_1)))
    return Kg

def print_results(disps, forces, el_disps, el_forces, crit_load_factor, crit_load_vec):
    res_nodes = Texttable()
    res_nodes_geo = Texttable()
    res_elements = Texttable()
    print('Node Displacements and Forces')
    res_nodes.add_rows([['Node'] + [str(i) for i in range(0,len(disps))],
                ['x'] + list(disps[:,0]), ['y'] + list(disps[:,1]), ['z'] + list(disps[:,2]),
                [chr(952) + 'x'] + list(disps[:,3]), [chr(952) + 'y'] + list(disps[:,4]), [chr(952) + 'z'] + list(disps[:,5]),
                ['Fx'] + list(forces[:,0]), ['Fy'] + list(forces[:,1]), ['Fz'] + list(forces[:,2]),
                ['Mx'] + list(forces[:,3]), ['My'] + list(forces[:,4]), ['Mz'] + list(forces[:,5])])
    print(res_nodes.draw())

    print('\nCritical Load Factor: ' + str(crit_load_factor))
    print('Critical Load Vector')
    res_nodes_geo.add_rows([['Node'] + [str(i) for i in range(0,len(disps))],
                ['eig_x'] + list(crit_load_vec[:,0]), ['eig_y'] + list(crit_load_vec[:,1]), ['eig_z'] + list(crit_load_vec[:,2]),
                ['eig_' + chr(952) + 'x'] + list(crit_load_vec[:,3]), ['eig_' + chr(952) + 'y'] + list(crit_load_vec[:,4]), ['eig_' + chr(952) + 'z'] + list(crit_load_vec[:,5])])
    print(res_nodes_geo.draw())

    print('\nLocal Coordinate Displacements and Internal Forces')
    res_elements.add_rows([['Element'] + [str(i) for i in range(0,len(el_disps))],
                ['dx'] + list(el_disps[:,0]), ['dy'] + list(el_disps[:,1]), ['dz'] + list(el_disps[:,2]),
                [chr(952) + 'x'] + list(el_disps[:,3]), [chr(952) + 'y'] + list(el_disps[:,4]), [chr(952) + 'z'] + list(el_disps[:,5]),
                ['Fx'] + list(el_forces[:,0]), ['Fy'] + list(el_forces[:,1]), ['Fz'] + list(el_forces[:,2]),
                ['Mx'] + list(el_forces[:,3]), ['My'] + list(el_forces[:,4]), ['Mz'] + list(el_forces[:,5])])
    print(res_elements.draw())


def plot_frame(inter_coords: np.array, elements):
    # Plot interpolated coordinates
    ax = plt.figure().add_subplot(projection='3d')
    clrs = ['b', 'g', 'r', 'm', 'k']
    for i_el in range(0, len(elements)): 
        ax.plot(inter_coords[:, i_el, 0], inter_coords[:, i_el, 1], inter_coords[:, i_el, 2], color=clrs[elements[i_el][2]])
    plt.show()
    return