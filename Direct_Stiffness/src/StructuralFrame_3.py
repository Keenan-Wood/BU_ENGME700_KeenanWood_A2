import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from texttable import Texttable
import copy as cp

class frame:

    def __init__(nodes: np.array, element_list: list, xsection_list: list, constraints: list = [[]]):
        if np.size(nodes, 1) < 6: nodes = np.hstack((nodes, np.zeros(len(nodes), 6 - np.size(nodes, 1))))
        self.nodes = nodes

        # Build constraint array
        self.fixed_dof = np.zeros((N_nodes, 6))
        for constraint in constraints: self.fixed_dof[constraint[0]] = constraint[1:7]

        # Create list of elements
        self.elements = []
        for el in element_list: self.elements.append(element(el[0], el[1], el[2], el[3]))

        # Assemble and Partition elastic stiffness matrix
        self.Ke = np.zeros((6*len(nodes), 6*len(nodes)))
        for el in self.elements:
            # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
            a_rng = np.arange(6*el.node_a, 6*(el.node_a + 1))
            b_rng = np.arange(6*el.node_b, 6*(el.node_b + 1))
            self.Ke[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += el.Ke
        self.free_ind = np.flatnonzero(self.fixed_dof == 0)
        self.fixed_ind = np.flatnonzero(self.fixed_dof)
        ff_ind = np.ix_(free_ind, free_ind)
        sf_ind = np.ix_(fixed_ind, free_ind)
        (self.Ke_ff, self.Ke_sf) = (self.Ke[ff_ind], self.Ke[sf_ind])

    def bending_shape_functions(L, N_pts):
        # Calculate Hermite shape functions for bending displacement interpolation
        x = np.outer(np.linspace(0,1,N_pts), L)
        L = np.outer(np.ones(N_pts), L)
        shape_mat = np.array([1 - (x/L)**2*(3 - 2*x/L), x*(1 - x/L)**2, (x/L)**2*(3 - 2*x/L), x**2/L*(x/L - 1)])
        shape_mat = np.permute_dims(shape_mat, (1, 2, 0))
        return shape_mat

    def load_frame(self, applied_forces: list, N_pts: int = 10):
        (N_nodes, N_elements) = (len(self.nodes), len(self.elements))

        #Build force array
        forces = np.zeros((N_nodes, 6))
        for force in applied_forces: forces[force[0], :] = force[1:7]

        # Calculate free displacements and support forces; Combine with known displacements and forces
        disps = np.zeros((N_nodes, 6))
        (disps_vec, forces_vec) = (disps.reshape(-1), forces.reshape(-1))
        free_disps = np.matmul(np.linalg.inv(self.Ke_ff), forces_vec[self.free_ind])
        support_forces = np.matmul(self.Ke_sf, free_disps)
        disps_vec[self.free_ind] = free_disps
        forces_vec[self.fixed_ind] = support_forces

        # Calculate element displacements and forces
        el_disps = np.array([np.hstack((
            el.gamma @ disps[el.node_a, 0:3], el.gamma @ disps[el.node_a, 3:6],
            el.gamma @ disps[el.node_b, 0:3], el.gamma @ disps[el.node_b, 3:6]
            )) for el in self.elements])
        el_forces = np.array([el.Ke @ el_disps[i, :] for i, el in enumerate(self.elements)])

        # Calculate local interpolated displacements with bending shape functions (and linear axial displacement)
        el_lengths = np.array([np.linalg.norm(nodes[el.node_b, 0:3] + disps[el.node_a, 0:3] - nodes[el.node_a, 0:3] - disps[el.node_b, 0:3]) for el in self.elements])
        shapefuncs = bending_shape_functions(el_lengths, N_pts)
        el_x = np.linspace(el_disps[:,0], el_lengths + el_disps[:,6], N_pts)
        el_y = np.sum(shapefuncs * el_disps[np.newaxis, :, [1,5,7,11]], axis=2)
        el_z = np.sum(shapefuncs * el_disps[np.newaxis, :, [2,4,8,10]], axis=2)
        local_coords = np.stack((el_x, el_y, el_z), 2)
        # Transform to global coordinates
        el_gam = np.array([el.gamma.T for el in self.elements])
        start_coords = np.array([nodes[el.node_a, 0:3] for el in self.elements])
        inter_coords = start_coords[np.newaxis, ...] + np.squeeze(el_gam[np.newaxis, ...] @ local_coords[..., np.newaxis])

        # Assemble and Partition geometric stiffness matrix
        Kg = np.zeros((6*N_nodes, 6*N_nodes))
        for el in self.elements:
            Kg_global = geo_local_stiffness(el, el_forces)
            # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
            a_rng = np.arange(6*el.node_a, 6*(el.node_a + 1))
            b_rng = np.arange(6*el.node_b, 6*(el.node_b + 1))
            Kg[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += Kg_global
        Kg_ff = Kg[ff_ind]

        # Calculate critical load factor and vector
        (eig_vals, eig_vects) = sp.linalg.eig(self.Ke_ff, -Kg_ff)
        crit_ind = np.where(np.real(eig_vals) > 0, np.real(eig_vals), np.inf).argmin()
        crit_load_factor = np.real(eig_vals[crit_ind])
        crit_load_vec = np.zeros(6*N_nodes)
        crit_load_vec[self.free_ind] = eig_vects[crit_ind, :]
        crit_load_vec = crit_load_vec.reshape((N_nodes, 6))

        return (disps, forces, el_disps, el_forces, inter_coords, crit_load_factor, crit_load_vec)

    def print_deformed_results(disps, forces, el_disps, el_forces, crit_load_factor, crit_load_vec):
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

    def plot_deformed(inter_coords: np.array, elements):
        # Plot interpolated coordinates
        ax = plt.figure().add_subplot(projection='3d')
        clrs = ['b', 'g', 'r', 'm', 'k']
        for i_el in range(0, len(elements)): 
            ax.plot(inter_coords[:, i_el, 0], inter_coords[:, i_el, 1], inter_coords[:, i_el, 2], color=clrs[elements[i_el][2]])
        plt.show()
        return

class element:

    def __init__(nodes, node_a: int, node_b: int, xsec: xsection, z_vec = []):
        self.a = node_a
        self.b = node_b
        self.xsec = xsec
        self.L = np.linalg.norm(nodes[node_b, 0:3] - nodes[node_a, 0:3])

        # Calculate coord transform and elastic stiffness matrices
        x_vec = (nodes[node_b, 0:3] - nodes[node_a, 0:3]) / self.L
        if len(z_vec) == 0:
            if x_vec[0] < 0.95: z_vec = np.cross(x_vec, np.array([1, 0, 0]))
            else: z_vec = np.cross(x_vec, np.array([0, 1, 0]))
        z_vec = z_vec / np.linalg.norm(z_vec)
        self.gamma = np.vstack((x_vec, np.cross(z_vec, x_vec), z_vec))
        self.calc_elastic_stiffness(self)

    def calc_elastic_stiffness(self):
        # Calculate elastic stiffness matrix (in local coordinates)
        xsec = self.xsec
        (E, A, Iy, Iz, Ip, J, v) = (xsec.mat.E, xsec.A, xsec.Iy, xsec.Iz, xsec.Ip, xsec.J, xsec.v)
        k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz])
        k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz])
        k_cross_1 = 6/L * np.array([Iz, -Iy, 0, -Iy, Iz])
        k_cross_2 = 6/L * np.array([Iz, -Iy, 0, Iy, -Iz])
        # Build the quadrants of the full matrix, then combine
        k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
        k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
        k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_2, -1))
        k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_1, -1))
        Ke_local = E/L * np.vstack((np.hstack((k_a, k_b)), np.hstack((k_c, k_d))))
        # Transform to global coordinates
        gamma_full = sp.sparse.block_diag((self.gamma, self.gamma, self.gamma, self.gamma)).toarray()
        self.Ke = np.matmul(gamma_full.T, np.matmul(Ke_local, gamma_full))

    def geo_local_stiffness(self, forces: np.array):
        # Calculate geometric stiffness matrix (in local coordinates)
        (L, A, Ip) = (self.L, self.xsec.A, self.xsec.Ip)
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
        # Transform to global coordinates
        gamma_full = sp.sparse.block_diag((self.gamma, self.gamma, self.gamma, self.gamma)).toarray()
        Kg_global = np.matmul(gamma_full.T, np.matmul(Kg, gamma_full))
        return Kg_global

class xsection:

    def __init__(material_values: list, geometric_values: list = []):
        (self.E, self.v) = tuple(material_values)
        if len(geometric_values) > 0:
            (self.A, self.Iy, self.Iz, self.Ip, self.J) = tuple(geometric_values)
    
    def make_circular(self, r):
        (self.A, self.Iy, self.Iz, self.Ip, self.J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)

    def make_rectangular(self, b, h, J):
        (self.A, self.Iy, self.Iz, self.Ip, self.J) = (b*h, h*b**3/12, b*h**3/12, b*h*(b**2+h**2)/12, J)