import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from texttable import Texttable

class frame:
    def __init__(self, nodes: np.array, xsection_list: list, element_list: list, constraints: list = [[]]):
        # Build node array, xsection list, element list, and constraint array
        if np.size(nodes, 1) < 6: nodes = np.hstack((nodes, np.zeros((len(nodes), 6 - np.size(nodes, 1)))))
        self.nodes = nodes
        self.xsecs = [xsection(i, xsec[0], xsec[1], xsec[2], xsec[3]) for i, xsec in enumerate(xsection_list)]
        self.elements = [element(nodes, el[0], el[1], self.xsecs[el[2]], el[3] if 3 < len(el) else []) for el in element_list]
        self.fixed_dof = np.zeros((len(self.nodes), 6))
        for constraint in constraints: self.fixed_dof[constraint[0]] = constraint[1:7]

        # Assemble and Partition elastic stiffness matrix
        self.Ke = np.zeros((6*len(nodes), 6*len(nodes)))
        for el in self.elements:
            # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
            a_rng = np.arange(6*el.a, 6*(el.a + 1))
            b_rng = np.arange(6*el.b, 6*(el.b + 1))
            self.Ke[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += el.Ke
        self.free_ind = np.flatnonzero(self.fixed_dof == 0)
        self.fixed_ind = np.flatnonzero(self.fixed_dof)
        self.Ke_ff = self.Ke[np.ix_(self.free_ind, self.free_ind)]
        self.Ke_sf = self.Ke[np.ix_(self.fixed_ind, self.free_ind)]

    def bending_shape_functions(self, L, N_pts):
        # Calculate Hermite shape functions for bending displacement interpolation
        x = np.outer(np.linspace(0,1,N_pts), L)
        L = np.outer(np.ones(N_pts), L)
        shape_mat = np.array([1 - (x/L)**2*(3 - 2*x/L), x*(1 - x/L)**2, (x/L)**2*(3 - 2*x/L), x**2/L*(x/L - 1)])
        return np.permute_dims(shape_mat, (1, 2, 0))

    def apply_load(self, applied_forces: list, N_pts: int = 10):
        #Build force array
        forces = np.zeros((len(self.nodes), 6))
        for force in applied_forces: forces[force[0], :] = force[1:7]

        # Calculate free displacements and support forces; Combine with known displacements and forces
        disps = np.zeros((len(self.nodes), 6))
        (disps_vec, forces_vec) = (disps.reshape(-1), forces.reshape(-1))
        free_disps = np.matmul(np.linalg.inv(self.Ke_ff), forces_vec[self.free_ind])
        support_forces = np.matmul(self.Ke_sf, free_disps)
        disps_vec[self.free_ind] = free_disps
        forces_vec[self.fixed_ind] = support_forces

        # Calculate element displacements and forces
        el_disps = np.array([np.hstack((
            el.gamma @ disps[el.a, 0:3], el.gamma @ disps[el.a, 3:6],
            el.gamma @ disps[el.b, 0:3], el.gamma @ disps[el.b, 3:6]
            )) for el in self.elements])
        el_forces = np.array([el.Ke_local @ el_disps[i, :] for i, el in enumerate(self.elements)])
        self.deformation = {"disps": disps, "forces": forces, "el_disps": el_disps, "el_forces": el_forces}

        # Calculate local interpolated displacements with bending shape functions (and linear axial displacement)
        el_lengths = np.array([np.linalg.norm(self.nodes[el.b, 0:3] + disps[el.a, 0:3] - self.nodes[el.a, 0:3] - disps[el.b, 0:3]) for el in self.elements])
        shapefuncs = self.bending_shape_functions(el_lengths, N_pts)
        el_x = np.linspace(el_disps[:,0], el_lengths + el_disps[:,6], N_pts)
        el_y = np.sum(shapefuncs * el_disps[np.newaxis, :, [1,5,7,11]], axis=2)
        el_z = np.sum(shapefuncs * el_disps[np.newaxis, :, [2,4,8,10]], axis=2)
        local_coords = np.stack((el_x, el_y, el_z), 2)
        # Transform to global coordinates
        el_gamma = np.array([el.gamma.T for el in self.elements])
        start_coords = np.array([self.nodes[el.a, 0:3] for el in self.elements])
        self.deformation["inter_coords"] = start_coords[np.newaxis, ...] + np.squeeze(el_gamma[np.newaxis, ...] @ local_coords[..., np.newaxis])

        # Assemble and Partition geometric stiffness matrix
        Kg = np.zeros((6*len(self.nodes), 6*len(self.nodes)))
        for i, el in enumerate(self.elements):
            Kg_global = el.geo_local_stiffness(el_forces[i,:])
            # Add quadrants of global stiffness matrix to full matrix in positions defined by nodes
            a_rng = np.arange(6*el.a, 6*(el.a + 1))
            b_rng = np.arange(6*el.b, 6*(el.b + 1))
            Kg[np.ix_(np.concatenate((a_rng,b_rng)), np.concatenate((a_rng,b_rng)))] += Kg_global
        Kg_ff = Kg[np.ix_(self.free_ind, self.free_ind)]

        # Calculate critical load factor and vector
        (eig_vals, eig_vects) = sp.linalg.eig(self.Ke_ff, -Kg_ff)
        crit_ind = np.where(np.real(eig_vals) > 0, np.real(eig_vals), np.inf).argmin()
        self.deformation["crit_load_factor"] = np.real(eig_vals[crit_ind])
        crit_load_vec = np.zeros(6*len(self.nodes))
        crit_load_vec[self.free_ind] = eig_vects[crit_ind, :]
        self.deformation["crit_load_vec"] = crit_load_vec.reshape((len(self.nodes), 6))

    def print_deformed_results(self):
        (disps, forces, el_disps, el_forces, crit_load_vec) = tuple([self.deformation[key] for key in ["disps", "forces", "el_disps", "el_forces", "crit_load_vec"]])
        (res_nodes, res_nodes_geo, res_elements) = (Texttable(max_width=0), Texttable(max_width=0), Texttable(max_width=0))
        print('Node Displacements and Forces')
        res_nodes.set_precision(8)
        res_nodes.set_cols_align(["c" for node in self.nodes])
        res_nodes.add_rows([['Node'] + [str(i) for i in range(0,len(disps))],
                    ['x'] + list(disps[:,0]), ['y'] + list(disps[:,1]), ['z'] + list(disps[:,2]),
                    [chr(952) + 'x'] + list(disps[:,3]), [chr(952) + 'y'] + list(disps[:,4]), [chr(952) + 'z'] + list(disps[:,5]),
                    ['Fx'] + list(forces[:,0]), ['Fy'] + list(forces[:,1]), ['Fz'] + list(forces[:,2]),
                    ['Mx'] + list(forces[:,3]), ['My'] + list(forces[:,4]), ['Mz'] + list(forces[:,5])])
        print(res_nodes.draw())
        
        print('\nCritical Load Factor: ' + str(self.deformation["crit_load_factor"]))
        print('Critical Load Vector')
        res_nodes_geo.set_precision(8)
        res_nodes_geo.set_cols_align(["c" for node in self.nodes])
        res_nodes_geo.add_rows([['Node'] + [str(i) for i in range(0,len(disps))],
                    ['eig_x'] + list(crit_load_vec[:,0]), ['eig_y'] + list(crit_load_vec[:,1]), ['eig_z'] + list(crit_load_vec[:,2]),
                    ['eig_' + chr(952) + 'x'] + list(crit_load_vec[:,3]), ['eig_' + chr(952) + 'y'] + list(crit_load_vec[:,4]), ['eig_' + chr(952) + 'z'] + list(crit_load_vec[:,5])])
        print(res_nodes_geo.draw())

        print('\nLocal Coordinate Displacements and Internal Forces')
        res_elements.set_precision(8)
        res_elements.set_cols_align(["c" for el in self.elements])
        res_elements.add_rows([['Element'] + [str(i) for i in range(0,len(el_disps))],
                    ['dx'] + list(el_disps[:,0]), ['dy'] + list(el_disps[:,1]), ['dz'] + list(el_disps[:,2]),
                    [chr(952) + 'x'] + list(el_disps[:,3]), [chr(952) + 'y'] + list(el_disps[:,4]), [chr(952) + 'z'] + list(el_disps[:,5]),
                    ['Fx'] + list(el_forces[:,0]), ['Fy'] + list(el_forces[:,1]), ['Fz'] + list(el_forces[:,2]),
                    ['Mx'] + list(el_forces[:,3]), ['My'] + list(el_forces[:,4]), ['Mz'] + list(el_forces[:,5])])
        print(res_elements.draw())

    def plot_deformed(self):
        inter_coords = self.deformation["inter_coords"]
        ax = plt.figure().add_subplot(projection='3d')
        clrs = ['b', 'g', 'r', 'm', 'k']
        for i, el in enumerate(self.elements): 
            ax.plot(inter_coords[:,i,0], inter_coords[:,i,1], inter_coords[:,i,2], color=clrs[el.xsec.id % len(clrs)])
        plt.show()

class element:
    def __init__(self, nodes, node_a: int, node_b: int, xsec, z_vec = []):
        self.a = node_a
        self.b = node_b
        self.L = np.linalg.norm(nodes[node_b, 0:3] - nodes[node_a, 0:3])
        if isinstance(xsec, list): self.xsec = xsection(xsec)
        else: self.xsec = xsec

        # Calculate coord transform and elastic stiffness matrices
        x_vec = (nodes[node_b, 0:3] - nodes[node_a, 0:3]) / self.L
        if len(z_vec) == 0:
            if x_vec[0] < 0.95: z_vec = np.cross(x_vec, np.array([1, 0, 0]))
            else: z_vec = np.cross(x_vec, np.array([0, 1, 0]))
        z_vec = z_vec / np.linalg.norm(z_vec)
        self.gamma = np.vstack((x_vec, np.cross(z_vec, x_vec), z_vec))
        self.calc_elastic_stiffness()

    def calc_elastic_stiffness(self):
        # Calculate elastic stiffness matrix (in local coordinates)
        xsec = self.xsec
        (L, E, v, A, Iy, Iz, Ip, J) = (self.L, xsec.E, xsec.v, xsec.A, xsec.Iy, xsec.Iz, xsec.Ip, xsec.J)
        k_diag_1 = np.array([A, 12*Iz/L**2, 12*Iy/L**2, J/(2*(1+v)), 4*Iy, 4*Iz])
        k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*Iy, 6*Iz])
        k_cross_1 = 6/L * np.array([Iz, -Iy, 0, -Iy, Iz])
        k_cross_2 = 6/L * np.array([Iz, -Iy, 0, Iy, -Iz])
        # Build the quadrants of the full matrix, then combine
        k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
        k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
        k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_2, -1))
        k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_1, -1))
        self.Ke_local = E/L * np.vstack((np.hstack((k_a, k_b)), np.hstack((k_c, k_d))))
        # Transform to global coordinates
        gamma_full = sp.sparse.block_diag((self.gamma, self.gamma, self.gamma, self.gamma)).toarray()
        self.Ke = np.matmul(gamma_full.T, np.matmul(self.Ke_local, gamma_full))

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
    def __init__(self, id: int, E: float, v: float, geometry_type: str, geometry_parameters: list):
        (self.id, self.E, self.v) = (id, E, v)
        match geometry_type:
            case 'circle':
                r = geometry_parameters[0]
                (self.A, self.Iy, self.Iz, self.Ip, self.J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
            case 'rectangle':
                if len(geometry_parameters) == 3:
                    (b, h, J) = tuple(geometry_parameters)
                else:
                    (b, h) = tuple(geometry_parameters)
                    (L2, L1) = tuple(sorted([b, h]))
                    J = L1*L2**3/16 * (16/3 - 3.36*L2/L1*(1 - (L2/L1)**4/12))
                (self.A, self.Iy, self.Iz, self.Ip, self.J) = (b*h, h*b**3/12, b*h**3/12, b*h*(b**2+h**2)/12, J)
            case _:
                (self.A, self.Iy, self.Iz, self.Ip, self.J) = tuple(geometry_parameters)