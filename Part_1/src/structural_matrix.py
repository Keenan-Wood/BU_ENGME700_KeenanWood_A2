import numpy as np
from dataclasses import dataclass

# A class to list all nodes
@dataclass
class node:
    id: int

    def __init__(self, coords: list = [0,0,0,0,0,0], supported: list = [], exisiting_nodes: list = [], id = None):
        N_nodes = len(existing_nodes)
        if N_nodes > 0:
            existing_node_ids = [existing_nodes[i_node].id for i_node in range(0, N_nodes)]
            if id in existing_node_ids:
                raise Exception('Non-unique node id provided')
            elif id is None:
                id = min(existing_node_ids, N_nodes)
        elif id is None: id = 0
        if not isinstance(id, int): raise Exception('Invalid node id type')
        elif id < 0: raise Exception('Node id must be positive')
        self.id = id
        self.add_support(self, supported)
            
    def add_support(self, support: int):
        if self.supported is None: self.supported = []
        self.supported = list(set(self.supported).union(support))

@dataclass
class material:
    E: float

    def __init__(self, name: str, E: float):
        self.name = name
        self.E = E

@dataclass
class xsection:
    A: float
    I_y: float
    I_z: float
    J: float

    def __init__(self, A: float, I_y: float, I_z: float, J: float, sub_xsection = None, sub_axis = None):
        self.A = A
        self.I_y = I_y
        self.I_z = I_z
        self.J = J
        self.sub_xsection = sub_xsection
        self.sub_axis = sub_axis

@dataclass
class element:
    L_0: float

    def __init__(self, mat: material, x_sec: xsection, node_a: node, node_b: node, L_0: float = None, x_sec_angle: float = 0):
        self.mat: mat
        self.x_sec = x_sec
        self.node_a = node_a
        self.node_b = node_b
        if L_0 is None: self.L_0 = np.linalg.norm(node_a.coords - node_b.coords)
        else: self.L_0 = L_0
        self.x_sec_angle = x_sec_angle

    def get_local_stiffness(self):
        E = self.mat.E
        L = self.L_0
        A = self.x_sec.A
        c_Iy = self.x_sec.I_y / A
        c_Iz = self.x_sec.I_z / A
        c_J = self.x_sec.J / A
        k_diag_1 = np.array([1, 12*c_Iz/L**2, 12*c_Iy/L**2, c_J/(2*(1+v)), 4*c_Iy, 4*c_Iz])
        k_diag_2 = -k_diag_1 + np.array([0, 0, 0, 0, 6*c_Iy, 6*c_Iz])
        k_cross_1 = 6/L**2 * np.array([c_Iz, -c_Iy, 0, -c_Iy, c_Iz])
        k_cross_2 = 6/L**2 * np.array([c_Iz, -c_Iy, 0, c_Iy, -c_Iz])
        k_a = np.diag(k_diag_1) + np.fliplr(np.diag(k_cross_1, -1))
        k_b = np.diag(k_diag_2) + np.fliplr(np.diag(k_cross_2, -1))
        k_c = np.diag(k_diag_2) - np.fliplr(np.diag(k_cross_1, -1))
        k_d = np.diag(k_diag_1) - np.fliplr(np.diag(k_cross_2, -1))
        return E*A/L * np.vstack(np.hstack(k_a, k_b), np.hstack(k_c, k_d))

    def calc_coord_transform(self):
        x_vec = self.node_b.coords[0:2] - self.node_a.coords[0:2]
        z_vec = self.z_vec
        y_vec = np.cross(z_vec, x_vec)
        gam_small = np.hstack(x_vec, y_vec, z_vec)
        gam_full = np.zeros(12, 12)
        for i in range(0,4):
            gam_full[3*i:3*i+2, 3*i:3*i+2] = gam_small
        self.Gam = gam_full

    def calc_global_stiffness(self):
        if self.Gam is None:
            self.calc_coord_transform()
        K_local = self.calc_local_stiffness()
        self.K_e = np.matmul(np.transpose(self.Gam), np.matmul(K_local, self.Gam))

@dataclass
class frame:

    def __init__(self, nodes: list = [], elements: list = []):
        if all(isinstance(input_node, node) for input_node in nodes):
            self.nodes = nodes
        else:
            self.add_nodes(self, nodes)
        if all(isinstance(input_element, element) for input_element in elements):
            self.elements = elements
        else:
            self.add_elements(self, elements)

    def add_nodes(self, node_arg_list: list):
        for i_node in range(0, len(node_arg_list)):
            self.nodes.append(node(*[node_arg_list[i_node]]))

    def add_elements(self, element_arg_list: list):
        for i_element in range(0, len(element_arg_list)):
            self.elements.append(element(*[element_arg_list[i_element]]))

    def deform(self, forces):
        # Update node positions and internal forces using the direct stiffness method
        self.nodes.coords

