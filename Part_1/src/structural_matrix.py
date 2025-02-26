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

    def __init__(self, name: str = 'anonymous', E: float = 0):
        self.name = name
        self.E = E

@dataclass
class xsection:

    def __init__(self, A: float = 0, I_y: float = 0, I_z: float = 0, J: float = 0, sub_xsection = None, sub_axis = None):
        self.A = A
        self.I_y = I_y
        self.I_z = I_z
        self.J = J
        self.sub_xsection = sub_xsection
        self.sub_axis = sub_axis

@dataclass
class element:

    def __init__(self, mat: material, x_sec: xsection, node_a: node, node_b: node, L_0 = None, z_vec = None):
        (L_0, z_vec) = handle_inputs_element(mat, x_sec, node_a, node_b, L_0, z_vec)
        self.mat = mat
        self.x_sec = x_sec
        self.node_a = node_a
        self.node_b = node_b
        self.L_0 = L_0
        self.z_vec = z_vec
        # Calculate and assign self.Gam and self.K_e
        self.calc_coord_transform()
        self.calc_global_stiffness()

    def handle_inputs_element(mat, x_sec, node_a, node_b, L_0, z_vec):
        # Check for valid material properties
        if mat.E <= 0: raise Exception('Material modulus of elasticity missing')

        # Check for valid cross section properties
        if x_sec.A <= 0: raise Exception('Invalid cross section area')
        if x_sec.I_y <= 0: raise Exception('Invalid cross section I_y')
        if x_sec.I_z <= 0: raise Exception('Invalid cross section I_z')
        if x_sec.J <= 0: raise Exception('Invalid cross section J')

        # Check rest length - set if not provided
        if L_0 is None: self.L_0 = np.linalg.norm(node_a.coords - node_b.coords)
        else:
            if not isinstance(L_0, (int, float)): raise Exception('Rest length L_0 must be a float')

        # Check cross section orientation vector (z_vec) - set if not provided; normalize
        if z_vec is None:
            z_vec = get_assumed_z_vec(self)
        else:
            if isinstance(z_vec, list): z_vec = np.array(z_vec)
            if not isinstance(z_vec, np.array): raise Exception('Cross section z_vector must be a numpy array')
            if len(z_vec) != 3: raise Exception('Cross section z_vector must be an array with 3 elements (x,y,z)')
        z_vec = z_vec / np.linalg.norm(z_vec)

        return (L_0, z_vec)

    def get_assumed_z_vec(self):
        x_vec = self.node_b.coords[0:2] - self.node_a.coords[0:2]
        z_vec = np.cross(x_vec, np.array([1, 1, 0]))
        if np.linalg.norm(z_vec) / np.linalg.norm(x_vec) < 0.01:
            z_vec = np.cross(x_vec, np.array([1, -1, 0]))
        return z_vec

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
        K_local = self.get_local_stiffness()
        self.K_e = np.matmul(np.transpose(self.Gam), np.matmul(K_local, self.Gam))

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

    def handle_inputs_compose_arrays(ar1, ar2, a, b):
        w2 = len(ar2)/2
        if w2 > floor(w2): raise Exception('Small array must have even dimensions')
        if len(ar2) != len(ar2[0]): raise Exception('Small array must be square')
        if len(ar1) != len(ar1[0]): raise Exception('Large array must be square')
        if len(ar1)/w2 > floor(len(ar1)/w2): raise Exception('Size of large array must be multiple of size of small array')
        if a < 0 or b < 0: raise Exception('Indices cannot be negative')
        if a == b: raise Exception('Indices must be different')
        if a >= len(ar1)/w2 or b >= len(ar1)/w2: raise Exception('Indices provided are too large')
        # Swap 'a' and 'b' if a > b
        if a > b:
            return (b, a)
        else:
            return (a, b)

    def compose_arrays(ar1: np.array, ar2: np.array, a: int, b: int):
        (a, b) = handle_inputs_compose_arrays(ar1, ar2, a, b)
        w2 = len(ar2)/2
        a_rng = range(a*w2, (a+1)*w2)
        b_rng = range(b*w2, (b+1)*w2)
        ar1[a_rng, a_rng] = ar1[a_rng, a_rng] + ar2[0:w2, 0:w2]
        ar1[a_rng, b_rng] = ar1[a_rng, b_rng] + ar2[0:w2, w2:-1]
        ar1[b_rng, a_rng] = ar1[b_rng, a_rng] + ar2[w2:-1, 0:w2]
        ar1[b_rng, b_rng] = ar1[b_rng, b_rng] + ar2[w2:-1, w2:-1]
        return ar1

    def assemble_stiffness_matrix(self):
        self.K_e = np.zeros(len(self.nodes), len(self.nodes))
        for elem in self.elements:
            self.K_e = compose_arrays(self.K_e, elem.K_e, elem.node_a.id, elem.node_b.id)




    def deform(self, forces):
        # Update node positions and internal forces using the direct stiffness method
        self.nodes.coords

