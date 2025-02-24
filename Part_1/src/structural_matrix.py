import numpy as np
from dataclasses import dataclass

# A class to list all nodes
@dataclass
class node:
    id: int

    def __init__(self, coords: list = [0,0,0], exisiting_nodes: list = [], id = None):
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
            
    def add_constraint(self, constraint_type: int):
        self.constraint_type = constraint_type

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
        self.nodes.coords = (0,0,0)
