import numpy as np
from dataclasses import dataclass

# A class to list all nodes
@dataclass
class node:
    x: float
    y: float
    z: float

    def __init__(self, name: str, E: float, H: float, yield_stress: float):
        self.name = name
        self.E = E
        self.H = H
        self.yield_stress = yield_stress

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

    def __init__(self, mat: material, x_sec: xsection, node_a: node, node_b: node, L_0: float):
        self.mat: mat
        self.x_sec = x_sec
        self.node_a = node_a
        self.node_b = node_b
        self.L_0 = L_0

@dataclass
class constraint:

    def __init__(self, restrained_node: node, restricted_DOF: int):
       self.restrained_node = restrained_node
       self.restricted_DOF = restricted_DOF

@dataclass
class force:

    def __init__(self, forced_node: node, applied_force):
       self.forced_node = forced_node
       self.applied_force = applied_force

@dataclass
class frame:

    def __init__(self, nodes, elements, constraints, forces):
        self.nodes = nodes
        self.elements = elements
        self.constraints = constraints
        self.forces = forces

    def deform():
        # Update node positions and internal forces using the direct stiffness method
        self.nodes.coords = (0,0,0)
