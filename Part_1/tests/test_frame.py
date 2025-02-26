import numpy as np
import matplotlib.pyplot as plt
from StructuralFrame import *

def test_example_1():
    # Build node initialization arguments
    #    node inputs-> coords: list = [0,0,0,0,0,0], fixed_dof: list = []
    (x_0, y_0, z_0) = (-1, 1, 0)
    (x_1, y_1, z_1) = (.8, .7, 0)
    (x_2, y_2, z_2) = (0, 0, 0)
    coords_1 = [x_0, y_0, z_0, 0, 0, 0]
    coords_2 = [x_1, y_1, z_1, 0, 0, 0]
    coords_3 = [x_2, y_2, z_2, 0, 0, 0]
    node_coords = [coords_1, coords_2, coords_3]
    fixed_nodes = [0]
    pinned_nodes = [2]
    fixed_dof = []
    for fixed_node in fixed_nodes:
        fixed_dof.append(6*fixed_node + range(0,6))
    for pinned_node in pinned_nodes:
        fixed_dof.append(6*pinned_node + range(0,3))
    nodes = [node_coords, fixed_dof]

    # Build element initialization arguments
    #    element inputs-> mat: material, x_sec: xsection, node_a: node, node_b: node
    # Create test material
    mat1 = material('test_material', 10**6)
    # Create rectangular cross section
    b = .01
    h = .02
    sec_A = b*h
    sec_I_y = b*h**3/12
    sec_I_z = h*b**3/12
    sec_J = b*h*(b**2 + h**2)/12
    xsec1 = xsection(sec_A, sec_I_y, sec_I_z, sec_J)
    # Assign node pairs and fill element args
    node_pairs = [[0,1], [1,2]]
    elements = []
    for node_pair in node_pairs:
        elements.append([mat1, xsec1, node_pair[0], node_pair[1]])

    # Create frame instance with node and element data
    test_frame = frame(nodes, elements)

    # Format applied loads
    F = 1
    M = 1
    forces = [[0, F, 0, 0, 0, M]]
    node_ids = [1]
    (all_disps, all_forces) = test_frame.calc_apply_load(forces, node_ids)
    assert True