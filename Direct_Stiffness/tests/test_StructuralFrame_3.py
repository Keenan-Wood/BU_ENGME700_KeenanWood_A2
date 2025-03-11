import numpy as np
#import sys, os
#sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'src'))
from StructuralFrame_3 import *

def test_load_frame_simple():
    # Frame geometry definition
    nodes = np.array([[0,0,10], [15,0,10], [15,0,0]])
    elements = [[0, 1, 0, [0,0,1]], [1, 2, 0, [1,0,0]]]

    # Cross section list
    (E, v) = (1000, 0.3)
    (b, h, J) = (.5, 1, 0.02861)
    xsection = [[E, v, 'rectangle', [b, h, J]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1], [2,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, -0.05, 0.075, 0.1, -0.05, 0.1, -0.25]]

    # Create frame, apply loads, and display results
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    result = simple_frame.deformation

    (all_disps, all_forces, crit_factor, crit_vec) = (result["disps"], result["forces"], result["crit_load_factor"], result["crit_load_vec"])

    # Check displacements
    DISPS_MATCH = True
    true_disps = np.zeros((len(nodes), 6))
    true_disps[1, :] = np.array([-0.00134143, 1.24465464, 0.00190536, -0.12134412, 0.0167805, 0.10196546])
    true_disps[2, 3:6] = np.array([-0.12602614, -0.00859146, 0.10196546])
    if np.linalg.norm(all_disps - true_disps) >= 10**-6: DISPS_MATCH = False

    # Check forces
    FORCES_MATCH = True
    true_forces = np.zeros((len(nodes), 6))
    true_forces[0, :] = np.array([0.04471418, -0.07109832, -0.00473182, 0.0890168, 0.02383551, -0.8164748])
    true_forces[1, :] = np.array([-0.05, 0.075, 0.1, -0.05, 0.1, -0.25])
    true_forces[2, 0:3] = np.array([0.00528582, -0.00390168, -0.09526818])
    if np.linalg.norm(all_forces - true_forces) >= 10**-6: FORCES_MATCH = False

    #np.set_printoptions(precision=10)
    #np.set_printoptions(suppress=True)
    #print('Calculated Displacement:')
    #print(all_disps)
    #print('\nTrue Displacements:')
    #print(true_disps)
    #print('\nCalculated Forces:')
    #print(all_forces)
    #print('\nTrue Forces:')
    #print(true_forces)

    assert DISPS_MATCH and FORCES_MATCH

def test_load_frame_simple_2():
    # Frame geometry definition
    nodes = np.array([[0,0,0], [-5,1,10], [-1,5,13], [-3,7,11], [6,9,5]])
    elements = [[0,1], [1,2], [2,3], [2,4]]

    # Cross section list
    (E, v) = (500, 0.3)
    r = 1
    xsection = [[E, v, 'circle', [r]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,0,0,1,0,0,0], [3,1,1,1,1,1,1], [4,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, 0.05, 0.05, -0.1, 0, 0, 0], [2, 0, 0, 0, -0.1, -0.1, 0.3]]

    # Create frame, apply loads, and display results
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    result = simple_frame.deformation

    (all_disps, all_forces, crit_factor, crit_vec) = (result["disps"], result["forces"], result["crit_load_factor"], result["crit_load_vec"])

    # Check displacements
    DISPS_MATCH = True
    true_disps = np.zeros((len(nodes), 6))
    true_disps[0, :] = np.array([.12500826, .06604393, 0, .00498268, -.00977281, .00237119])
    true_disps[1, :] = np.array([.02913522, .00351598, -.04180902, .00523327, -.00851986, .00237119])
    true_disps[2, :] = np.array([8.66002551e-5, 7.23622353e-4, 4.38264195e-4, 1.8077038e-3, -2.97736985e-3, 2.52180509e-3])
    true_disps[4, 3:6] = np.array([-.00257306, .00062871, .00049162])
    if np.linalg.norm(all_disps - true_disps) >= 10**-6: DISPS_MATCH = False

    # Check forces
    FORCES_MATCH = True
    true_forces = np.zeros((len(nodes), 6))
    for force in forces:
        true_forces[force[0], :] = force[1:7]
    true_forces[0, 2] = .01753352
    true_forces[3, :] = np.array([-.03892948, -.02361856, .10543765, -.22303444, .1294369, -.28470812])
    true_forces[4, 0:3] = np.array([-.01107052, -.02638144, -.02297117])
    if np.linalg.norm(all_forces - true_forces) >= 10**-6: FORCES_MATCH = False

    #np.set_printoptions(precision=10)
    #np.set_printoptions(suppress=True)
    #print('Calculated Displacement:')
    #print(all_disps)
    #print('\nTrue Displacements:')
    #print(true_disps)
    #print('\nCalculated Forces:')
    #print(all_forces)
    #print('\nTrue Forces:')
    #print(true_forces)

    assert DISPS_MATCH and FORCES_MATCH

def test_critical_load():
    # Frame geometry definition
    nodes = np.array([[0,0,0], [30,40,0]])
    elements = [[0,1]]

    # Cross section list
    (E, v) = (1000, 0.3)
    r = 1
    xsection = [[E, v, 'circle', [r]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, -3/5, -4/5, 0, 0, 0, 0]]

    # Create frame, apply loads, and display results
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    result = simple_frame.deformation

    crit_factor = result["crit_load_factor"]
    true_crit_factor = .7809879011060754
    CRIT_FACTORS_MATCH = abs(crit_factor - true_crit_factor) < 10**-6
    
    assert CRIT_FACTORS_MATCH

def test_z_vec_creation_printing_plotting():
    NOEXCEPTION = True
    # Frame geometry definition
    (L1, L2, L3, L4) = (11, 23, 15, 13)
    x = [0, L1, L1, 0, 0, L1, L1, 0, 0, L1, L1, 0]
    y = [0, 0, L2, L2, 0, 0, L2, L2, 0, 0, L2, L2]
    z = [0,0,0,0, L3,L3,L3,L3, L3+L4,L3+L4,L3+L4,L3+L4]
    nodes = np.array([np.array([x[i], y[i], z[i], 0, 0, 0]) for i in range(0, 12)])
    elements = [[i, i+4, 0] for i in range(0, 8)]
    elements.extend([4*lvl + i, 4*lvl + (i+1)%4, 1] for i in range(0, 4) for lvl in [1,2])

    # Cross section list
    (E1, v1, E2, v2) = (10000, 0.3, 50000, 0.3)
    (r, b, h) = (1, 0.5, 1)
    xsection = [[E1, v1, 'circle', [r]], [E2, v2, 'rectangle', [b, h]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[i,1,1,1,1,1,1] for i in range(0,4)]

    # Force list (node_id, forces on each DOF)
    forces = [[i,0,0,-1,0,0,0] for i in range(8,12)]

    # Create frame, apply loads, and display results
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    print("\nTechnical Correctness 1 - Problem 3:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed()
    divided_frame = simple_frame.subdivide(2)
    simple_frame.plot_deformed("buckled")

    assert NOEXCEPTION

# Debugging tests:
#test_load_frame_simple()
#test_load_frame_simple_2()
#test_critical_load()