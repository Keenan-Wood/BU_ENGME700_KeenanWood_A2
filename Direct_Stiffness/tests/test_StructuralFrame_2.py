import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'src'))
from StructuralFrame_2 import load_frame

def test_load_frame_simple():
    # Frame geometry definition
    nodes = np.array([[0,0,10,0,0,0], [15,0,10,0,0,0], [15,0,0,0,0,0]])
    zvec = np.array([[0,0,1], [1,0,0]])
    elements = [[0,1,0,zvec[0,:]], [1,2,0,zvec[1,:]]]

    # Cross section list
    E = 1000
    (b, h) = (.5, 1)
    (A, I_y, I_z, I_p, J) = (b*h, h*b**3/12, b*h**3/12, b*h*(b**2+h**2)/12, .02861)
    v = .3
    xsection = [[E, A, I_y, I_z, I_p, J, v]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1], [2,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, -0.05, 0.075, 0.1, -0.05, 0.1, -0.25]]

    (all_disps, all_forces, crit_factor, crit_vec) = load_frame(nodes, elements, xsection, constraints, forces)

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

    np.set_printoptions(precision=10)
    np.set_printoptions(suppress=True)
    print('Calculated Displacement:')
    print(all_disps)
    print('\nTrue Displacements:')
    print(true_disps)
    print('\nCalculated Forces:')
    print(all_forces)
    print('\nTrue Forces:')
    print(true_forces)

    assert DISPS_MATCH and FORCES_MATCH

def test_load_frame_simple_2():
    # Frame geometry definition
    nodes = np.array([[0,0,0,0,0,0], [-5,1,10,0,0,0], [-1,5,13,0,0,0],[-3,7,11,0,0,0],[6,9,5,0,0,0]])
    elements = [[0,1,0,[]], [1,2,0,[]], [2,3,0,[]], [2,4,0,[]]]

    # Cross section list
    E = 500
    r = 1
    (A, I_y, I_z, I_p, J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
    v = .3
    xsection = [[E, A, I_y, I_z, I_p, J, v]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,0,0,1,0,0,0], [3,1,1,1,1,1,1], [4,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, 0.05, 0.05, -0.1, 0, 0, 0], [2, 0, 0, 0, -0.1, -0.1, 0.3]]

    (all_disps, all_forces, crit_factor, crit_vec) = load_frame(nodes, elements, xsection, constraints, forces)

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

    np.set_printoptions(precision=10)
    np.set_printoptions(suppress=True)
    print('Calculated Displacement:')
    print(all_disps)
    print('\nTrue Displacements:')
    print(true_disps)
    print('\nCalculated Forces:')
    print(all_forces)
    print('\nTrue Forces:')
    print(true_forces)

    assert DISPS_MATCH and FORCES_MATCH

def test_critical_load():
    # Frame geometry definition
    nodes = np.array([[0,0,0,0,0,0], [30,40,0,0,0,0]])
    elements = [[0,1,0,[]]]

    # Cross section list
    E = 1000
    r = 1
    (A, I_y, I_z, I_p, J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
    v = .3
    xsection = [[E, A, I_y, I_z, I_p, J, v]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, -3/5, -4/5, 0, 0, 0, 0]]

    (all_disps, all_forces, crit_factor, crit_vec) = load_frame(nodes, elements, xsection, constraints, forces)
    true_crit_factor = .7809879011060754
    CRIT_FACTORS_MATCH = abs(crit_factor - true_crit_factor) < 10**-6
    
    assert CRIT_FACTORS_MATCH