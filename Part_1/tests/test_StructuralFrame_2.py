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
    (A, I_y, I_z, J) = (b*h, h*b**3/12, b*h**3/12, b*h*(b**2+h**2)/12)
    v = .3
    xsection = [[E, A, I_y, I_z, J, v]]

    # Constraint list (node_id, type) - 1 for pinned, 2 for fixed
    constraints = [[0,2], [2,1]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, -0.05, 0.075, 0.1, -0.05, 0.1, -0.25]]

    (all_disps, all_forces) = load_frame(nodes, elements, xsection, constraints, forces)

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

    return DISPS_MATCH and FORCES_MATCH

print(str(test_load_frame_simple()))

# Debug Notes:

# Disp Errors: 
# Node 1: y, tx, tz
# Node 2:    tx, tz

# Force Errors:
# Node 1: y, Mx, Mz
# Node 2:    Mx, Mz
