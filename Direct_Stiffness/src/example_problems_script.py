# %%
#%matplotlib ipympl
import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'src'))
import StructuralFrame_2 as sf
from StructuralFrame_3 import *

# %%

# Example 1 - Part 1 Code Review 1
def solve_CR1_P1_ex1():
    # Frame geometry definition
    nodes = np.array([[0,0,10], [15,0,10], [15,0,0]])
    elements = [[0, 1, 0, [0,0,1]], [1, 2, 0, [1,0,0]]]

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

    N_pts = 30
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = sf.load_frame(nodes, elements, xsection, constraints, forces, N_pts)
    sf.print_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    #sf.plot_frame(inter_coords, elements)

    # Test StructuralFrame_3
    N_pts = 30
    xsection = [[E, v]]
    elements = [[0, 1, 0, [0,0,1]], [1, 2, 0]]
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.xsecs[0].make_rectangular(b, h, J)
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = simple_frame.load_frame(forces, N_pts)
    print_deformed_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    simple_frame.plot_deformed(inter_coords)

solve_CR1_P1_ex1()

# %%

# Example 2 - Part 1 Code Review 1
def solve_CR1_P1_ex2():
    # Frame geometry definition
    nodes = np.array([[0,0,0], [-5,1,10], [-1,5,13], [-3,7,11], [6,9,5]])
    elements = [[0,1,0,[]], [1,2,0,[]], [2,3,0,[]], [2,4,0,[]]]

    # Cross section list
    E = 500
    r = 1
    (A, I_y, I_z, I_p, J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
    v = 0.3
    xsection = [[E, A, I_y, I_z, I_p, J, v]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,0,0,1,0,0,0], [3,1,1,1,1,1,1], [4,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, 0.05, 0.05, -0.1, -0.05, 0.1, -0.25]]
#############
    N_pts = 30
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = sf.load_frame(nodes, elements, xsection, constraints, forces, N_pts)
    sf.print_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    #sf.plot_frame(inter_coords, elements)

solve_CR1_P1_ex2()

# %%
# Problem 1 and 2 - Technical Correctness 1
def solve_T1_problem_1_2():
    # Frame geometry definition
    (x, y, z) = (np.linspace(0, 25, 7), np.linspace(0, 50, 7), np.linspace(0, 37, 7))
    nodes = np.array([np.array([x[i], y[i], z[i], 0, 0, 0]) for i in range(0, 7)])
    elements = [[i, i+1, 0, []] for i in range(0, 6)]

    # Cross section list
    E = 10000
    r = 1
    (A, I_y, I_z, I_p, J) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
    v = .3
    xsection = [[E, A, I_y, I_z, I_p, J, v]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1]]

    # Problem 1 - Force list (node_id, forces on each DOF)
    forces = [[6, 0.05, -0.1, 0.23, 0.1, -0.025, -0.08]]

    N_pts = 30
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = sf.load_frame(nodes, elements, xsection, constraints, forces, N_pts)
    sf.print_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    sf.plot_frame(inter_coords, elements)

    # Display Problem 1 results
    print("\nProblem 1:\n")
    print("Reaction force at Node 0:")
    print(all_forces[0, 0:3])
    print("Reaction moment at Node 0:")
    print(all_forces[0, 3:6])
    print("Displacement at Node 3:")
    print(all_disps[3, 0:3])
    print("Rotation at Node 3:")
    print(all_disps[3, 3:6])
    print("Displacement at Node 6:")
    print(all_disps[6, 0:3])
    print("Rotation at Node 6:")
    print(all_disps[6, 3:6])

    # Problem 2 - Force list (node_id, forces on each DOF)
    P = 1
    L = np.linalg.norm(np.array([25, 50, 37]))
    Fx = -25*P/L
    Fy = -50*P/L
    Fz = -37*P/L
    forces_2 = [[6, Fx, Fy, Fz, 0, 0, 0]]

    (all_disps_2, all_forces_2, el_disps_2, el_forces_2, inter_coords_2, crit_factor_2, crit_vec_2) = sf.load_frame(nodes, elements, xsection, constraints, forces_2, N_pts)
    sf.plot_frame(inter_coords, elements)

    # Display Problem 2 Result
    print("\nProblem 2:\n")
    print("Critical Load Factor:")
    print(crit_factor_2)

solve_T1_problem_1_2()

# %%
# Problem 3 - Technical Correctness 1
def solve_T1_problem_3():
    # Frame geometry definition
    (L1, L2, L3, L4) = (11, 23, 15, 13)
    x = [0, L1, L1, 0, 0, L1, L1, 0, 0, L1, L1, 0]
    y = [0, 0, L2, L2, 0, 0, L2, L2, 0, 0, L2, L2]
    z = [0,0,0,0, L3,L3,L3,L3, L3+L4,L3+L4,L3+L4,L3+L4]
    nodes = np.array([np.array([x[i], y[i], z[i], 0, 0, 0]) for i in range(0, 12)])
    z_vec2 = np.array([0, 0, 1])
    elements = [[i, i+4, 0, []] for i in range(0, 8)]
    elements.extend([4*lvl + i, 4*lvl + (i+1)%4, 1, z_vec2] for i in range(0, 4) for lvl in [1,2])

    # Cross section list
    (E1, E2) = (10000, 50000)
    (r, b, h) = (1, 0.5, 1)
    (v1, v2) = (0.3, 0.3)
    (A1, I_y1, I_z1, I_p1, J1) = (np.pi*r**2, np.pi*r**4/4, np.pi*r**4/4, np.pi*r**4/2, np.pi*r**4/2)
    (A2, I_y2, I_z2, I_p2, J2) = (b*h, h*b**3/12, b*h**3/12, b*h*(b**2 + h**2)/12, 0.028610026041666667)
    xsection = [[E1, A1, I_y1, I_z1, I_p1, J1, v1], [E2, A2, I_y2, I_z2, I_p2, J2, v2]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[i,1,1,1,1,1,1] for i in range(0,4)]

    # Force list (node_id, forces on each DOF)
    forces = [[i,0,0,-1,0,0,0] for i in range(8,12)]

    N_pts = 30
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = sf.load_frame(nodes, elements, xsection, constraints, forces, N_pts)
    sf.plot_frame(inter_coords, elements)

    # Display Problem 3 Result
    print("\nProblem 3:\n")
    print("Critical Load Factor:")
    print(crit_factor)

solve_T1_problem_3()
# %%
