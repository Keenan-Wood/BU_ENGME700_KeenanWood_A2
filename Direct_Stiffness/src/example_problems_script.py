# %%
#%matplotlib ipympl
import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), 'src'))
from StructuralFrame_3 import *

# Functions defining each problem, running solver, and printing/plotting results
# %%
# Example 1 - Part 1 Code Review 1
def solve_CR1_P1_ex1():
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
    N_pts = 30
    simple_frame = frame(nodes, xsection, elements, constraints)
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = simple_frame.load_frame(forces, N_pts)
    print("\nCode Review 1 - Problem 1:\n")
    simple_frame.print_deformed_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    simple_frame.plot_deformed(inter_coords)

# Example 2 - Part 1 Code Review 1
def solve_CR1_P1_ex2():
    # Frame geometry definition
    nodes = np.array([[0,0,0], [-5,1,10], [-1,5,13], [-3,7,11], [6,9,5]])
    elements = [[0,1,0,[]], [1,2,0,[]], [2,3,0,[]], [2,4,0,[]]]

    # Cross section list
    (E, v) = (500, 0.3)
    r = 1
    xsection = [[E, v, 'circle', [r]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,0,0,1,0,0,0], [3,1,1,1,1,1,1], [4,1,1,1,0,0,0]]

    # Force list (node_id, forces on each DOF)
    forces = [[1, 0.05, 0.05, -0.1, -0.05, 0.1, -0.25]]

    # Create frame, apply loads, and display results
    N_pts = 30
    simple_frame = frame(nodes, xsection, elements, constraints)
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = simple_frame.load_frame(forces, N_pts)
    print("\nCode Review 1 - Problem 2:\n")
    simple_frame.print_deformed_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    simple_frame.plot_deformed(inter_coords)

# Problem 1 and 2 - Technical Correctness 1
def solve_T1_problem_1_2():
    # Frame geometry definition
    (x, y, z) = (np.linspace(0, 25, 7), np.linspace(0, 50, 7), np.linspace(0, 37, 7))
    nodes = np.array([np.array([x[i], y[i], z[i], 0, 0, 0]) for i in range(0, 7)])
    elements = [[i, i+1, 0, []] for i in range(0, 6)]

    # Cross section list
    (E, v) = (10000, 0.3)
    r = 1
    xsection = [[E, v, 'circle', [r]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[0,1,1,1,1,1,1]]

    # Problem 1 - Force list (node_id, forces on each DOF)
    forces = [[6, 0.05, -0.1, 0.23, 0.1, -0.025, -0.08]]

    # Create frame, apply loads, and display results
    N_pts = 30
    simple_frame = frame(nodes, xsection, elements, constraints)
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = simple_frame.load_frame(forces, N_pts)
    print("\nTechnical Correctness 1 - Problem 1:\n")
    simple_frame.print_deformed_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    simple_frame.plot_deformed(inter_coords)

    # Problem 2 - Force list (node_id, forces on each DOF)
    P = 1
    L = np.linalg.norm(np.array([25, 50, 37]))
    (Fx, Fy, Fz) = (-25*P/L, -50*P/L, -37*P/L)
    forces_2 = [[6, Fx, Fy, Fz, 0, 0, 0]]

    # Create frame, apply loads, and display results
    (all_disps_2, all_forces_2, el_disps_2, el_forces_2, inter_coords_2, crit_factor_2, crit_vec_2) = simple_frame.load_frame(forces_2, N_pts)
    print("\nTechnical Correctness 1 - Problem 2:\n")
    simple_frame.print_deformed_results(all_disps_2, all_forces_2, el_disps_2, el_forces_2, crit_factor_2, crit_vec_2)
    simple_frame.plot_deformed(inter_coords_2)

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
    (E1, v1, E2, v2) = (10000, 0.3, 50000, 0.3)
    (r, b, h, J2) = (1, 0.5, 1, 0.028610026041666667)
    xsection = [[E1, v1, 'circle', [r]], [E2, v2, 'rectangle', [b, h, J2]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[i,1,1,1,1,1,1] for i in range(0,4)]

    # Force list (node_id, forces on each DOF)
    forces = [[i,0,0,-1,0,0,0] for i in range(8,12)]

    # Create frame, apply loads, and display results
    N_pts = 30
    simple_frame = frame(nodes, xsection, elements, constraints)
    (all_disps, all_forces, el_disps, el_forces, inter_coords, crit_factor, crit_vec) = simple_frame.load_frame(forces, N_pts)
    print("\nTechnical Correctness 1 - Problem 3:\n")
    simple_frame.print_deformed_results(all_disps, all_forces, el_disps, el_forces, crit_factor, crit_vec)
    simple_frame.plot_deformed(inter_coords)

# %%
solve_CR1_P1_ex1()
solve_CR1_P1_ex2()
solve_T1_problem_1_2()
solve_T1_problem_3()