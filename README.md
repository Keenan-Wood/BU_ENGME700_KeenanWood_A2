# Matrix Structural Analysis

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1)
[![tests](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions/workflows/tests.yml/badge.svg)](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions)

---

### Table of Contents
* [The Method](#algo)
* [Conda environment, installation, and testing](#install)
* [Documentation & Examples](#tutorial)
* [More Information](#more)

---

### Conda environment, install, and testing <a name="install"></a>

To install this package, please begin by setting up a conda environment and activating it. For example:
```bash
conda create --name me700-env python=3.12
conda activate me700-env
```

Navigate to the project directory (./DirectStiffness) and create an editable install of the code:
```bash
pip install -e .
```

Test that the code is working with pytest:
```bash
pytest -v --cov=StructuralFrame_3 --cov-report term-missing
```

If you are using VSCode to run this code, don't forget to set VSCode virtual environment to the newly-activated environment.

---

#### **Documentation**

**Frame Class**  

    def __init__(self, nodes: np.array, xsection_list: list, element_list: list, constraints: list = [[]]):
        Inputs:    
        1. nodes - a numpy array (# Nodes x 6) of node coordinates  
                - the row number of a set of coordinates is the node's id  
                - initializing with an array of fewer than 6 columns fills the rest with 0s  
        2. xsection_list - A nested list of parameters defining each cross section (see XSection Class)  
        3. element_list - A nested list of parameters defining each element (see Element Class)  
        4. constraints - A nested list of the freedom of each DOF for each constrained node (1 for fixed, 0 for free)  
                - ie. for pinned 'node_id1' and fixed 'node_id2':  
                - constraints = [[node_id1, 1, 1, 1, 0, 0, 0], [node_id2, 1, 1, 1, 1, 1, 1]]  

    def bending_shape_functions(self, L, N_pts):
        # Calculate Hermite shape functions for bending displacement interpolation  
        Inputs:  
        1. L - length(s) of element(s), including end displacements  
        2. N_pt - Number of interpolation points  
        Evaluates the shape functions on a mesh defined by L and N_pts to allow returning a matrix of values for all elements at once

    def apply_load(self, applied_forces: list, N_pts: int = 10, scale = 1):  
        # Applies a load to the frame  
        Inputs:  
        1. applied_forces - nested list of applied forces, with a similar format to the constraint specification  
                -ie. applied_forces = [[node_id, Fx, Fy, Fz, Mx, My, Mz], ...]  
        2. N_pts - number of interpolation points along each element  
        3. scale - multiplies the applied forces to magnify the resulting displacements  
        Results:  
        Rewrites the frame property "deformation", a dictionary of results:  
        "disps" - numpy array (#Nodes x 6), displacement of all nodes  
        "forces" - numpy array (#Nodes x 6), forces at all nodes  
        "el_disps" - numpy array (#Elements x 12), local displacements at each end of each element  
        "el_forces" - numpy array (#Elements x 12), local forces at each end of each element  
        "inter_coords" - numpy array (N_pts x #Elements x 6), interpolated displacements due to applied forces  
        "crit_load_factor" - float, smallest positive eigenvalue representing approximate buckling load factor  
        "crit_load_vec" - numpy array (#Nodes x 6), primary buckling mode shape (given as displacement for each DOF)  
        "inter_coords_buckled" - numpy array (N_pts x #Elements x 6), interpolated displacements due to buckling mode displacments  

    def interpolate_displacements(self, node_disps, el_disps, N_pts):  
        # Interpolates displacements with bending shape functions and linear axial displacement  

    def subdivide(self, N_div: int):  
        # Returns new frame with elements divided into N_div subelements (with same cross sections)
        # Keeps old nodes for easy constraint and force specification, and appends new intermediate nodes

    def print_deformed_results(self):  
        # Using *Texttable*, prints an ascii table of all deformation results, for each node and for each element  

    def plot_deformed(self, disp_type = ""):  
        # Plots the deformed frame with a different color for each cross section  
        # Setting disp_type to "buckled" plots the primary buckling mode shape instead of the displacement due to the applied load  

**Element Class**  

    def __init__(self, parent_frame, node_a: int, node_b: int, xsec = 0, z_vec = []):  
    Inputs:  
    1. parent_frame - The frame instance that the element belongs to  
    2. node_a, node_b - The endpoints of the element, given as node ids  
    4. xsec - The element's cross section, given by its id, its instance, or a set of parameters to create a new instance  
    5. z_vec - The orientation vector of the element's cross section  
            - if not provided, assumed to be the normalized vector of cross(x_vec, [1,0,0]) (or cross(x_vec, [0,1,0]), if x_vec is nearly parallel to [1,0,0])  
    Results:  
    The element's coordination transformation matrix is assigned to property '.gamma', and its length to '.L' - note: gamma implemented more concisely than in provided math.utils  
    The element's elastic stiffness matrix is calculated and assigned to property '.Ke_local' in local coordinates, and '.Ke' in global coordinates  

    def calc_elastic_stiffness(self):  
        # Calculate elastic stiffness matrix - note: A more concise version than in provided math.utils  

    def calc_geometric_stiffness(self, forces: np.array):  
        # Calculate geometric stiffness matrix - note: A more concise version than in provided math.utils  

**XSection Class**

    def __init__(self, id: int, E: float, v: float, geometry_type: str, geometry_parameters: list):  
    Inputs:  
    1. id - a unique integer for the cross section (within a given frame)
    2. E - the modulus of elasticity
    3. v - nu
    4. geometry_type, geometry_parameters - the cross section shape and defining parameters
            - "circle", [radius]
            - "rectangle", [base, height, J (optional)]
            - other, [Area, I_y, I_z, I_rho, J]





---

#### **Tutorial/Examples**

The following examples are from problems posted for in-class code reviews in ENG ME700 (Spring 2025).
They can also be found in Direct_Stiffness/src/example_problems_script.py for convenient use.

##### 1.

Here is an example of how to create a frame, apply a load to it, and see the results.
This example corresponds to the first one presented in "Assignment 2 - Code Review 1 - Example Problems":  
![image](A2_ex1.png)

There are three nodes, two elements with given z-directions, and one cross section shared by all elements.

Each row of the 'nodes' array contains the coordinate of a node - the index of the row is the node's id.

A list of elements is created, where the parameters for each element are given as a list of two node ids, the cross section id, and the z-vector. Both the z-vector and the cross section can be omitted if the z-vector can be assumed and if there is only one cross section in the frame, respectively.

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
    print("\nCode Review 1 - Problem 1:\n")
    simple_frame.print_deformed_results()

*Output:  
![image](A2_ex1_outA.png)  
![image](A2_ex1_outB.png)
![image](A2_ex1_outC.png)  
![image](A2_ex1_outD.png)  

##### 2. 

This example corresponds to the first example presented in "Assignment 2 - Code Review 1 - Example Problems":
![image](A2_ex2.png)

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
    print("\nCode Review 1 - Problem 2:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed()

#### 3.

This example corresponds to the first example presented in "Assignment_2_Code_Review_Part_2_.pdf":
![image](A2_ex3.png)

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
    print("\nCode Review 2 - Problem 1:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed()

#### 4.

This example corresponds to the second example presented in "Assignment_2_Code_Review_Part_2_.pdf":
![image](A2_ex4.png)

    # Frame geometry definition
    (Lx, Ly, Lz) = (10, 20, 25)
    x = [0, Lx, Lx, 0, 0, Lx, Lx, 0]
    y = [0, 0, Ly, Ly, 0, 0, Ly, Ly]
    z = [0, 0, 0, 0, Lz, Lz, Lz, Lz]
    nodes = np.array([np.array([x[i], y[i], z[i], 0, 0, 0]) for i in range(0, 8)])
    elements = [[i, i+4] for i in range(0, 4)]
    elements.extend([4 + i, 4 + (i+1) % 4] for i in range(0, 4))

    # Cross section list
    (E, v) = (500, 0.3)
    r = 0.5
    xsection = [[E, v, 'circle', [r]]]

    # Constraint list (node_id, fixed DOF)
    constraints = [[i,1,1,1,1,1,1] for i in range(0, 4)]

    # Force list (node_id, forces on each DOF)
    forces = [[i,0,0,-1,0,0,0] for i in range(4, 8)]

    # Create frame, apply loads, and display results
    simple_frame = frame(nodes, xsection, elements, constraints)
    divided_frame = simple_frame.subdivide(2)
    divided_frame.apply_load(forces, 20, 5)
    print("\nCode Review 2 - Problem 2:\n")
    divided_frame.print_deformed_results()
    divided_frame.plot_deformed("buckled")

#### 5.

This example corresponds to the first and second problem presented in "Assignment_2_Final_Code_Review.pdf":
![image](A2_ex5.png)

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
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    print("\nTechnical Correctness 1 - Problem 1:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed()

    # Problem 2 - Force list (node_id, forces on each DOF)
    P = 1
    L = np.linalg.norm(np.array([25, 50, 37]))
    (Fx, Fy, Fz) = (-25*P/L, -50*P/L, -37*P/L)
    forces_2 = [[6, Fx, Fy, Fz, 0, 0, 0]]

    # Create frame, apply loads, and display results
    
    simple_frame.apply_load(forces_2, 30)
    print("\nTechnical Correctness 1 - Problem 2:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed()

#### 6.

This example corresponds to the third problem presented in "Assignment_2_Final_Code_Review.pdf":
![image](A2_ex6.png)

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
    simple_frame = frame(nodes, xsection, elements, constraints)
    simple_frame.apply_load(forces, 30)
    print("\nTechnical Correctness 1 - Problem 3:\n")
    simple_frame.print_deformed_results()
    simple_frame.plot_deformed("buckled")

---

### More information <a name="more"></a>
More information can be found here:
* https://digitalcommons.bucknell.edu/cgi/viewcontent.cgi?article=1006&context=books
* https://learnaboutstructures.com/Matrix-Structural-Analysis-Introduction
