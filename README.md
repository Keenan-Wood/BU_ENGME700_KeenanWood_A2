# Matrix Structural Analysis

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1)
[![tests](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions/workflows/tests.yml/badge.svg)](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions)

## Apologies, the plotting functionallity for part 2 is still untested
### Please see the example below (under Core Implementation) for how a problem should be set up, and how applying a load should be done
### Advice on interpolating displacement (see functions in StructuralFrame_2.py) is very appreciated - thanks!
---

<details>
    <summary>1. Core Implementation</summary>
    
# Core Implementation

### Table of Contents
* [The Method](#algo)
* [Conda environment, installation, and testing](#install)
* [Documentation & Examples](#tutorial)
* [More Information](#more)

---

### Matrix Structural Analysis <a name="algo"></a>

**Matrix Structural Analysis** (To be written)

---

### Conda environment, install, and testing <a name="install"></a>

To install this package, please begin by setting up a conda environment and activating it. For example:
```bash
conda create --name me700-env python=3.12
conda activate me700-env
```

Navigate to the project directory (./Part_1) and create an editable install of the code:
```bash
pip install -e .
```

Test that the code is working with pytest:
```bash
pytest -v --cov=newtonmethod --cov-report term-missing
```

If you are using VSCode to run this code, don't forget to set VSCode virtual environment to the newly-activated environment.

---

#### **Documentation**

**load_frame function**
*(all_disps, all_forces, crit_factor, crit_vec) = load_frame(nodes, elements, xsection, constraints, forces)*

Inputs:
1. nodes - 2D numpy array of node coordinates (#Nodes x 6)
    ie. for two nodes n0 and n1:
    nodes = np.array([[n0x, n0y, n0z, n0tx, n0ty, n0tz], [n1x, n1y, n1z, n1tx, n1ty, n1tz]])
    where tx, ty, and tz represent the angular coordinates of the nodes (typically 0 to start)

2. elements - nested list containing element information (#Elements x 4)
    ie. for two elements el0 and el1:
    elements = [[el0_node_a_id, el0_node_b_id, el0_section_id, el0_zvec], [el1_node_a_id, el1_node_b_id, el1_section_id, el1_zvec]]

3. xsection - nested list of section properties (#Different Sections x 7)
    xsection = [[E, A, I_y, I_z, I_rho, J, nu]] creates a section with ID 0 with the given properties

4. constraints - nested list of fixed DOF (#Constrained Nodes x 7)
    ie. for node 0 with fixed z, node 3 completely fixed, and node 4 pinned:
    constraints = [[0,0,0,1,0,0,0], [3,1,1,1,1,1,1], [4,1,1,1,0,0,0]]

5. forces - nested list of forces on each DOF for indicated nodes (#Forced Nodes x 7)
    ie. forces and moments applied to node 1:
    forces = [[1, -0.05, 0.075, 0.1, -0.05, 0.1, -0.25]]

---

### Tutorial <a name="tutorial"></a>

---

#### **Examples**

##### 1.

Here is an example of how to setup and apply a load to a frame.
This example corresponds to the first one presented in "Assignment 2 - Code Review 1 - Example Problems":
![image](A2_ex1.png)

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

##### 2. 

This example corresponds to the first example presented in "Assignment 2 - Code Review 1 - Example Problems":
![image](A2_ex2.png)

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

#### 3.

This example corresponds to the first example presented in "Assignment_2_Code_Review_Part_2_.pdf":
![image](A2_ex3.png)

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

---

### More information <a name="more"></a>
More information can be found here:
* https://learnaboutstructures.com/Matrix-Structural-Analysis-Introduction

</details>


<details>
    <summary>2. Visualization</summary>

# Visualization

</details>
