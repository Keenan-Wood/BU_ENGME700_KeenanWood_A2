# Matrix Structural Analysis

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1)
[![tests](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions/workflows/tests.yml/badge.svg)](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1/actions)

## Apologies, the code has not been successfully tested yet
### Please see the example below for how a problem should be set up, and how applying a load should be done
### Advice on code (see StructuralFrame.py and test_frame.py) is very appreciated - thanks!
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

### Tutorial <a name="tutorial"></a>

#### Documentation**


---

#### **Examples**

##### 1.

Here is an example of how to setup and apply a load to a frame.
This example corresponds to the first one presented in "Assignment 2 - Code Review 1 - Example Problems":
![image](A2_ex1.png)

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


##### 2. 

---

### More information <a name="more"></a>
More information can be found here:
* https://learnaboutstructures.com/Matrix-Structural-Analysis-Introduction

</details>


<details>
    <summary>2. Visualization</summary>

# Visualization

</details>
