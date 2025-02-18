# Introductory Numerical Methods

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod)
[![tests](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/actions/workflows/tests.yml/badge.svg)](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/actions)

---

<details>
    <summary>1. Bisection Method</summary>
    
# Bisection Method Implementation

</details>

<details>
    <summary>2. Newton's Method</summary>

# Newton's Method Implementation

### Table of Contents
* [Getting Started](#gs)
* [Newton Method Algorithm](#algo)
* [Conda environment, installation, and testing](#install)
* [Tutorial](#tutorial)
* [More Information](#more)

---

### Getting Started

To be written

---

### Newton Method Algorithm <a name="algo"></a>

**Newton's Method** is a numerical technique to find roots of a continuous function f(x) whose jacobian **$J$** is continuous. Given an initial point **$x_0$** and its resultant **$R(x_0) = f(x_0)$**, Newton's method generates a more accurate estimate of the zero of f, **$x_1 = x_0 - J(x_0)^{-1} R(x_0)$**. Iteration produces a sequence of positions which for most well-behaved functions converges quadtratically to a root of f.

**Advantages of Newton's Method**:
1. **Fast**: The method in most cases converges quadtratically.
2. **Efficient Evaluation**: Higher-order derivatives of f, which may be expensive to evaluate, do not need to be evaluated.
3. **Robustness**: It works well for a wide range of functions.

**Limitations of Newton's Method**:
1. **Non-convergence**: Certain combinations of functions and initial points can diverge or oscillate ad infinitum.
2. **Identifies single root**: If convergent, only one root of the function is identified; In particular, if f takes an N-dimensional input and outputs a P-dimensional vector, the zero set of f typically has dimension N-P.
3. **Function Continuity and Differentiability Required**: f and its jacobian must be continuous.
4. **Unbounded domain**: In its simplest implementation, Newton's method does not utilize information on domain bounds to improve convergence.

---

### Conda environment, install, and testing <a name="install"></a>

To install this package, please begin by setting up a conda environment and activating it. For example:
```bash
conda create --name newton-method-env python=3.12
conda activate newton-method-env
```

Navigate to the project directory and create an editable install of the code:
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

#### **What Does the Function Do?**

The `NewtonMethod` class instantiates with:
- A symbolic expression, or array of expressions, **fun**
- A list of the symbolic variables, **vars_indep**
- A starting point (numpy array), **start_pt**
- The jacobian, as a matrix of symbolic expressions, **J**: Default **None**
- The solver tolerance (on the Residual) as a float, **tol**: Default **$10**-12$**
- The maximum number of iterations as an integer, **max_iter**: Default **$10**3$**

If no Jacobian is provided, it is calculated during initialization.

The created object has **.pt** and **.ptVal** properties, representing the current estimate of the zero's position, and the value of the residual, respectively.

Iterating over the initialized object until finished results in a final object which either converged, with the **.pt** attribute giving the position of the root, or failed to converge, in which case the number of iterations given by the property **.num_iter** will equal the property **.max_iter**.

---

### **Summary of Errors and Their Causes**

To be written

---

#### **Examples**

    x, y, z = sp.symbols('x y z')
    fun_expr = sp.Matrix([1 + x + y, x**2 - y**3 + z, x*y - z])
    fun_vars = [x, y, z]
    start_pt = np.array([1,1,1])
    newton = NewtonMethod(fun_expr, fun_vars, start_pt)
    subNewton = None
    for subNewton in newton:
        pass
    if not subNewton is None: 
        print("x=", subNewton.pt, ";", "# of iterations=", subNewton.num_iter)

   Output:
   x= [-2.32471796  1.32471796 -3.07959562] ; # of iterations= 12

---

### More information <a name="more"></a>
More information can be found here:
* https://en.wikipedia.org/wiki/Newton%27s_method

</details>

<details>
    <summary>3. Elasto-Plastic Strain Hardening Simulation</summary>

# Elasto-Plastic Strain Hardening Simulation

### Table of Contents
* [Elasto-Plastic Strain Hardening with Predictor-Corrector](#algo)
* [Conda environment, installation, and testing](#install)
* [Tutorial](#tutorial)
* [More Information](#more)

---

### Elasto-Plastic Strain Hardening with Predictor-Corrector <a name="algo"></a>

The **Predictor-Corrector** approach consists in using data to approximate a solution, and then using that same data to add a correcting term to the approximation to get a much better approximation. Its primary benefit is to reduce the memory and computation needed for large problems that would otherwise require the calculation of two potentially large and computationally-intensive sets of data to calculate a good approximation, instead of the one.

The elasto-plastic model considers materials to be elastic until they reach their yield stress, at which point plastic flow starts occuring. In conjuction with the isotropic hardening model, the material's yield surface expands under plastic loading, while in the kinematic hardening model, the yield surface remains the same size, but its center, defined by a *back-stress*, shifts in the direction of the applied load.

The module presented here establishes a *material* and an *ElastoPlastic* class. The *material* class holds the basic material properties (the elastic and plastic moduli, and the yield strength before loading), and the ElastoPlastic class generates an object on which a *stretch()* method can be calculate its stress and strain after applying a series of strain increments. Two parameters to the stretch function determine the hardening behavior - if either is set to 0 and the other to 1, then the material will harden either fully isotropically or kinematically. The values can be adjusted to calculate the results of mixed hardening.

Instead of checking if the difference of the stress and the yield stress is negative to see if the material is in the elastic regime, the *stretch()* method takes the *min()* of 0 and the difference, so that no branching is necessary (ie. the value being 0 instead of negative naturally calculates elastic behavior).

# To-Do
- Debug pytest's file-finding/system path issue (see errors section)
- Add matplotlib functionallity for easy data visualization
- Verify the results with commercial solvers
- Expand the tutorial examples with rich descriptions and examples of how small changes effect the results

---

### Conda environment, install, and testing <a name="install"></a>

To install this package, please begin by setting up a conda environment and activating it. For example:
```bash
conda create --name elastoplastic-env python=3.12
conda activate elastoplastic-env
```

Navigate to the project directory (the *ElastoPlastic* folder) and create an editable install of the code:
```bash
pip install -e .
```

Test that the code is working with pytest:
```bash
pytest -v --cov=elastoplastic --cov-report term-missing
```

If you are using VSCode to run this code, don't forget to set VSCode virtual environment to the newly-activated environment.

---

### Tutorial <a name="tutorial"></a>

#### **Class Structure**

The `Material` class instantiates with:
- A name, as a string
- The elastic modulus, as a float
- The plastic modulus, as a float
- The yield strength, as a float

*Note*: Units are not currently supported - convert all values to consistent units for accurate results.

The `ElastoPlastic` class instantiates with:
- A material (defined above)
- The current strain, as a float
- The current stress, as a float
- The current back stress, as a float

Once instantiated, an ElastoPlastic object can be acted upon with the **stretch()** function, with arguments:
- The strain increments to apply, as a numpy array
- A scaling parameter relating to the isotropic hardening (0 for none, 1 for full)
- A scaling parameter relating to the kinematic hardening (0 for none, 1 for full)

---

### **Summary of Errors and Their Causes**

Known issue: pytest does not recognize elastoplastic module (likely system path/reference issue)

---

#### **Examples**

# Steel - Isotropic Hardening
    steel = material('steel', 210, 2.10, 0.250)
    isotropic_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    isotropic_steel.stretch(set_strain, 1, 0)

# Aluminum - Kinematic Hardening
    alum = material('aluminum', 70, 0.07, 0.095)
    kinematic_alum = ElastoPlastic(alum, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    kinematic_alum.stretch(set_strain, 1, 0)

# Copper - Isotropic and Kinematic Hardening
    copper = material('copper', 117, 1.17, 0.070)
    elastoplastic_copper = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    elastoplastic_copper.stretch(set_strain, 1, 1)

# Pre-strained Nylon - Isotropic Hardening
    nylon = material('nylon6', 3, 0.003, 0.045)
    isotropic_nylon = ElastoPlastic(nylon, .3, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    isotropic_nylon.stretch(set_strain, 1, 0)

# Silicon-Carbide - Kinematic Hardening with Back Stress
    silicon_carbide = material('silicon_carbide', 450, 4.50, 3.440)
    kinematic_carbide = ElastoPlastic(silicon_carbide, .3, 100)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    kinematic_carbide.stretch(set_strain, 0, 1)

---

### More information <a name="more"></a>
More information can be found here:
* https://innovationspace.ansys.com/courses/wp-content/uploads/sites/5/2020/12/Lesson-3-Hardening-of-Plasticity.pdf

</details>