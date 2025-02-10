# Introductory Numerical Methods

<details>
    <summary>1. Bisection Method</summary>
    
# Bisection Method Implementation

</details>

<details>
    <summary>2. Newton's Method</summary>

# Newton Method Implementation

[![python](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
![os](https://img.shields.io/badge/os-ubuntu%20|%20macos%20|%20windows-blue.svg)
[![license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/sandialabs/sibl#license)

[![codecov](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/graph/badge.svg?token=p5DMvJ6byO)](https://codecov.io/gh/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod)
[![tests](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/actions/workflows/tests.yml/badge.svg)](https://github.com/Keenan-Wood/BU_ENGME700_KeenanWood_A1_NewtonMethod/actions)

---

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
