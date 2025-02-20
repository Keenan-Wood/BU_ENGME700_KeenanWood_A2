import numpy as np
import sympy as sp
import pytest
from newtonmethod import NewtonMethod

def test_not_symbolic_jacobian():
    with pytest.raises(Exception) as e_info:
        x, y = sp.symbols('x y')
        fun_expr = sp.Matrix([x + y, x - y])
        fun_jacobian = lambda a, b: [[1, 1], [1, -1]]
        fun_vars = [x, y]
        start_pt = np.array([1,1])
        newton = NewtonMethod(fun_expr, fun_vars, start_pt, fun_jacobian)


def test_singular_jacobian():
    x, y = sp.symbols('x y')
    fun_expr = sp.Matrix([x + y, x - y])
    fun_vars = [x, y]
    start_pt = np.array([1,1])
    newton = NewtonMethod(fun_expr, fun_vars, start_pt)
    subNewton = None
    for subNewton in newton:
        pass
    if not subNewton is None: 
        print("x=", subNewton.pt, ";", "# of iterations=", subNewton.num_iter)
    assert subNewton is None

def test_unevaluatable_function():
    x, y = sp.symbols('x y')
    fun_expr = sp.Matrix([1/x + y, x - y])
    fun_vars = [x, y]
    start_pt = np.array([0,1])
    newton = NewtonMethod(fun_expr, fun_vars, start_pt)
    subNewton = None
    for subNewton in newton:
        pass
    if not subNewton is None: 
        print("x=", subNewton.pt, ";", "# of iterations=", subNewton.num_iter)
    assert subNewton is None

def test_simple_1D_fun():
    x = sp.symbols('x')
    fun_expr = sp.Matrix([(x-1)**2])
    J_expr = sp.Matrix([2*(x-1)])
    fun_vars = [x]
    start_pt = np.array([7])
    newton = NewtonMethod(fun_expr, fun_vars, start_pt, J_expr)
    subNewton = None
    for subNewton in newton:
        pass
    if not subNewton is None: 
        print("x=", subNewton.pt, ";", "# of iterations=", subNewton.num_iter)
    assert abs(subNewton.ptVal) < 10**-12

def test_simple_3D_fun():
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
    assert np.linalg.norm(subNewton.ptVal) < 10**-10