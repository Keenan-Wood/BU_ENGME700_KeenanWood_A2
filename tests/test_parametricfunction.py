import numpy as np
import sympy as sp
from zerosetnewton import ParametricFunction

def test_simple_function():
    x, y = sp.symbols('x y')
    fun_expr = sp.Matrix([x + y, x - y])
    fun_vars = [x, y]
    pt_a = np.array([1,1])
    pt_b = np.array([2,2])
    param_fun = ParametricFunction(fun_expr, fun_vars, pt_a, pt_b)
    assert param_fun.fun_scalar == 1

test_simple_function()
v = 2