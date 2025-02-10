import numpy as np
import sympy as sp
from parametricfunction import ParametricFunction

def test_simple_function():
    x, y = sp.symbols('x y')
    fun_expr = sp.Matrix([x + y, x - y])
    fun_vars = [x, y]
    pt_a = np.array([1,1])
    pt_b = np.array([2,2])
    param_fun = ParametricFunction(fun_expr, fun_vars, pt_a, pt_b)
    param_fun_pt_a = param_fun.fun_param(0)
    test_fun_pt_a_vec = fun_expr.subs([(x, 1), (y, 1)]).n()
    test_fun_pt_a = sum(element**2 for element in test_fun_pt_a_vec)
    assert abs(param_fun_pt_a - test_fun_pt_a) < 10**-10

test_simple_function()
v = 2