import numpy as np
import sympy as sp
from src.zerosetnewton import ParametricFunction

def test_simple_function():
    x, y = sp.symbols('x y')
    fun_expr = [x, y]
    pt_a = np.ndarray([1,1])
    pt_b = np.ndarray([2,2])
    param_fun = ParametricFunction(fun_expr, pt_a, pt_b)
    assert param_fun.fun_scalar == 1