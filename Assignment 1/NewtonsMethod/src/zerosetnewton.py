import numpy as np
import sympy as sp
from sympy.abc import x
from dataclasses import dataclass

# Class storing a vector-valued vector function, its Jacobian, and two points;
#   as well as the scalar-valued function defined as the square of the norm of the vector function, and its gradient;
#   as well as the parameterized function defined as the scalar-valued function evaluated on the line between the two points,
#   and the gradient projected onto the line between the two points (ie. the derivative of the parameterized function)
#     Purpose: To collapse a bracketed N-dimensional root finding problem onto one dimension
@dataclass
class ParametricFunction:

    def __init__(self, fun, pt_a: np.ndarray, pt_b: np.ndarray, J = None):
        self.fun = fun
        self.J = J
        self.pt_a = pt_a
        self.pt_b = pt_b
        self.to_scalar()
        self.parameterize()
    
    # Validate Instance Inputs
    @property
    def fun(self):
        return self._fun
    
    @fun.setter
    def fun(self, input_fun):
        if all(isinstance(element, sp.Expr) for element in input_fun):
            self.vars = list(input_fun.free_symbols)
            self.dim = len(self.vars)
            if self.dim == 0: raise Exception("Symbolic function has no independent variables")
        else: raise Exception("Function is not a symbolic expression")
        self._fun = input_fun

    @property
    def J(self):
        return self._J
    
    @J.setter
    def J(self, input_J):
        if input_J is None:
            self.J = self.fun.jacobian(self.vars)
        else:
            if all(isinstance(element, sp.Expr) for element in input_J):
                n_dim = len(input_J.free_symbols)
                if n_dim != self.dim: raise Exception("Jacobian has different number of variables than function")
            else: raise Exception("Jacobian is not a symbolic expression")
            self._J = input_J

    @property
    def pt_a(self):
        return self._pt_a
    
    @pt_a.setter
    def pt_a(self, input_pt_a):
        if len(input_pt_a) != self.dim: raise Exception("not enough coordinates given for point a")
        self._pt_a = input_pt_a

    @property
    def pt_b(self):
        return self._pt_b
    
    @pt_b.setter
    def pt_b(self, input_pt_b):
        if len(input_pt_b) != self.dim: raise Exception("not enough coordinates given for point b")
        if self.pt_a == input_pt_b: raise Exception("point a and point b must be different")
        self._pt_b = input_pt_b

    # Function to create a scalar function and its gradient from a given vector function and its Jacobian
    # The scalar function created is equal to the magnitude of the vector function squared (ie. ||fun||^2)
    def to_scalar(self):
        # if (isinstance(self.fun, sp.Basic)):
        self.fun_scalar = self.fun.dot(self.fun)
        self.grad = 2 * self.J.T * self.fun
        # else:
        #     self.fun_scalar = lambda x: np.dot(self.fun(x), self.fun(x))
        #     self.grad = lambda x: 2 * np.matmul(self.J(x).transpose(), self.fun(x))
        return self
    
    # Function to parameterize the scalar function and its derivative (the projected gradient) along the line from pt_a to pt_b
    # The scalar value x ranges from 0 to 1 and defines the position in space as (pt_a + x*(pt_b-pt_a))
    def parameterize(self):
        vector = self.pt_b - self.pt_a
        # if (isinstance(self.fun, sp.Basic)):
        #sym_param = self.fun_scalar.subs(
        self.fun_param = sp.lambdify([x], self.fun_scalar(self.pt_a + x*vector), "numpy")
        self.grad_param = sp.lambdify([x], self.grad(self.pt_a + x*vector))
        # else:
        #     self.fun_param = lambda x: self.fun_scalar(self.pt_a + x*vector)
        #     self.grad_param = lambda x: np.dot(vector, self.grad(self.pt_a + x*vector))
        return self

@dataclass
class ZeroSetNewton:
    tol: float
    max_iter: int

    def __init__(self, fun, fun_d1, tol: float = 10**-12, max_iter: int = 10**3):
        self.fun = fun
        self.fun_d1 = fun_d1
        self.collapsed_fun = None
        self.collapsed_J = None
        self.tol = tol
        self.max_iter = max_iter
        self.num_iter = 0

    # Validate Instance Inputs
    @property
    def fun(self):
        return self._fun
    
    @fun.setter
    def fun(self, input_fun):
        if not callable(input_fun): raise Exception("Function to evaluate must be a function")
        if not (input_fun.__code__.co_argcount == 1): raise Exception("Function to evaluate must take only 1 argument")
        self._fun = input_fun

    @property
    def fun_d1(self):
        return self._fun_d1
    
    @fun_d1.setter
    def fun_d1(self, input_fun_d1):
        if not callable(input_fun_d1): raise Exception("Function to evaluate must be a function")
        if not (input_fun_d1.__code__.co_argcount == 1): raise Exception("Function to evaluate must take only 1 argument")
        self._fun_d1 = input_fun_d1

    @property
    def tol(self):
        return self._tol

    @tol.setter
    def tol(self, input_tol):
        if input_tol <= 0: raise Exception("Tolerance must be greater than zero")
        self._tol = input_tol

    @property
    def max_iter(self):
        return self._max_iter

    @max_iter.setter
    def maxiter(self, input_max_iter):
        if input_max_iter <= 0: raise Exception("Max number of iterations must be greater than zero")
        self._max_iter = input_max_iter

    def general_newton_nd(self):
        return self
    
    def bracketed_newton_1d(self):
        return self
    
    def fit_quadtratic(self):
        return self
    

    

    # Calculate coordinate of zero as if evaluated function were linear
    # -- Used to partition interval, instead of midpoint as in the standard bisection method
    # y(x) = slope * (x-x_0) + y(x_0)
    # y(x) = (pt_vals[1]-pt_vals[0])/(pt_coords[1]-pt_coords[0]) * (x - pt_coords[0]) + pt_vals[0]
    # calculate x so that y(x) = 0
    def calc_zero(self):
        scaledX = - self.vals[0] / (self.vals[1]-self.vals[0])
        if scaledX < self.xBuffer:
          scaledX = self.xBuffer
        elif scaledX > 1 - self.xBuffer:
          scaledX = 1 - self.xBuffer
        return self.xDomain[0] + scaledX * (self.xDomain[1]-self.xDomain[0])
    
    def __iter__(self):
        return self
    
    def __next__(self):
        self.num_iter = self.num_iter + 1
        self.pt = self.calcMidpoint()
        try:
            self.ptVal = self.f(self.pt)
        except:
            print("Function could not be evaluated at x=", self.pt)
            raise StopIteration
        if abs(self.ptVal) <= self.tol or self.num_iter > self.max_iter:
            raise StopIteration
        if self.ptVal*np.sign(self.vals[1]) < 0:
            self.xDomain[0] = self.pt
            self.vals[0] = self.ptVal
        else:
            self.xDomain[1] = self.pt
            self.vals[1] = self.ptVal
        return self