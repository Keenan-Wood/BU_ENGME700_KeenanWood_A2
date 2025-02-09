import numpy as np
import sympy as sp
from dataclasses import dataclass

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