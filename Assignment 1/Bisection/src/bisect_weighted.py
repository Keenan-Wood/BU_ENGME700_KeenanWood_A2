import numpy as np
import math
from dataclasses import dataclass

@dataclass
class BisInterval:
    xDomain: np.ndarray
    TOL: float
    nIteMax: int
    xBuffer: float

    def __init__(self, xDomain: np.ndarray, f, TOL: float = 10**-12, nIteMax: int = 10**3, xBuffer: float = .02):
        self.xDomain = xDomain
        self.f = f
        self.TOL = TOL
        self.nIteMax = nIteMax
        self.xBuffer = xBuffer
        self.vals = np.array([f(self.xDomain[0]), f(self.xDomain[1])])
        self.pt = sum(xDomain)/2
        self.ptVal = f(self.pt)
        self.nIte = 0

    # Validate Instance Inputs
    @property
    def xDomain(self):
        return self._xDomain

    @xDomain.setter
    def xDomain(self, input_xDomain):
        if not (len(input_xDomain) == 2): raise Exception("Specified domain must have two values")
        self._xDomain = input_xDomain.astype(np.float64)

    @property
    def f(self):
        return self._f
    
    @f.setter
    def f(self, input_f):
        if not callable(input_f): raise Exception("Function to evaluate must be a function")
        if not (input_f.__code__.co_argcount == 1): raise Exception("Function to evaluate must take only 1 argument")
        try: SAMESIGN = np.sign(input_f(self.xDomain[0]))*np.sign(input_f(self.xDomain[1])) >= 0
        except: raise Exception("Function could not be evaluated at bounds of provided domain")
        if SAMESIGN: raise Exception("The function should have opposite signs when evaluated at each boundary")
        self._f = input_f

    @property
    def TOL(self):
        return self._TOL

    @TOL.setter
    def TOL(self, input_TOL):
        if input_TOL <= 0: raise Exception("Tolerance must be greater than zero")
        self._TOL = input_TOL

    @property
    def nIteMax(self):
        return self._nIteMax

    @nIteMax.setter
    def nIteMax(self, input_nIteMax):
        if input_nIteMax <= 0: raise Exception("Max number of iterations must be greater than zero")
        self._nIteMax = input_nIteMax

    @property
    def xBuffer(self):
        return self._xBuffer

    @xBuffer.setter
    def xBuffer(self, input_xBuffer):
        if input_xBuffer < 0: raise Exception("xBuffer must be >= 0")
        if input_xBuffer >= .5: raise Exception("xBuffer must be < 0.5")
        self._xBuffer = input_xBuffer

    # Calculate coordinate of zero as if evaluated function were linear
    # -- Used to partition interval, instead of midpoint as in the standard bisection method
    # y(x) = slope * (x-x_0) + y(x_0)
    # y(x) = (pt_vals[1]-pt_vals[0])/(pt_coords[1]-pt_coords[0]) * (x - pt_coords[0]) + pt_vals[0]
    # calculate x so that y(x) = 0
    def calcMidpoint(self):
        scaledX = - self.vals[0] / (self.vals[1]-self.vals[0])
        if scaledX < self.xBuffer:
          scaledX = self.xBuffer
        elif scaledX > 1 - self.xBuffer:
          scaledX = 1 - self.xBuffer
        return self.xDomain[0] + scaledX * (self.xDomain[1]-self.xDomain[0])
    
    def __iter__(self):
        return self
    
    def __next__(self):
        self.nIte = self.nIte + 1
        self.pt = self.calcMidpoint()
        try:
            self.ptVal = self.f(self.pt)
        except:
            print("Function could not be evaluated at x=", self.pt)
            raise StopIteration
        if abs(self.ptVal) <= self.TOL or self.nIte > self.nIteMax:
            raise StopIteration
        if self.ptVal*np.sign(self.vals[1]) < 0:
            self.xDomain[0] = self.pt
            self.vals[0] = self.ptVal
        else:
            self.xDomain[1] = self.pt
            self.vals[1] = self.ptVal
        return self