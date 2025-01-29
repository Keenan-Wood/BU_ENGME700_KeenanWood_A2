import numpy as np
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

## Ex. 1 - Simple cubic
f = lambda x: (x-.101)**3
xBounds = np.array([-.1,1])
bisectInterval = BisInterval(xBounds, f)

## Ex. 2 - Spring problem from lecture 1
#k = 1
#L = 1
#F = 0.25
#f = lambda x: 2*k*x*(np.sqrt(L**2 + x**2) - L)/np.sqrt(L**2 + x**2) - F
#xBounds = np.array([0,2])
#bisectInterval = BisInterval(xBounds, f, 10**-9, 100, .01)

## Ex. 3 - Pendulum and spring
    # Consider a pendulum with length L and point mass M
    # Attach a spring (constant k, rest length R) to the point mass and to a fixed point (a,b)
    # For pendulum angle "t" at equilibrium:
    # cA = ((a-L*np.sin(t))*L*np.sin(t) + (b-L*np.cos(t))*L*np.cos(t)) / (np.sqrt((a-L*np.sin(t))**2 + (b-L*np.cos(t))**2) * L)
    # F_tangent = 0 = M*g*np.sin(t) + k*(np.sqrt( (a-L*np.sin(t))**2 +  (b-L*np.cos(t))**2 ) - R) * (1-cA**2)
#L = 1.0
#M = 1.0
#a = 1.5
#b = 1.5
#R = 1.0
#k = 10
#g = 9.8
#cA = lambda t: ((a-L*np.sin(t))*L*np.sin(t) + (b-L*np.cos(t))*L*np.cos(t)) / (np.sqrt((a-L*np.sin(t))**2 + (b-L*np.cos(t))**2) * L)
#f = lambda t: -M*g*np.sin(t) + k*(np.sqrt( (a-L*np.sin(t))**2 +  (b-L*np.cos(t))**2 ) - R) * (1-cA(t)**2)
#xBounds = np.array([0,1.6])
#bisectInterval = BisInterval(xBounds, f)

## Ex. 4 - Transcendental equation
#f = lambda x: np.cos(x) - x
#xBounds = np.array([0,2])
#bisectInterval = BisInterval(xBounds, f)

## Ex. 5 - Trig function - many roots
f = lambda x: np.sin(x)
xBounds = np.array([-.1,6.4])
bisectInterval = BisInterval(xBounds, f)


subInterval = None
for subInterval in bisectInterval:
    pass
if not subInterval is None: 
    print("x=", subInterval.pt, ";", "nIte=", subInterval.nIte)