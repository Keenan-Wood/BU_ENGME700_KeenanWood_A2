import numpy as np
from bisect_weighted import BisInterval

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

#f = lambda x: 0
#for n in range(1,200):
#    f = lambda x: f(x) + np.sin(2*math.pi*n**2*x) / n**2
#xBounds = np.array([math.pi/2+.1, 3])
#bisectInterval = BisInterval(xBounds, f)


subInterval = None
for subInterval in bisectInterval:
    pass
if not subInterval is None: 
    print("x=", subInterval.pt, ";", "nIte=", subInterval.nIte)
