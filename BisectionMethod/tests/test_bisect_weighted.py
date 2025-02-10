import numpy as np
import math
from bisect_weighted import BisInterval

def test_simple_function():
    f = lambda x: x - .5
    xBounds = np.array([-1,2])
    bisectInterval = BisInterval(xBounds, f)
    assert abs(bisectInterval.pt - .5) < 10**-8

#def test_Riemann():
#    f = lambda x: 0
#    for n in range(1,200):
#        f = lambda x: f(x) + np.sin(2*math.pi*n**2*x) / n**2
#    xBounds = np.array([math.pi/2, 3])
#    bisectInterval = BisInterval(xBounds, f)
#    assert abs(bisectInterval.pt) < 10**-8
