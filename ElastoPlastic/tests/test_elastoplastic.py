import numpy as np
import matplotlib.pyplot as plt
from elastoplastic import material
from elastoplastic import ElastoPlastic

steel = material('steel', 210, 2.10, 0.250)
alum = material('aluminum', 70, 0.07, 0.095)
nylon = material('nylon6', 3, 0.003, 0.045)
silicon_carbide = material('silicon_carbide', 450, 4.50, 3.440)
copper = material('copper', 117, 1.17, 0.070)

def test_isotropic_steel():
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 1, 0)
    assert true

def test_kinematic_steel():
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 1, 0)
    assert true

def test_isotropic_kinematic_steel():
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 1, 1)
    assert true

def test_isotropic_nylon_prestrain():
    test_nylon = ElastoPlastic(nylon, .3, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_nylon.stretch(set_strain, 1, 0)
    assert true

def test_kinematic_nylon_backstress():
    test_nylon = ElastoPlastic(nylon, .3, 100)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_nylon.stretch(set_strain, 0, 1)
    assert true