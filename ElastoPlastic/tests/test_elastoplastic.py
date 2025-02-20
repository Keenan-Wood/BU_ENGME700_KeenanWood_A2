import numpy as np
import matplotlib.pyplot as plt
from elastoplastic import material
from elastoplastic import ElastoPlastic

# Misc. additional material properties
# alum = material('aluminum', 70, 0.07, 0.095)
# silicon_carbide = material('silicon_carbide', 450, 4.50, 3.440)
# copper = material('copper', 117, 1.17, 0.070)

stress_tol = 10**-6
strain_tol = 10**-6

def test_isotropic_steel():
    steel = material('steel', 210, 210/9, 0.250)
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 1, 0)
    correct_strain = abs(test_steel.strain + .63) < strain_tol
    correct_stress = abs(test_steel.stress + 18.5310144) < stress_tol
    correct_yield_stress = abs(test_steel.yield_stress - 18.5310144) < stress_tol
    correct_back_stress = abs(test_steel.back_stress - 0) < stress_tol
    assert correct_strain and correct_stress and correct_yield_stress and correct_back_stress

def test_kinematic_steel():
    steel = material('steel', 210, 210/9, 0.250)
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 0, 1)
    correct_strain = abs(test_steel.strain + .63) < strain_tol
    correct_stress = abs(test_steel.stress + 13.455) < stress_tol
    correct_yield_stress = abs(test_steel.yield_stress - .25) < stress_tol
    correct_back_stress = abs(test_steel.back_stress + 13.205) < stress_tol
    assert correct_strain and correct_stress and correct_yield_stress and correct_back_stress

def test_isotropic_kinematic_steel():
    steel = material('steel', 210, 210/9, 0.250)
    test_steel = ElastoPlastic(steel, 0, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_steel.stretch(set_strain, 1, 1)
    correct_strain = abs(test_steel.strain + .63) < strain_tol
    correct_stress = abs(test_steel.stress + 17.7433163415) < stress_tol
    correct_yield_stress = abs(test_steel.yield_stress - 18.8656619765) < stress_tol
    correct_back_stress = abs(test_steel.back_stress + 12.7285204065) < stress_tol
    assert correct_strain and correct_stress and correct_yield_stress and correct_back_stress

def test_isotropic_nylon_prestrain():
    nylon = material('nylon6', 3, 3/10, 0.045)
    test_nylon = ElastoPlastic(nylon, .3, 0)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_nylon.stretch(set_strain, 1, 0)
    correct_strain = abs(test_nylon.strain + .33) < strain_tol
    correct_stress = abs(test_nylon.stress + .2577037087630625) < stress_tol
    correct_yield_stress = abs(test_nylon.yield_stress - .2577037087630626) < stress_tol
    correct_back_stress = abs(test_nylon.back_stress - 0) < stress_tol
    assert correct_strain and correct_stress and correct_yield_stress and correct_back_stress

def test_kinematic_nylon_backstress():
    nylon = material('nylon6', 3, 3/10, 0.045)
    test_nylon = ElastoPlastic(nylon, .3, 100)
    set_strain = np.array([.01, .01, .01, -.03, .05, -.08, .1, -.7])
    test_nylon.stretch(set_strain, 0, 1)
    correct_strain = abs(test_nylon.strain + .33) < strain_tol
    correct_stress = abs(test_nylon.stress - 8.878181818181819) < stress_tol
    correct_yield_stress = abs(test_nylon.yield_stress - .045) < stress_tol
    correct_back_stress = abs(test_nylon.back_stress - 8.923181818181817) < stress_tol
    assert correct_strain and correct_stress and correct_yield_stress and correct_back_stress