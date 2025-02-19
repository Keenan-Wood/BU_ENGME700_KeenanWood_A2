import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt

@dataclass
class material:
    E: float
    H: float
    yield_stress: float

    def __init__(self, name: str, E: float, H: float, yield_stress: float):
        self.name = name
        self.E = E
        self.H = H
        self.yield_stress = yield_stress

@dataclass
class ElastoPlastic:

    def __init__(self, mat: material, strain: float, stress: float, back_stress: float = 0):
        self.mat = mat
        self.yield_stress = self.mat.yield_stress
        self.strain = strain
        self.stress = stress
        self.back_stress = back_stress

    # Validate Instance Inputs
 #   @property
 #   def set_strain(self):
 #       return self._set_strain
    
 #   @set_strain.setter
 #   def set_strain(self, input_set_strain):
 #       if len(input_set_strain) < 1: raise Exception("Array of strain increments should have at least one element")
 #       self._set_strain = input_set_strain

    def stretch(self, set_strain: np.array, hard_iso: float, hard_kin: float):
        # Validate inputs
        if set_strain.size < 1: raise Exception("Array of strain increments should have at least one element")
        if hard_iso < 0 or hard_iso > 1: raise Exception("Isotropic hardening parameter must be between 0 (off) and 1 (on)")
        if hard_kin < 0 or hard_kin > 1: raise Exception("Kinematic hardening parameter must be between 0 (off) and 1 (on)")
        
        if set_strain.size == 1: set_strain = [set_strain]
        for strain_incr in set_strain:
            self.strain = self.strain + strain_incr

            # Update stress, back stress, and yield surface (Predictor-Corrector in 1 step using min())
            elastic_stress = self.stress + self.mat.E*strain_incr
            phi = max(0, abs(elastic_stress - self.back_stress) - self.yield_stress)
            plastic_strain = phi / (self.mat.E + self.mat.H)
            signed_strain = plastic_strain * np.sign(elastic_stress - self.back_stress)
            self.stress = elastic_stress - self.mat.E * signed_strain

            # Isotropic Part - Von-Mises Yield Criterion -> d(yield_stress) = 3/2 * h * d(plastic_strain)
            self.yield_stress = self.yield_stress + hard_iso * self.mat.H * plastic_strain
            # Kinematic Part - Prager's Rule -> d(back_stress) = c * d(plastic_strain)
            self.back_stress = self.back_stress + hard_kin * self.mat.H * signed_strain

            ## Compare with lecture notes - confusing update steps

## Test/Debugging script
if 1 == 0:
    alum = material('aluminum', 70, 70/9, 0.095)
    elasto_alum = ElastoPlastic(alum, 0, 0)
    # strain_incr_array = np.concatenate((np.ones(100) * .002/100, -np.ones(100) * .002/100), axis=0)
    #strain_incr_array = np.ones(100) * .002/100
    strain_incr_array = np.array([.002])
    N_steps = strain_incr_array.size
    strain = np.zeros(N_steps + 1)
    stress = np.zeros(N_steps + 1)
    back_stress = np.zeros(N_steps + 1)
    yield_stress = np.zeros(N_steps + 1)
    strain[0] = elasto_alum.strain
    stress[0] = elasto_alum.stress
    back_stress[0] = elasto_alum.back_stress
    yield_stress[0] = elasto_alum.yield_stress
    for i_step in range(1, N_steps + 1):
        elasto_alum.stretch(strain_incr_array[i_step-1], 1, 0)
        strain[i_step] = elasto_alum.strain
        stress[i_step] = elasto_alum.stress
        back_stress[i_step] = elasto_alum.back_stress
        yield_stress[i_step] = elasto_alum.yield_stress
    plt.figure()
    plt.plot(strain*10**3, stress, "-k", label="Stress")
    plt.plot(strain*10**3, back_stress, ":b", label="Back Stress")
    plt.plot(strain*10**3, yield_stress, "--r", label="Yield Stress")
