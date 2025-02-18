import numpy as np
from dataclasses import dataclass

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
    tol: float
    max_iter: int

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
        if len(set_strain) < 1: raise Exception("Array of strain increments should have at least one element")
        if hard_iso < 0 or hard_iso > 1: raise Exception("Isotropic hardening parameter must be between 0 (off) and 1 (on)")
        if hard_kin < 0 or hard_kin > 1: raise Exception("Kinematic hardening parameter must be between 0 (off) and 1 (on)")
        
        for strain_incr in set_strain:
            self.strain = self.strain + strain_incr

            # Update stress, back stress, and yield surface (Predictor-Corrector in 1 step using min())
            elastic_stress = self.stress + self.mat.E*strain_incr
            phi = min(0, abs(elastic_stress - self.back_stress) - self.yield_stress)
            plastic_strain = phi / (self.mat.E + self.mat.H)
            signed_strain = plastic_strain * np.sign(elastic_stress - self.back_stress)
            self.stress = elastic_stress - self.mat.E * signed_strain

            # Isotropic Part - Von-Mises Yield Criterion -> d(yield_stress) = 3/2 * h * d(plastic_strain)
            self.yield_stress = self.yield_stress + hard_iso * self.mat.H * plastic_strain
            # Kinematic Part - Prager's Rule -> d(back_stress) = c * d(plastic_strain)
            self.back_stress = self.back_stress + hard_kin * self.mat.H * signed_strain

            ## Compare with lecture notes - confusing update steps
        


        