import numpy as np
from dataclasses import dataclass

@dataclass
class material:
    E: float
    H: float
    Y_0: float

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

    def stretch(self, set_strain: np.array, hard_iso: float, hard_kin):
        # Validate inputs
        if len(set_strain) < 1: raise Exception("Array of strain increments should have at least one element")
        if abs(hard_iso-.5) > .5: raise Exception("Isotropic hardening parameter must be between 0 (off) and 1 (on)")
        if abs(hard_kin-.5) > .5: raise Exception("Kinematic hardening parameter must be between 0 (off) and 1 (on)")
        
        for strain_incr in set_strain:
            # Update stress and strain (Predictor-Corrector in 1 step)
            self.yield_stress = self.yield_stress + hard_iso*self.mat.H*self.strain
            self.strain = self.strain + strain_incr

            elastic_stress = self.stress + self.mat.E*strain_incr - self.back_stress
            phi = min(0, abs(elastic_stress) - self.yield_stress)
            plastic_strain = phi / (self.mat.E + self.mat.H)
            directed_strain = plastic_strain * np.sign(elastic_stress)
            self.stress = elastic_stress - self.mat.E * directed_strain
            self.back_stress = self.back_stress + hard_kin * self.mat.H * directed_strain

            ## Compare with lecture notes - confusing update steps
        


        