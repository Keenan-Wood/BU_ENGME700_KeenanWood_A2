import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt

# A "material" class is useful to store basic material properties.
# Properties: Elastic Modulus "E", Plastic Modulus "H", and un-deformed yield stress "yield_stress"
# This eliminates the need for specifying these properties many times for multiple objects of the same material
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

# An "ElastoPlastic" class represents a deformable material, composed of the assigned "material".
# Properties: material, current yield stress, current stress, current strain, and current back stress

@dataclass
class ElastoPlastic:

    def __init__(self, mat: material, strain: float, stress: float, back_stress: float = 0):
        self.mat = mat
        self.yield_stress = self.mat.yield_stress
        self.strain = strain
        self.stress = stress
        self.back_stress = back_stress

    # The stretch() method deforms the object according to the provided strain increments
    # behavior in the plastic regime is determined by isotropic and kinematic hardening parameters
    # An (isotropic, kinematic) parameter pair of (1,0) gives standard isotropic hardening behavior, and vice versa,
    # while intermediate values representing mixed hardening are permitted
    def stretch(self, set_strain: np.array, hard_iso: float, hard_kin: float):
        # Validate inputs
        if set_strain.size < 1: raise Exception("Array of strain increments should have at least one element")
        if hard_iso < 0 or hard_iso > 1: raise Exception("Isotropic hardening parameter must be between 0 (off) and 1 (on)")
        if hard_kin < 0 or hard_kin > 1: raise Exception("Kinematic hardening parameter must be between 0 (off) and 1 (on)")
        
        # If set_strain has only one value, it needs to be made into a list to be iterable
        if set_strain.size == 1: set_strain = [set_strain]

        for strain_incr in set_strain:
            self.strain = self.strain + strain_incr

            # Update stress, back stress, and yield surface (Predictor-Corrector without branching using max())
            # The elastic_stress is the stress if the deformation from the strain increment is entirely elastic
            elastic_stress = self.stress + self.mat.E*strain_incr
            phi = max(0, abs(elastic_stress - self.back_stress) - self.yield_stress)
            plastic_strain = phi / (self.mat.E + self.mat.H)
            signed_strain = plastic_strain * np.sign(elastic_stress - self.back_stress)
            self.stress = elastic_stress - self.mat.E * signed_strain

            # Isotropic Hardening - expand the yield surface
            self.yield_stress = self.yield_stress + hard_iso * self.mat.H * plastic_strain
            # Kinematic update - shift the yield surface's center (the back stress)
            self.back_stress = self.back_stress + hard_kin * self.mat.H * signed_strain