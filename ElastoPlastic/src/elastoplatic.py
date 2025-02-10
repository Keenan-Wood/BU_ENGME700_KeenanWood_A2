import numpy as np
from dataclasses import dataclass

@dataclass
class ElastoPlastic:
    tol: float
    max_iter: int

    def __init__(self, fun, vars_indep, start_pt, J = None, tol: float = 10**-12, max_iter: int = 10**3):
        self.vars_indep = vars_indep
        self.dim = len(self.vars_indep)
        self.fun = fun
        self.pt = start_pt
        self.J = J
        self.tol = tol
        self.max_iter = max_iter
        self.num_iter = 0

    # Validate Instance Inputs
    @property
    def vars_indep(self):
        return self._vars_indep
    
    @vars_indep.setter
    def vars_indep(self, input_vars_indep):
        if not all(isinstance(var, sp.core.symbol.Symbol) for var in input_vars_indep): raise Exception("Function is not a symbolic expression")
        if len(input_vars_indep) != len(set(input_vars_indep)): raise Exception("Duplicate variables provided")
        self._vars_indep = input_vars_indep