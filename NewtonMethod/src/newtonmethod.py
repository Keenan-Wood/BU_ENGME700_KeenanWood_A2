import numpy as np
import sympy as sp
from dataclasses import dataclass

# Class of an iterable object "NewtonMethod" representing an estimate of the zero of a provided symbolic function
# Properties: Function (fun), Jacobian (J), independent variable vector (vars_indep)
#             x-coordinate (pt), value at x-coordinate (ptVal)
# Solver parameters: convergence criteria (tol) and maximum number of iterations (max_iter)
@dataclass
class NewtonMethod:
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

    @property
    def fun(self):
        return self._fun
    
    @fun.setter
    def fun(self, input_fun):
        if all(isinstance(element, sp.Expr) for element in input_fun):
            if self.dim == 0: raise Exception("Symbolic function has no independent variables")
            fun_vars = [list(fun_component.free_symbols) for fun_component in input_fun]
            if not all(fun_var in self.vars_indep for fun_var in fun_vars[0]): raise Exception("Function contains unrecognized variable")
        # ADD AND TEST >> else: raise Exception("Function is not a symbolic expression")
        self._fun = input_fun

    @property
    def J(self):
        return self._J
    
    @J.setter
    def J(self, input_J):
        if input_J is None:
            self.J = self.fun.jacobian(self.vars_indep)
        else:
            if not all(isinstance(element, sp.Expr) for element in input_J): raise Exception("Jacobian is not a symbolic expression")
            # ADD AND TEST >> code to verify Jacobian size/shape and variables
            self._J = input_J

    @property
    def tol(self):
        return self._tol

    @tol.setter
    def tol(self, input_tol):
        if input_tol <= 0: raise Exception("Tolerance must be greater than zero")
        self._tol = input_tol

    @property
    def max_iter(self):
        return self._max_iter

    @max_iter.setter
    def max_iter(self, input_max_iter):
        if input_max_iter <= 0: raise Exception("Max number of iterations must be greater than zero")
        self._max_iter = input_max_iter

    
    # Define iteration method for the "NewtonMethod" object so that
    # each iteration produces the next estimate of the zero (and the corresponding function evaluation) 
    #   - Stop iteration on error (function evaluation issue or singular Jacobian), 
    #     or if a solver termination criteria is met (convergence or maximum number of iterations)  
    def __iter__(self):
        return self
    
    def __next__(self):
        self.num_iter = self.num_iter + 1
        try:
            sub_vars = [(old_var, new_val) for (old_var, new_val) in zip(self.vars_indep, self.pt)]
            self.ptVal = np.array(self.fun.subs(sub_vars)).astype(np.float64)
            self.ptVal = self.ptVal[:,0]
        except:
            print("Function could not be evaluated at x=", self.pt)
            raise StopIteration
        self.JVal = np.array(self.J.subs(sub_vars)).astype(np.float64)
        if np.linalg.det(self.JVal) < 10**-8:
            print("Singular Jacobian reached at x=", self.pt)
            raise StopIteration
        if self.dim == 1:
            self.pt = self.pt - self.ptVal[0] / self.JVal[0]
        else:
            self.pt = self.pt - np.matmul(np.linalg.inv(self.JVal), self.ptVal)
        if np.linalg.norm(self.ptVal) <= self.tol or self.num_iter > self.max_iter:
            raise StopIteration
        return self