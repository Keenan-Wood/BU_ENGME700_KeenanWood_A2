import numpy as np
import sympy as sp
from dataclasses import dataclass

# Class storing a vector-valued vector function, its Jacobian, and two points;
#   as well as the scalar-valued function defined as the square of the norm of the vector function, and its gradient;
#   as well as the parameterized function defined as the scalar-valued function evaluated on the line between the two points,
#   and the gradient projected onto the line between the two points (ie. the derivative of the parameterized function)
#     Purpose: To collapse a bracketed N-dimensional root finding problem onto one dimension
@dataclass
class ParametricFunction:

    def __init__(self, fun, vars_indep: list, pt_a: np.array, pt_b: np.array, J = None):
        self.vars_indep = vars_indep
        self.dim = len(vars_indep)
        self.fun = fun
        self.pt_a = pt_a
        self.pt_b = pt_b
        self.J = J
        self.to_scalar()
        self.parameterize()
    
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
        else: raise Exception("Function is not a symbolic expression")
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
            # Add code to verify Jacobian size/shape and variables
            self._J = input_J

    @property
    def pt_a(self):
        return self._pt_a
    
    @pt_a.setter
    def pt_a(self, input_pt_a):
        if input_pt_a.size != self.dim: raise Exception("not enough coordinates given for point a")
        self._pt_a = input_pt_a

    @property
    def pt_b(self):
        return self._pt_b
    
    @pt_b.setter
    def pt_b(self, input_pt_b):
        if input_pt_b.size != self.dim: raise Exception("not enough coordinates given for point b")
        if all(self.pt_a == input_pt_b): raise Exception("point a and point b must be different")
        self._pt_b = input_pt_b

    # Function to create a scalar function and its gradient from a given vector function and its Jacobian
    # The scalar function created is equal to the magnitude of the vector function squared (ie. ||fun||^2)
    def to_scalar(self):
        # if (isinstance(self.fun, sp.Basic)):
        self.fun_scalar = self.fun.dot(self.fun)
        self.grad = 2 * self.J.T * self.fun
        # else:
        #     self.fun_scalar = lambda x: np.dot(self.fun(x), self.fun(x))
        #     self.grad = lambda x: 2 * np.matmul(self.J(x).transpose(), self.fun(x))
        return self
    
    # Function to parameterize the scalar function and its derivative (the projected gradient) along the line from pt_a to pt_b
    # The scalar value x ranges from 0 to 1 and defines the position in space as (pt_a + x*(pt_b-pt_a))
    def parameterize(self):
        vector = self.pt_b - self.pt_a
        # if (isinstance(self.fun, sp.Basic)):
        x = sp.symbols('x')
        param_vars = [(old_var, new_val) for (old_var, new_val) in zip(self.vars_indep, self.pt_a + x*vector)]
        sym_fun_param = self.fun_scalar.subs(param_vars)
        sym_grad_param = self.grad.subs(param_vars)
        self.fun_param = sp.lambdify([x], sym_fun_param, "numpy")
        self.grad_param = sp.lambdify([x], sym_grad_param, "numpy")
        # else:
        #     self.fun_param = lambda x: self.fun_scalar(self.pt_a + x*vector)
        #     self.grad_param = lambda x: np.dot(vector, self.grad(self.pt_a + x*vector))
        return self