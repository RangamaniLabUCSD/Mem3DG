from Problem import problem 
from scipy.optimize import minimize as mini

class optimizer(object):
    
    def __init__(self,penalty,ex_force,mesh,nu_grad_size = 1e-8):
        self.prob = problem(penalty,ex_force,mesh,nu_grad_size) 
        
    def minimize_newton(self):
        
        x_new = mini(self.prob.obj_func, \
        self.prob.x0, args=(), method='Newton-CG', \
        jac=self.prob.gradient, hess=self.prob.hessian, hessp=None, tol=None, callback=None, \
        options={'xtol': 1e-05, 'eps': 1.4901161193847656e-08, 'maxiter': None, 'disp': True, 'return_all': None});
           
    
        # x_new = mini(self.prob.obj_func, \
        # self.prob.x0, args=(), method='BFGS', \
        #  jac=self.prob.gradient, tol=None, callback=None, \
        #  options={'gtol': 1e-05, 'norm': 2, 'eps': 1.4901161193847656e-08, 'maxiter': None, 'disp': True, 'return_all': True})    
            
        return x_new.x 