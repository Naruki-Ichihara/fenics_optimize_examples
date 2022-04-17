from dolfin import *
from dolfin_adjoint import *
import numpy as np
from fenics_optimize as op

m = 0.3    # Target rate of the material amount
p = 5      # Penalty parameter
eps = 1e-3 # Material lower bound
R = 0.1    # Helmholtz filter radius
n = 250    # Resolution

def k(a):
    return eps + (1 - eps) * a ** p

mesh = UnitSquareMesh(n, n)
X = FunctionSpace(mesh, 'CG', 1)

class Left(SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/n + eps
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary

bc = DirichletBC(X, Constant(0.0), Left())

U = FunctionSpace(mesh, 'CG', 1)
T = Function(U, name='Temperture')
t = TrialFunction(U)
dt = TestFunction(U)

@op.with_derivative([X])
def forward(x):
    a = op.helmholtzFilter(x[0], X, R=r)
    a.rename('label', 'control')
    a = inner(grad(t), k(a)*grad(dt))*dx
    L = f*dt*dx
    A, b = assemble_system(a, L, bcs)
    T = op.AMGsolver(A, b).solve(T, U, False)
    J = assemble(inner(grad(T), k(a)*grad(T))*dx)
    return J

@op.with_derivative([X])
def constraint(xs):
    a_bulk = project(Constant(1.0), X)
    a_0 = assemble(a_bulk*dx)
    a_f = assemble(xs[0]*dx)
    rel = a_f/a_0
    return rel - m

problemSize = Function(X).vector().size()
x0 = np.ones(problemSize) * m

op.MMAoptimize(problemSize, x0, forward, [constraint], constraints_targets=[0.0], maxeval=100, bounds=[0, 1], rel=1e-20)