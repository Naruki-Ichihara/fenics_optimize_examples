---
sidebar_position: 1
---

# A First example

## Theory
Let's start with your first optimization example. If you are the first of FEniCS computing, please see the [tutorial of the FEniCS](https://fenicsproject.org/documentation/).
This example is based on the dolfin-adjoint's documentation. Please also see [this page](http://www.dolfin-adjoint.org/en/latest/documentation/poisson-topology/poisson-topology.html).

Here, we are concerned about minimizing the following compliance:

$$
\int_\Omega fT\,dx
$$

that subject to the Poisson equation with following boundary constraints

$$
-\nabla\cdot(k(a)\nabla T) = f \,\, \mathrm{in}\,\Omega
$$
$$
T = 0\,\,\mathrm{on}\,\Gamma_D
$$
$$
k(a)\nabla T = 0\,\,\mathrm{on}\,\Gamma_N
$$

where, $\Omega$ is the design domain, $T$ is the scalar field, $a$ is the control, $f$ is the force term, and $k(a)$ is the Solid isotropic material with penalization (SIMP) parameterization. Physically, this optimization problem is to find material distribution $a(x)$ that minimizes the temperature $T$ under the material limits (see for [reference](https://link.springer.com/article/10.1007/s00158-005-0584-3)). The total material amount bounds the control $a$:

$$
0\leq a(x) \leq 1
$$
$$
\int_\Omega a(x)\,dx \leq m
$$

here, $m$ is the target material amount. In this example, the Helmholtz regularization is used to avoid numerical instability (see for [article](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.3072)).

## Implimentation 

First, import the `fenics_optimize` with `dolfin`, `dolfin_adjoint` and `numpy` packages.
```python
from dolfin import *
from dolfin_adjoint import *
import numpy as np
from fenics_optimize as op
```
Note that, `dolfin_adjoint` must be called after the `dolfin`.
Next, we define the parameters

```python
m = 0.3    # Target rate of the material amount
p = 5      # Penalty parameter
eps = 1e-3 # Material lower bound
R = 0.1    # Helmholtz filter radius
n = 250    # Resolution
```

Next we define the SIMP function

```python
def k(a):
    return eps + (1 - eps) * a ** p
```

Then the unit square design space is defined by

```python
mesh = UnitSquareMesh(n, n)
X = FunctionSpace(mesh, 'CG', 1)
```

where, design space is discretized by the first order Continuous Galerkin space.

Next, we define the sub-domain for the Dirichlet boundary condition:

```python
class Left(SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/n + eps
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary

bc = DirichletBC(X, Constant(0.0), Left())
```

```
image is here
```

Next, we define the function space of the model:

```python
U = FunctionSpace(mesh, 'CG', 1)
T = Function(U, name='Temperture')
t = TrialFunction(U)
dt = TestFunction(U)
```

Now we can define the `forward` function with the weak formation of the model that return the cost:

```python
def forward(x):
    a = op.helmholtzFilter(x, X, R=r)
    a = inner(grad(t), k(a)*grad(dt))*dx
    L = f*dt*dx
    A, b = assemble_system(a, L, bcs)
    T = op.AMGsolver(A, b).solve(T, U, False)
    J = assemble(inner(grad(T), k(a)*grad(T))*dx)
    return J
```

The `op.AMGsolver` provides efficient and robust FE solvers for multiphysics. 

In the case of the topology optimization, the Jacobian is needed to optimize because the number of controls is enormous.
Generally, the Jacobian is calculated by the analytical approach. However, the Jacobian will be calculated automatically by using the decorator, `with_derivative`.

To ready the automatic derivative, decorate the `forward` function by the `with_derivative([X])`.

```python
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
```

Next we define the mass volume amount function with derivative:

```python
@op.with_derivative([X])
def constraint(xs):
    a_bulk = project(Constant(1.0), X)
    a_0 = assemble(a_bulk*dx)
    a_f = assemble(xs[0]*dx)
    rel = a_f/a_0
    return rel - m
```

Now we ready for the optimization: Finally, the initial value is defined

```python
problemSize = Function(X).vector().size()
x0 = np.ones(problemSize) * m
```

and start optimization using `MMAoptimize`

```python
op.MMAoptimize(problemSize, x0, forward, [constraint], [0.0], maxeval=100, bounds=[0, 1], rel=1e-20)
```

```
Here is result image (TODO).
```
## References
1. A. Gersborg-Hansen et. al., [Topology optimization of heat conduction problems using the finite volume method](https://link.springer.com/article/10.1007/s00158-005-0584-3), *Structural and Multidisciplinary Optimization*, **31** (2006), 251-259
1. B. S. Lazarov et. al., [Filters in topology optimization based on Helmholtz-type differential equations](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.3072), *International Journal for Numerical Methods in Engineering*, **28** (2010), 765-781

## Complete code
```python
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

op.MMAoptimize(problemSize, x0, forward, [constraint], maxeval=100, bounds=[0, 1], rel=1e-20)
```


