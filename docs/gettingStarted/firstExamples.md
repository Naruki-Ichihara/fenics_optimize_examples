---
sidebar_position: 1
---

# A First example
Let's start your first optimization example. If you are first of FEniCS computing, please see the tutorial of the FEniCS.
This example is based on the dolfin-adjoint's documentation. Please also see this page.

Here, we concern minimizing following compliance:

$$
\int_\Omega fT\,dx
$$

that subject to the Poisson equation with boundary constraints

$$
-\nabla\cdot(k(a)\nabla T) = f \,\, \mathrm{in}\,\Omega
$$
$$
T = 0\,\,\mathrm{on}\,\Gamma_D
$$
$$
k(a)\nabla T = 0\,\,\mathrm{on}\,\Gamma_N
$$

where, $\Omega$ is the design domain, $T$ is the schalar field, $a$ is the control, $f$ is the force term, and $k(a)$ is the Solid isotropic material with penalisation (SIMP) parameterisation. Physically, this optimization problem is to find material distribution $a(x)$ that minimizes the total of the temperture $T$ under the material limits. The control $a$ is bounded by the total material amount:

$$
0\leq a(x) \leq 1
$$
$$
\int_\Omega a(x)\,dx \leq m
$$

here, $m$ is the target material amount. In this example, the Helmholtz regularization is used to avoid the numerical instability (see for example).

## Implimentation 

First, import the `fenics_optimize` with `dolfin` and `dolfin_adjoint` packages.
```python
from dolfin import *
from dolfin_adjoint import *
from fenics_optimize as op
```
Note that, `dolfin_adjoint` must be called after the `dolfin`.
Next we define the parameters

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

where, design space is descritezed by the first order Continuous Galerkin space.

Next we define the sub-domain for the Dirichlet boundaty condition:

```python
class Left(SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/n + eps
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary

bc = DirichletBC(X, Constant(0.0), Left())
```

this is the following boundary setting.

```
image is here
```

Next we define the function space of the mathematical model:

```python
U = FunctionSpace(mesh, 'CG', 1)
T = Function(U, name='Temperture')
t = TrialFunction(U)
dt = TestFunction(U)
```

Now we can define the `forward` function that return the cost function:

```python
def forward(x):
    a = op.helmholtzFilter(x, X, R=r)
    a.rename('label', 'control')
    a = inner(grad(t), k(a)*grad(dt))*dx
    L = f*dt*dx
    A, b = assemble_system(a, L, bcs)
    T = op.AMGsolver(A, b).solve(T, U, False)
    J = assemble(inner(grad(T), k(a)*grad(T))*dx)
    return J
```

In the case of the topology optimization, the Jacobian is needed to optimize because the number of controls is huge.
Generaly, the Jacobian is calculated by the analytical approach. However, by using the `fenics_optimize`, the jacobian will be calucalated automatically.

To ready the automatic derivative, decorate the `forward` function by the `with_derivative([X])`.

```python
@op.with_derivative([X])
def forward(x):
    a = op.helmholtzFilter(x, X, R=r)
    a.rename('label', 'control')
    a = inner(grad(t), k(a)*grad(dt))*dx
    L = f*dt*dx
    A, b = assemble_system(a, L, bcs)
    T = op.AMGsolver(A, b).solve(T, U, False)
    J = assemble(inner(grad(T), k(a)*grad(T))*dx)
    return J
```

Next we define the mass volume amount function by the `fenics_optimize` scope:

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
op.MMAoptimize(problemSize, x0, forward, constraint, maxeval=1000, bounds=[0, 1], rel=1e-20)
```


