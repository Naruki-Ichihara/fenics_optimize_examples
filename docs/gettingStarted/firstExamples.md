---
sidebar_position: 1
---

# A First example

## Theory
Start with your first optimization example. If you are the first of FEniCS computing, please see the [tutorial of the FEniCS](https://fenicsproject.org/documentation/).
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

First, import the `optfx` with `numpy` package.
```python
import numpy as np
import optfx as opt
```
Next, we define the parameters

```python
opt.parameters["form_compiler"]["optimize"] = True
opt.parameters["form_compiler"]["cpp_optimize"] = True
opt.parameters['form_compiler']['quadrature_degree'] = 5

comm = opt.MPI.comm_world
m = 0.30   # Target rate of the material amount
p = 5      # Penalty parameter
eps = 1e-3 # Material lower bound
R = 0.01   # Helmholtz filter radius
n = 256    # Resolution
f = opt.interpolate(opt.Constant(1e-2), X)
```
Notably, `opt.MPI.comm_world` enables the parallel computing with MPI.

Next we define the SIMP function

```python
def k(a):
    return eps + (1 - eps) * a ** p
```

Then the unit square design space is defined by

```python
mesh = opt.UnitSquareMesh(comm, n, n)
X = opt.FunctionSpace(mesh, 'CG', 1)
```

where, design space is discretized by the first order Continuous Galerkin space.

Next, we define the sub-domain for the Dirichlet boundary condition:

```python
class Left(opt.SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/250 + 1e-5
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary
```

Next, we define the function space of the model:

```python
U = opt.FunctionSpace(mesh, 'CG', 1)
t = opt.TrialFunction(U)
dt = opt.TestFunction(U)
bc = opt.DirichletBC(U, opt.Constant(0.0), Left())
```

Now we can define the `Problem` function with the weak formation of the model that return the cost:

```python
class PoissonProblem():
    def problem(self, controls):
        rho = controls[0]
        rho = opt.helmholtzFilter(rho, X, R)
        a = opt.inner(opt.grad(t), k(rho)*opt.grad(dt))*opt.dx
        L = opt.inner(f, dt)*opt.dx
        Th = opt.Function(U, name='Temperture')
        opt.solve(a==L, Th, bc)
        J = opt.assemble(opt.inner(opt.grad(Th), k(rho)*opt.grad(Th))*opt.dx)
        rho_bulk = opt.project(opt.Constant(1.0), X)
        rho_0 = opt.assemble(rho_bulk*opt.dx)
        rho_total = opt.assemble(controls[0]*opt.dx)
        rel = rho_total/rho_0
        self.volumeFraction = rel
        return J
    def constraint_volume(self):
        return self.volumeFraction - m
```

In the case of the topology optimization, the sensitivity is needed to optimize because the number of controls is enormous.
Generally, the sensitivity is calculated by the analytical approach. However, the sensitivity will be calculated automatically by using the optfx.

To ready the automatic derivative, `opt.Module` can be used like as:

```python
class PoissonProblem(opt.Module):
    def problem(self, controls):
        rho = controls[0]
        rho = opt.helmholtzFilter(rho, X, R)
        a = opt.inner(opt.grad(t), k(rho)*opt.grad(dt))*opt.dx
        L = opt.inner(f, dt)*opt.dx
        Th = opt.Function(U, name='Temperture')
        opt.solve(a==L, Th, bc)
        J = opt.assemble(opt.inner(opt.grad(Th), k(rho)*opt.grad(Th))*opt.dx)
        rho_bulk = opt.project(opt.Constant(1.0), X)
        rho_0 = opt.assemble(rho_bulk*opt.dx)
        rho_total = opt.assemble(controls[0]*opt.dx)
        rel = rho_total/rho_0
        self.volumeFraction = rel
        return J
    def constraint_volume(self):
        return self.volumeFraction - m
```

Next we define the mass volume amount function with derivative:
Now we ready for the optimization: Finally, the initial value is defined

```python
class Initial_density(opt.UserExpression):
    def eval(self, value, x):
        value[0] = m
x0 = opt.Function(X)
x0.interpolate(Initial_density())
```

and start optimization using `optimize`:

```python
N = opt.Function(X).vector().size()
min_bounds = np.zeros(N)
max_bounds = np.ones(N)

setting = {'set_lower_bounds': min_bounds,
           'set_upper_bounds': max_bounds,
           'set_maxeval': 1000
          }
params = {'verbosity': 1}

problem = PoissonProblem()
solution = opt.optimize(problem, [x0], [0], setting, params)
```

## References
1. A. Gersborg-Hansen et. al., [Topology optimization of heat conduction problems using the finite volume method](https://link.springer.com/article/10.1007/s00158-005-0584-3), *Structural and Multidisciplinary Optimization*, **31** (2006), 251-259
1. B. S. Lazarov et. al., [Filters in topology optimization based on Helmholtz-type differential equations](https://onlinelibrary.wiley.com/doi/full/10.1002/nme.3072), *International Journal for Numerical Methods in Engineering*, **28** (2010), 765-781

## Complete code
```python
import numpy as np
import optfx as opt

opt.parameters["form_compiler"]["optimize"] = True
opt.parameters["form_compiler"]["cpp_optimize"] = True
opt.parameters['form_compiler']['quadrature_degree'] = 5

comm = opt.MPI.comm_world
m = 0.30   # Target rate of the material amount
p = 5      # Penalty parameter
eps = 1e-3 # Material lower bound
R = 0.01   # Helmholtz filter radius
n = 256    # Resolution

def k(a):
    return eps + (1 - eps) * a ** p

mesh = opt.UnitSquareMesh(comm, n, n)
X = opt.FunctionSpace(mesh, 'CG', 1)
f = opt.interpolate(opt.Constant(1e-2), X)

class Initial_density(opt.UserExpression):
    def eval(self, value, x):
        value[0] = m

class Left(opt.SubDomain):
    def inside(self, x, on_boundary):
        gamma = 1/250 + 1e-5
        return x[0] == 0.0 and 0.5 - gamma < x[1] < 0.5 + gamma and on_boundary

U = opt.FunctionSpace(mesh, 'CG', 1)
t = opt.TrialFunction(U)
dt = opt.TestFunction(U)
bc = opt.DirichletBC(U, opt.Constant(0.0), Left())

class PoissonProblem(opt.Module):
    def problem(self, controls):
        rho = controls[0]
        rho = opt.helmholtzFilter(rho, X, R)
        a = opt.inner(opt.grad(t), k(rho)*opt.grad(dt))*opt.dx
        L = opt.inner(f, dt)*opt.dx
        Th = opt.Function(U, name='Temperture')
        opt.solve(a==L, Th, bc)
        J = opt.assemble(opt.inner(opt.grad(Th), k(rho)*opt.grad(Th))*opt.dx)
        rho_bulk = opt.project(opt.Constant(1.0), X)
        rho_0 = opt.assemble(rho_bulk*opt.dx)
        rho_total = opt.assemble(controls[0]*opt.dx)
        rel = rho_total/rho_0
        self.volumeFraction = rel
        return J
    def constraint_volume(self):
        return self.volumeFraction - m

x0 = opt.Function(X)
x0.interpolate(Initial_density())
N = opt.Function(X).vector().size()
min_bounds = np.zeros(N)
max_bounds = np.ones(N)

setting = {'set_lower_bounds': min_bounds,
           'set_upper_bounds': max_bounds,
           'set_maxeval': 1000
          }
params = {'verbosity': 1}

problem = PoissonProblem()
solution = opt.optimize(problem, [x0], [0], setting, params)
```


