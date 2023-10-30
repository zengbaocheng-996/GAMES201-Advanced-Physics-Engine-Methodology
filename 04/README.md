# Finite Elements and Topology Optimization

### Finite element method

Finite element method (FEM) belongs to the family of Galerkin methods. In FEM, continuous PDEs are converted to discrete linear systems. Understanding FEM is important for many other discretization methods, including the Material Point Method (MPM) later in this course.

###### Typical steps

1. Convert strong-form PDEs to weak forms, using a test function w
2. Integrate by parts to redistribute gradient operators
3. Use the divergence theorem to simplify equations and enforce Neumann Boundary conditions (BCs)
4. Discretization (build the stiffness matrix and right-hand side)
5. Solve the (discrete) linear system

###### 2D Poisson's equation

pressure projection in fluid simulations
$$
\nabla\cdot\nabla u=0
$$
Dirichlet boundary 第一类边界条件
$$
u(x)=f(x),x\in \partial\Omega
$$
Neumann boundary 第二类边界条件
$$
\nabla u(x)\cdot n = g(x),x\in \partial\Omega
$$

### Discretizing Poisson's equation

###### Weak formulation

Arbitrary 2D test function w(x)
$$
\nabla\cdot\nabla u=0\Longleftrightarrow \forall w,\iint_\Omega w(\nabla\cdot\nabla u)dA=0
$$

###### Getting rid of second-order derivatives

$$
\begin{align}
\nabla\cdot\nabla u&=0\\
\nabla w\cdot\nabla u+w\nabla\cdot\nabla u&= \nabla\cdot(w\nabla u)\\
\nabla w\cdot\nabla u&=\nabla\cdot(w\nabla u)\

\end{align}
$$

$$
\nabla\cdot\nabla=0\Longleftrightarrow\forall w,\iint_{\Omega}\nabla w\cdot\nabla udA=\iint_{\Omega}\nabla\cdot(w\nabla u)dA
$$

$$
\iint_{\Omega}\nabla w\cdot\nabla udA=\oint_{\partial\Omega}w\nabla u\cdot dn
$$

###### Applying discretization Basis function

$$
\begin{align}
u(x)&=\sum_{j}u_j\phi_j(x)\\
\iint_{\Omega}\nabla w\cdot\nabla udA&=\oint_{\partial\Omega}w\nabla u\cdot dn\\
\forall w,\iint_{\Omega}\nabla w\cdot\nabla(\sum_{j}u_j\phi_j)dA&=\oint_{\partial\Omega}w\nabla u\cdot dn\\
\forall i,\iint_{\Omega}\nabla \phi_i\cdot\nabla(\sum_{j}u_j\phi_j)dA&=\oint_{\partial\Omega}\phi_i\nabla u\cdot dn\\
\forall i,\sum_{j}(\iint_{\Omega}\nabla \phi_i\cdot\nabla \phi_jdA)u_j&=\oint_{\partial\Omega}\phi_i\nabla u\cdot dn\\
Ku&=f\\
K_{ij}&=\iint_{\Omega}\nabla \phi_i\cdot\nabla \phi_jdA
\end{align}
$$

K: stiffness matrix

u: degree of freedoms / solution vector

f: load vector

###### Boundary Conditions

1. Dirichlet boundary conditions
   $$
   u(x)=f(x),x\in\partial\Omega
   $$

2. Neumann boundary conditions
   $$
   \nabla u(x)\cdot n=g(x),x\in\partial\Omega
   $$

### Discretizing linear elasticity

###### Linear elasticity FEM

$$
\frac{Dv}{Dt}=\frac{1}{\rho}\nabla\cdot\sigma+g
$$

v: velocity

rho: density

sigma: Cauchy stress tensor (symmetric 2/3D "matrix")

g: body force (e.g., gravity)

Quasistatic state (v = 0), constant density, no gravity
$$
\nabla\cdot\sigma=0
$$
Degree of freedom: displacement u sigma = sigma(u)

Infinitesimal deformation: Lagrangian/Eulerian classification does not make sense

###### Index notation

$$
\frac{Dv}{Dt}=\frac{1}{\rho}\nabla\cdot\sigma+g\Longleftrightarrow\frac{1}{\rho}\sum_{\beta}\sigma_{\alpha\beta,\beta}+g_\alpha
$$

###### Discretize Cauchy momentum equation using FEM

$$
\sum_\beta\sigma_{\alpha\beta,\beta}w_\alpha=0\\
\sum_\beta\sigma_{\alpha\beta,\beta}w_\alpha+\sum_\beta\sigma_{\alpha\beta}w_{\alpha,\beta}=\sum_\beta(\sigma_{\alpha\beta}w_\alpha)_{,\beta}\Rightarrow\sum_{\beta}\sigma_{\alpha\beta}w_{\alpha,\beta}=\sum_\beta(\sigma_{\alpha\beta}w_\alpha)_{,\beta}
$$

###### Divergence theorem

$$
\forall\alpha\forall w,\iint_{\Omega}\sum_{\beta}\sigma_{\alpha\beta}w_{\alpha,\beta}dA=\oint_{\partial\Omega}\sum_\beta(\sigma_{\alpha\beta}w_\alpha)dn_\beta
$$

###### Discretization

$$
\begin{align}
\forall\alpha\forall w,\iint_{\Omega}\sum_{\beta}\sigma_{\alpha\beta}w_{\alpha,\beta}dA&=\oint_{\partial\Omega}\sum_\beta(\sigma_{\alpha\beta}w_\alpha)dn_\beta\\
w_\alpha(x)=\sum_i w_{i\alpha}\phi_{i\alpha}(x)&,u_\alpha(x)=\sum_j u_{j\alpha}\phi_{j\alpha}(x)\\
\forall\alpha\forall i,\iint_{\Omega}\sum_{\beta}[\sigma(u(x))]_{\alpha\beta}\phi_{i\alpha}(x)dA&=\oint_{\partial\Omega}\sum_\beta(\sigma_{\alpha\beta}\phi_{i\alpha})dn_\beta
\end{align}
$$

###### Relating sigma to u

Strain tensor
$$
e=\frac{1}{2}(\nabla u+(\nabla u)^T)
$$
Cauchy stress tensor
$$
\sigma=\lambda tr(e)I+2\mu e
$$
Cauchy stress tensor (index notation)
$$
\begin{align}
e_{\alpha\beta}&=\frac{1}{2}(u_{\alpha,\beta}+u_{\beta,\alpha})\\
\sigma_{\alpha\beta}&=\lambda\delta_{\alpha\beta}\sum_{\alpha}e_{\alpha\alpha}+2\mu e_{\alpha\beta}\\
\delta_{\alpha\beta}&=\begin{cases}
1 \ \ \ {\mbox{if}\ \ \ \alpha = \beta}\\0 \ \ \ {\mbox{if}\ \ \ \alpha \neq \beta}
\end{cases}
\end{align}
$$

###### Building the linear system

$$
\forall\alpha\forall i,\iint_{\Omega}\sum_{\beta}[\sigma(u(x))]_{\alpha\beta}\phi_{i\alpha}(x)dA=\oint_{\partial\Omega}\sum_\beta(\sigma_{\alpha\beta}\phi_{i\alpha})dn_\beta
$$



### Topology optimization

Keywords: Solid Isotropic Material with Penalization (SIMP), Optimility Criterion (OC)

The minimal compliance topology optimization problem can be formulated as
$$
\begin{align}
min\ L(\rho)&=u^TK(\rho)u\\
s.t. \ K(\rho)u&=f\\
\sum_e\rho_e&\leq cV,\\
\rho_e&\in[\rho_{min},1]
\end{align}
$$
L: measure of deformation energy, or the loss function

c: volume fraction (e.g., 0.3)

\rho_e: material occupancy (0 = empty, 1 = filled) of cell e

V: total volume
