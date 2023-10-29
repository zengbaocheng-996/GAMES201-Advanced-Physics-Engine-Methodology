# Eulerian fluid simulation

For each time step,

###### Advection: "move" the fluid field. Solve for u* using u^t

Key: Use a advection scheme with low numerical viscosity (e.g., MacCormack/BFECC/Particle advection)

$$
\frac{Du}{Dt}=0,\ \ \ \frac{D\alpha}{Dt}=0
$$

###### Projection: make velocity field u^{t+1} divergence-free based on u*

Key: Use a fast linear solver (e.g., MGPCG)

$$
\frac{\partial u}{\partial t}=-\frac{1}{\rho}\nabla p\ \ \ \ \ \ s.t.\ \ \ \ \ \ \nabla\cdot u^{t+1}=0
$$

###### Other

SIGGRAPH 2015 Combining advection with reflection: IVOCK

SIGGRAPH 2018 Advection-Relfection solver

###### Possible extensions

1. Going to 3D
2. Accurate boundary conditions and fluid-solid coupling
3. Two phase fluid simulation
4. Handling free surfaces (level sets)
5. Vortex methods

### Material Derivatives

$$
\frac{D}{Dt}:=\frac{\partial}{\partial{t}}+u\cdot\nabla
$$

### Navier-Stokes equations (Incompressible)

$$
\begin{align}
\rho\frac{Du}{Dt}&=-\nabla{p}+\mu\nabla^2u+\rho{g}\\
\frac{Du}{Dt}&=-\frac{1}{\rho}\nabla{p}+v\nabla^2u+g\\
\nabla\cdot{u}&=0
\end{align}
$$

###### Operator splitting

$$
\begin{align}
\frac{Du}{Dt}&=-\frac{1}{\rho}\nabla{p}+g\\
\frac{Du}{Dt}&=0,\frac{D\alpha}{Dt}=0\\
\frac{\partial{u}}{\partial{t}}&=g\\
\frac{\partial{u}}{\partial{t}}&=-\frac{1}{\rho}\nabla{p}\ \ \ \ \ s.t.\ \ \ \ \nabla\cdot{u}=0
\end{align}
$$

###### Eulerian fluid simulation cycle

Time discretization with splitting: for each time step,

1. Advection: "move" the fluid field. Solve u* using u^t

   $$
   \frac{Du}{Dt}=0,\frac{D\alpha}{Dt}=0
   $$

2. External forces (optional): evaluate u\** using u*

   $$
   \begin{align}
   \frac{\partial{u}}{\partial{t}}=g
   \end{align}
   $$

3. Projection: make velocity field u^{t+1} divergence-free based on u\**

   $$
   \frac{\partial{u}}{\partial{t}}=-\frac{1}{\rho}\nabla{p}\ \ \ \ \ s.t.\ \ \ \ \nabla\cdot{u^{t+1}}=0
   $$



### BFECC & MacCormack

###### Back and Forth Error Compensation and Correction

$$
\begin{align}
x^*&=SL(x,\Delta t)\\
x^{**}&=SL(x^*,-\Delta t)\\
x^{error}&=\frac{1}{2}(x^{**}-x)\\
x^{final}&=x^*-x^{error}
\end{align}
$$

### Grid

###### Spatial discretization using cell-centered grids

###### Spatial discretization using staggered grids

###### Bilinear interpolation

### Advection schemes

1. Semi-Lagrangian advection
2. MacCormack/BFECC
3. "BiMocq^2"
4. Particle advection (PIC/FLIP/APIC/PolyPIC)
5. ...

### Projection

Chorin-style projection

$$
\begin{align}
u*-u&=-\Delta t
\end{align}
$$

# Solving large-scale linear systems

$$
Ax=b
$$

###### How to solve it

1. Direct solvers (e.g., PARDISO)

2. Iterative solvers

   Gauss-Seidel

   (Damped) Jacobi

   (Preconditioned) Krylov-subspace solvers (e.g., conjugate gradients)

###### How to store A

1. As a dense matrix (e.g., float A\[1024]\[1024] doesn't scale but works)
2. As a sparse matrix (various sparse matrix formats: CSR, COO, ...)
3. Don't store it at all (aka. Matrix-free, often the ultimate solution...)

Modern computer architecture: memory bandwidth is expensive but FLOPs are free. So compute matrix entries on-the-fly (instead of fetching valeus from memory) can sometimes be good to performance.

### Krylov-subspace solvers

Krylov-subspace solvers are among most efficient linear system solvers. The most well-known version: conjugate gradients (CG).

Less frequently used (in graphics):

1. Conjugate residuals (CR)
2. Generalized minimal residual method (GMRES)
3. Biconjugate gradient stabilized (BiCGStab)

###### Conjugate gradients

###### Eigenvalues and condition numbers

$$
\begin{align}
Ax&=\lambda x\\
\kappa(A)&=\lambda_{max}/\lambda_{min}
\end{align}
$$

In general: a smaller condition numbers means faster convergence

(Note that condition numbers have many different definitions)

### Iterative solver trick: Warm starting

If you start with an initial guess that is close to the solution, very likely fewer iterations are needed.

F11"Warm starting": use the p from last frame as the initial guess of the current frame.

In practice works well for (damped) Jacobi/Gauess-Seidel/CG, but for MGPCG it doesn't work well.

### Preconditioning

$$
\begin{align}
Ax&=b\\
M^{-1}Ax&=M^{-1}b
\end{align}
$$

Intuition: M^{-1}A may have a smaller condition number (closer to identity) or better eigenvalue clustering than A itself.

###### Common preconditioners

1. Jacobi (diagonal) preconditioner M = drag(A)

2. Poisson preconditioner

3. (Incomplete) Cholesky decomposition

4. Multigrid: M = very complex linear operator that almost inverts A

   Geometirc multigrid

   Algebraic multigrid

5. Fast multipole method (FMM)

### Multigrid methods (Geometric)

Multigrid V-Cycle: Solving PHI in PDE f(PHI) = F

The Multigrid design space

1. Restriction/prolongation
2. Cycle (V/W/F cycle)
3. Smoothers ((red-black) Gauss-Seidel, Jacobi, damped Jacobi, etc.)
4. Number of levels (e.g., coarsen until the bottom level has < 50k voxels)
5. Bottom level solver (Brute-force Jacobi or direct solvers)
6. Number of pre/post iterations (usually, 2-5)
7. Coarsening and boundary handling (e.g., Galerkin coarsening, semi-algebraic multigrid)
8. ...

### Multigrid preconditioned conjugate gradients (MGPCG)

When solving Poisson's equation in graphics, people usually use geometric multigrid as the preconditioner for conjugate gradients.

# 代码结构

```python
def step(mouse_data):
    advect(velocities_pair.cur, velocities_pair.cur, velocities_pair.nxt)
    advect(velocities_pair.cur, dyes_pair.cur, dyes_pair.nxt)
    velocities_pair.swap()  # 速度场交换新旧缓冲
    dyes_pair.swap()  # 颜色场交换新旧缓冲
    apply_impulse(velocities_pair.cur, dyes_pair.cur, mouse_data)
    divergence(velocities_pair.cur)
    if curl_strength:
        vorticity(velocities_pair.cur)
        enhance_vorticity(velocities_pair.cur, velocity_curls)
    if use_sparse_matrix:
        solve_pressure_sp_mat()
    else:
        solve_pressure_jacobi()
    subtract_gradient(velocities_pair.cur, pressures_pair.cur)
    if debug:
        divergence(velocities_pair.cur)
        div_s = np.sum(velocity_divs.to_numpy())
        print(f"divergence={div_s}")
```

## Advection 流动

###### advect

根据上一帧速度场 velocity_pair.cur 模拟流体流动，得到当前帧初始速度场 velocities_pair.nxt 和初始颜色场 dyes_pair.nxt

对于每个位置p，根据位置p的速度，反向寻找位置 p' = p - dt * v_p ，使用位置 p' 的场的值来作为位置 p 的场的新值，这种方法称为半拉格朗日法 semi_lagrangian

```python
@ti.kernel
def advect(vf: ti.template(), qf: ti.template(), new_qf: ti.template()):
    for i, j in vf:
        p = ti.Vector([i, j]) + 0.5
        p = backtrace(vf, p, dt)
        new_qf[i, j] = bilerp(qf, p) * dye_decay
```

同时，根据阶数不同可分为 RK1 RK2 RK3

1. Runge Kutta 1 前向欧拉法

   ```python
   p -= dt * velocity(p)
   ```

2. Runge Kutta 2 显式中点法

   ```python
   p_mid = p - 0.5 * dt *velocity(p)
   p -= dt * velocity(p_mid)
   ```

3. Runge Kutta 3

   ```python
   v1 = velocity(p)
   p1 = p - 0.5 * dt * v1
   v2 = velocity(p1)
   p2 = p - 0.75 * dt * v2
   v3 = velocity(p2)
   ```

4. BFECC

###### apply_impulse

根据鼠标输入，对速度场和颜色场中某些位置进行设置

```python
@ti.kernel
def apply_impulse(vf: ti.template(), dyef: ti.template(), imp_data: ti.types.ndarray()):
    g_dir = -ti.Vector([0, 9.8]) * 300
    for i, j in vf:
        omx, omy = imp_data[2], imp_data[3]
        mdir = ti.Vector([imp_data[0], imp_data[1]])
        dx, dy = (i + 0.5 - omx), (j + 0.5 - omy)
        d2 = dx * dx + dy * dy
        # dv = F * dt
        factor = ti.exp(-d2 / force_radius)
        dc = dyef[i, j]
        a = dc.norm()
        momentum = (mdir * f_strength * factor + g_dir * a / (1 + a)) * dt
        v = vf[i, j]
        vf[i, j] = v + momentum
        # add dye
        if mdir.norm() > 0.5:
            dc += ti.exp(-d2 * (4 / (res / 15) ** 2)) * ti.Vector([imp_data[4], imp_data[5], imp_data[6]])
        dyef[i, j] = dc
```

## Projection 投影

对应 divergence pressure_jacobi subtract_gradient 三个步骤。模拟时一般将流体视为不可压缩流体，在经过 advect 和 apply_impulse 步骤后，流体的速度场处于一个不平衡的状态，即“流入 != 流出”，而 Projection 阶段的目的是重新实现“流入=流出”，即不可压缩性。

计算思路：对“流入!=流出”状态的流体求各个位置的速度的散度，通过速度的散度更新压强，根据压强的梯度对流体进行加速以更新流体的速度

$$
\nabla\cdot\nabla p =\frac{\rho}{\Delta t}\nabla\cdot u\\
\nabla\cdot求散度\\
\nabla p \ 压强的梯度
$$

###### divergence

根据速度场计算每个位置的速度散度 divergence 函数为求速度的散度，公式如下

$$
(\frac{\rho}{\Delta t}\nabla\cdot u)_{ij}=\frac{\rho}{2\Delta t\Delta x}(u^x_{i+1,j}-u^x_{i-1,j}+u^y_{i,j+1}-u^y_{i,j-1})
$$

```python
@ti.kernel
def divergence(vf: ti.template()):
    for i, j in vf:
        vl = sample(vf, i - 1, j)
        vr = sample(vf, i + 1, j)
        vb = sample(vf, i, j - 1)
        vt = sample(vf, i, j + 1)
        vc = sample(vf, i, j)
        if i == 0:
            vl.x = -vc.x
        if i == res - 1:
            vr.x = -vc.x
        if j == 0:
            vb.y = -vc.y
        if j == res - 1:
            vt.y = -vc.y
        velocity_divs[i, j] = (vr.x - vl.x + vt.y - vb.y) * 0.5
```

###### pressure_jacobi

根据速度的散度更新各个位置的压力，多次迭代贴近平衡状态 pressure_jacobi 函数为从速度的散度对压强进行更新，公式如下

$$
\begin{align}
(\nabla\cdot\nabla p)_{i,j}&=\frac{1}{\Delta x^2}(-4p_{i,j}+p_{i+1,j}+p_{i-1,j}+p_{i,j+1}+p_{i,j-1})\\
p_{i,j}&=\frac{1}{4}(p_{i+1,j}+p_{i-1,j}+p_{i,j+1}+p_{i,j-1}+\Delta x^2(\nabla\cdot\nabla p)_{i,j})
\end{align}
$$
由于每个位置的压强与其周围的压强相关，因此需要多次迭代以接近稳定状态

```python
@ti.kernel
def pressure_jacobi(pf: ti.template(), new_pf: ti.template()):
    for i, j in pf:
        pl = sample(pf, i - 1, j)
        pr = sample(pf, i + 1, j)
        pb = sample(pf, i, j - 1)
        pt = sample(pf, i, j + 1)
        div = velocity_divs[i, j]
        new_pf[i, j] = (pl + pr + pb + pt - div) * 0.25
```

###### subtract_gradient

压强差会对流体施加力，使流体具有加速度，使用加速度来更新流体的速度 subtract_gradient 函数为使用梯度来更新速度，公式如下

$$
\begin{align}
\frac{\partial u}{\partial t}=-\frac{\Delta p}{\rho}\\
u_{t+1}=u_t-\Delta t\frac{\Delta p}{\rho}\\
\Delta p\ 周围各方向压强差
\end{align}
$$

```python
@ti.kernel
def subtract_gradient(vf: ti.template(), pf: ti.template()):
    for i, j in vf:
        pl = sample(pf, i - 1, j)
        pr = sample(pf, i + 1, j)
        pb = sample(pf, i, j - 1)
        pt = sample(pf, i, j + 1)
        vf[i, j] -= 0.5 * ti.Vector([pr - pl, pt - pb])
```

