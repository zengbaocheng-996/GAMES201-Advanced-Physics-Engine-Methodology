# Lagrangian View

### Time integration

1. Explicit (forward Euler, symplectic Euler, RK, ...)
   
   $$
   \begin{align}
   v_{t+1}&=v_t+\Delta{t}\frac{f_t}{m}\\
   x_{t+1}&=x_t+\Delta{tv_t}\\
   \Delta{t}&\leq{c}\sqrt{\frac{m}{k}}
   \end{align}
   $$

2. Semi-implicit Euler (aka. symplectic Euler, explicit)
   
   $$\\
   \begin{align}
   v_{t+1}&=v_t+\Delta{t}\frac{f_t}{m}\\
   x_{t+1}&=x_t+\Delta{tv_{t+1}}
   \end{align}\\
   $$

3. Implicit (backward Euler, middle-point, ...)

### Mass-spring systems

$$
\begin{align}
f_{ij}&=-k(||x_i-x_j||_{2}-l_{ij})(\widehat{x_i-x_j})\\
f_i&=\sum^{j\neq{i}}_jf_{ij}\\
\frac{\partial{v_i}}{\partial{t}}&=\frac{1}{m_i}f_i\\
\frac{\partial{x_i}}{\partial{t}}&=v_i
\end{align}
$$

###### Implicit time integration

$$
\begin{align}
x_{t+1}&=x_t+\Delta{tv_{t+1}}\\
v_{t+1}&=v_t+\Delta{t}M^{-1}f(x_{t+1})\\
v_{t+1}&=v_t+\Delta{t}M^{-1}f(x_t+\Delta{tv_{t+1}})\\
v_{t+1}&=v_t+\Delta{t}M^{-1}[f(x_t)+\frac{\partial{f}}{\partial{x}}(x_t)\Delta{tv_{t+1}})]\\
[I-\Delta{t^2M^{-1}}\frac{\partial{f}}{\partial{x}}(x_t)]v_{t+1}&=v_t+\Delta{t}M^{-1}f(x_t)
\end{align}
$$

### Solver

$$
\begin{align}
A&=[I-\Delta{t^2M^{-1}}\frac{\partial{f}}{\partial{x}}(x_t)]\\
b&=v_t+\Delta{t}M^{-1}f(x_t)\\
Av_{t+1}&=b
\end{align}
$$

$$
[I-\beta\Delta{t^2M^{-1}}\frac{\partial{f}}{\partial{x}}(x_t)]v_{t+1}=v_t+\Delta{t}M^{-1}f(x_t)
$$

1. β=0 forward/semi-implicit Euler (explicit)
2. β=1/2 middle-point (implicit)
3. β=1 backward Euler (implicit)

###### Jacobi iterations

```python
import taichi as ti
import random

ti.init()

n = 20

A = ti.field(dtype=ti.f32, shape=(n, n))
x = ti.field(dtype=ti.f32, shape=n)  # 未知数
new_x = ti.field(dtype=ti.f32, shape=n)
b = ti.field(dtype=ti.f32, shape=n)  # 右端项


@ti.kernel
def iterate():
    for i in range(n):
        r = b[i]
        for j in range(n):
            if i != j:
                r -= A[i, j] * x[j]
        new_x[i] = r / A[i, i]
    for i in range(n):
        x[i] = new_x[i]


@ti.kernel
def residual() -> ti.f32:
    res = 0.0
    for i in range(n):
        r = b[i] * 1.0
        for j in range(n):
            r -= A[i, j] * x[j]
        res += r * r
    return res


for i in range(n):
    for j in range(n):
        A[i, j] = random.random() - 0.5
    A[i, i] += n * 0.1
    b[i] = random.random() * 100

for i in range(66):
    iterate()
    print(f'iter {i} {residual():.10f}')

for i in range(n):
    lhs = 0.0
    for j in range(n):
        lhs += A[i, j] * x[j]
    assert abs(lhs - b[i]) < 1e-4
```
