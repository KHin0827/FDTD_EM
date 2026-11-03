# Yee Scheme for the 2D TM Electromagnetic Wave Model

A MATLAB implementation of the **Yee finite-difference time-domain (FDTD) scheme** for the 2D transverse-magnetic (TM) Maxwell system with perfect-electric-conductor (PEC) boundary conditions.

The entire project is contained in a single driver script, `FDTD.m`, which runs two experiments:

1. **Accuracy and convergence test** against a manufactured exact solution.
2. **Wave propagation from a localized Gaussian source**, with snapshots of $E_z$ at several times.

---

## Model Problem

The dimensionless 2D TM system is

$$
\begin{aligned}
\frac{\partial E_z}{\partial t} &= \frac{\partial H_y}{\partial x} - \frac{\partial H_x}{\partial y} - \sin(\omega t)J_z, \\
\frac{\partial H_x}{\partial t} &= -\frac{\partial E_z}{\partial y}, \\
\frac{\partial H_y}{\partial t} &= \frac{\partial E_z}{\partial x},
\end{aligned}
$$

with PEC boundary condition

$$
E_z = 0 \quad \text{on } \partial\Omega.
$$

The unknowns are placed on a staggered Yee grid:

- $E_z$ at integer nodes $(i\Delta x,\,j\Delta y)$.
- $H_x$ at $(i\Delta x,\,(j+\tfrac12)\Delta y)$.
- $H_y$ at $((i+\tfrac12)\Delta x,\,j\Delta y)$.

---
The explicit update implemented in `FDTD.m` is

$$
\frac{(E_z)_{i,j}^{n+1} - (E_z)_{i,j}^{n}}{\Delta t}=\frac{(H_y)_{i+\frac12,j}^{n+\frac12} - (H_y)_{i-\frac12,j}^{n+\frac12}}{\Delta x}-
\frac{(H_x)_{i,j+\frac12}^{n+\frac12} - (H_x)_{i,j-\frac12}^{n+\frac12}}{\Delta y}-
\sin\\left(\omega t^{n+\frac12}\right)(J_z)_{i,j},
$$

$$
\frac{(H_x)_{i,j+\frac12}^{n+\frac12} - (H_x)_{i,j+\frac12}^{n-\frac12}}{\Delta t}=-\frac{(E_z)_{i,j+1}^{n} - (E_z)_{i,j}^{n}}{\Delta y},
\\
\frac{(H_y)_{i+\frac12,j}^{n+\frac12} - (H_y)_{i+\frac12,j}^{n-\frac12}}{\Delta t}=\frac{(E_z)_{i+1,j}^{n} - (E_z)_{i,j}^{n}}{\Delta x}.
$$

On the boundary $E_z$ is set directly to zero; $H_x$ and $H_y$ on the boundary are updated by the same stencil. The CFL condition for 2D is

$$
\Delta t \le \frac{1}{c}\left(\frac{1}{\Delta x^{2}} + \frac{1}{\Delta y^{2}}\right)^{-1/2}.
$$

The code uses

$$
\Delta t = a\Delta x, \qquad a = 0.25,
$$

---

## Accuracy and Convergence

The manufactured source is

$$
J_z =
16\omega x^{2}(x-1)^{2}y^{2}(y-1)^{2}
+
\frac{32}{\omega}\Big[
(6x^{2}-6x+1) y^{2}(y-1)^{2}
+
(6y^{2}-6y+1) x^{2}(x-1)^{2}
\Big],
$$

with the exact solution

$$
E_z(x,y,t) = 16x^{2}(x-1)^{2}y^{2}(y-1)^{2}\cos(\omega t),
$$

which satisfies the PEC condition on the boundary of the unit square.

**Parameters used in the code**

| Quantity | Value |
|---|---|
| Domain | unit square $[0,1]^{2}$ |
| Final time | $T = 0.8$ |
| $\Delta t$ | $0.25\,\Delta x$, then adjusted to hit $T$ |
| Frequencies | $\omega = 25,\ 50,\ 100$ |
| Meshes | $N = 40,\ 80,\ 160,\ 320,\ 640,\ 1280$ |
| Initial $E_z$ | $16x^{2}(x-1)^{2}y^{2}(y-1)^{2}$ |
| Initial $H_x, H_y$ | $0$ |

For each $(\omega, N)$ the code computes the maximum error at $T = 0.8$ over interior grid points:

$$
\mathrm{err}(N) = \max
\left| (E_z)_{i,j}^{M} - E_z(x_i, y_j, T) \right|.
$$

The results are plotted on a log–log scale with $N$ on the horizontal axis. For reference, one line of slope $-1$ and one of slope $-2$ are drawn next to the data. The code also prints the observed convergence rates

$$
\text{rate}_{k} = \frac{\log \left(\mathrm{err}(N_k)/\mathrm{err}(N_{k+1})\right)}{\log 2},
$$

which should be close to $2$ once the mesh is fine enough.

**Expected behaviour:** The Yee scheme is second-order accurate in space and time, so $\mathrm{err}(N) \sim \mathcal{O}(N^{-2})$ when $N$ is large enough. On coarse meshes and high $\omega$, the scheme shows the expected pre-asymptotic effects from dispersion and phase error.

---

## Localized Gaussian Source

The second experiment uses the source

$$
J_z = -\omega \exp\Big(-\sigma\big[(x-0.01)^{2} + (y-0.015)^{2}\big]\Big),
\qquad
\sigma = \max(36, \omega^{2}),
$$

on the domain $[-1,1]^{2}$ with zero initial data:

$$
E_z\big|_{t=0} = H_x\big|_{t=0} = H_y\big|_{t=0} = 0.
$$

**Parameters used in the code**

| Quantity | Value |
|---|---|
| Domain | $[-1,1]^{2}$ |
| Final time | $T = 10$ |
| $\Delta t$ | $0.25\,\Delta x$, then adjusted to hit $T$ |
| Frequencies | $\omega = 10,\ 50,\ 100$ |
| Meshes | $N = 40,\ 80,\ 160$ |
| Snapshot times | $t = 2.5,\ 5,\ 7.5,\ 10$ |

For each $\omega$ the code produces one figure with a $3 \times 4$ grid of subplots: rows correspond to $N \in \{40, 80, 160\}$, and columns to the four snapshot times. Each panel shows $E_z(x, y, t)$ over the domain.

**Observations**

- **Circular radiation and interference.** The source \(J_z\) centred at (0.01,0.015) excites circularly spreading wavefronts; once they reach the PEC walls they reflect, and the later snapshots show interference between the outgoing and reflected waves. \(E_z\) stays zero on the boundary throughout, as required by the PEC condition.

- **Numerical dispersion.** On coarse meshes the wavefront develops trailing ripples and a phase error. This is the classical FDTD numerical dispersion: the discrete dispersion relation differs from the exact, and the error grows with the grid spacing measured in wavelengths.

- **Frequency dependence.** Larger $\omega$ means shorter $\lambda$ , so more oscillations are packed into the same grid. At $\omega=10$ the long wavelength is well resolved and the pattern is clearly visible; at $\omega=50$ and 100 the wavelength becomes comparable to $\Delta x$, the field looks like a uniform green background, and the interference detail is hard to see.

- **Resolution dependence.** Increasing $N$ from 40 to 160 reduces $\Delta x/\lambda$; the numerical solution becomes more realistic — wavefronts sharpen, ringing decreases, and the phase error diminishes, so the high $\omega$ cases are visibly better resolved.

