# Numerical Schemes for Gradient-Flow Models

This repository contains MATLAB implementations of numerical schemes for gradient-flow equations, based on two published works. The focus is on reproducing and comparing time discretization methods for phase-field models.

---

## Models

- Allen–Cahn equation (1D)
- Cahn–Hilliard equation (2D)

---

## Implementations

### 1. Allen–Cahn Equation (1D)

File:
AllenCahn/AC_IMEX_1D_TangYang2016.m

This script implements the implicit–explicit (IMEX) schemes proposed in:

T. Tang and J. Yang,  
Implicit–Explicit Scheme for the Allen–Cahn Equation Preserves the Maximum Principle,  
Journal of Computational Mathematics, 2016.

Schemes included:

- Scheme (2.6): standard IMEX method  
- Scheme (3.1): stabilized IMEX scheme with parameter β  

Features:

- Finite difference discretization with Neumann boundary conditions  
- Maximum norm evolution  
- Discrete energy decay  
- Stability comparison under different time step sizes and β values  

---

### 2. Cahn–Hilliard Equation (2D)

File:
CahnHilliard/CH_ETDRK2_2D_FuYang2022.m

This script implements a Fourier spectral method combined with the ETDRK2 scheme, following:

X. Fu and J. Yang,  
Energy-decreasing exponential time differencing Runge–Kutta methods for phase-field models,  
Journal of Computational Physics, 2022.

Features:

- Fourier spectral discretization (periodic boundary conditions)  
- ETDRK2 time integration  
- Snapshot visualization at multiple time instances  
- Free-energy evolution  

Numerical setup:

- Domain: (0, 2π) × (0, 2π)  
- Initial condition: u₀ = 0.05 sin(x) sin(y)  
- Grid size: N = 512  
- Parameter: ε = 0.1  
- Time step: Δt = 1e-3  
- Final time: T = 8  

---

## Language

- MATLAB

---

## Notes

This repository reproduces key numerical results and stability properties from the referenced papers. The codes are intended for numerical experiments and comparison of time-stepping schemes for gradient-flow PDEs.
