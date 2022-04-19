# Quadrotor H∞ Control Simulator

H∞ control is a robust control algorithm arised in 1970s. It suppresses the maximum frequency response of the disturbance to ensure robustness of the system. Mathematicians shows that the optimal control input can be obtained by solving the Algebraic Riccati Equation, which is strongly related to the Hamiltonian matrix.

This project provides a design example of using H∞ control to stablize the quadrotor. We also feature that the code is fully standalone without relying on other control toolbox.

## Run simulation

```
execute quadrotor_sim.m
```

## Code overview

``quadrotor_sim.m``: Main program of the simulation.

``dynamics.m``: The code implements the dynamics model of the quadrotor.

``hinf_norm.m``: The code to calculate the H-infinity norm of a state-space system.

``rigidbody_visualize.m``: A 3D animation utility for quadrotor visualization.

``se3_math.m``: Some math functions for calculating the rotation matrix, quaternion and Euler angles.

``care_sda.m``: A fast Continuous-time Algebraic Riccati Equation (CARE) solver based on the Structure-Preserving Doubling Algorithm (SDA).

``hinf_syn.m``: A H-infinity control synthesizer based on the bisection method and the SDA.

## References


1. [Geometric tracking control of a quadrotor UAV on SE(3)](https://ieeexplore.ieee.org/document/5717652)
2. [Quadrotor control: modeling, nonlinear
control design, and simulation](https://www.kth.se/polopoly_fs/1.588039.1600688317!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf)
3. [A Bisection Method for Computing the H∞ Norm of a Transfer Matrix and Related Problems
](https://web.stanford.edu/~boyd/papers/bisection_hinfty.html)
4. [On the computation of the optimal H∞ norms for two feedback control problems](https://www.infona.pl/resource/bwmeta1.element.elsevier-93e70f8a-b5a8-37fe-a5b2-c6f8551ddfdb)
5. [A structure-preserving doubling algorithm for
continuous-time algebraic Riccati equations](https://jupiter.math.nycu.edu.tw/~wwlin/papers_new/2005CFLp.pdf)
