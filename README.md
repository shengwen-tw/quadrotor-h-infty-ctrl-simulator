# Quadrotor H∞ Control Simulator

> Check our Python implementation here: https://github.com/shengwen-tw/rotor-bench 

H∞ control is a robust control strategy that emerged in the 1970s. It focuses on minimizing the worst-case impact of disturbances on system performance by suppressing the maximum singular value of the transfer function from disturbance to output over all frequencies. This ensures that the system remains stable and performs acceptably even under model uncertainties and external disturbances. The optimal control law is derived by solving an Algebraic Riccati Equation, which is closely linked to the Hamiltonian matrix.

This project provides a design example of using H∞ control to stablize the quadrotor. We also feature that the code is fully standalone without relying on other control toolbox.

<img src="https://github.com/shengwen-tw/quadrotor-h-infty-ctrl-simulator/blob/master/simulation.png?raw=true" width="90%" height="90%">


## Run simulation

```
execute quadrotor_sim.m
```

## Code overview

``quadrotor_sim.m``: The main program of the simulation.

``dynamics.m``: The quadrotor dynamics model.

``hinf_norm.m``: The function for calculating H∞ norm.

``rigidbody_visualize.m``: The 3D animation for quadrotor visualization.

``se3_math.m``: Some math functions for calculating rotation matrix, quaternion and Euler angles.

``care_sda.m``: Structure-Preserving Doubling Algorithm (SDA), a numerical method for solving Algebraic Riccati Equation.

``hinf_syn.m``: A H∞ control solver based on the bisection method and SDA.

## Publication

This work is published at the 2022 International Automatic Control Conference (CACS) with the title of [Robust State-Feedback H∞ Control of Quadrotor](https://ieeexplore.ieee.org/document/9969787), please feel free to cite the paper.

## References


1. [Geometric tracking control of a quadrotor UAV on SE(3)](https://ieeexplore.ieee.org/document/5717652)
2. [Quadrotor control: modeling, nonlinear
control design, and simulation](https://www.kth.se/polopoly_fs/1.588039.1600688317!/Thesis%20KTH%20-%20Francesco%20Sabatino.pdf)
3. [A Bisection Method for Computing the H∞ Norm of a Transfer Matrix and Related Problems
](https://web.stanford.edu/~boyd/papers/bisection_hinfty.html)
4. [On the computation of the optimal H∞ norms for two feedback control problems](https://www.infona.pl/resource/bwmeta1.element.elsevier-93e70f8a-b5a8-37fe-a5b2-c6f8551ddfdb)
5. [A structure-preserving doubling algorithm for
continuous-time algebraic Riccati equations](https://jupiter.math.nycu.edu.tw/~wwlin/papers_new/2005CFLp.pdf)
