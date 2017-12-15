# restrictedCircular3BodyProblem

Toolbox for computation of orbits in the RC3BP. 
MATLAB OOP implementation of primaries system, numerical methods for orbit propagation, state transition matrices.
Abstraction of targeting algorithms and continuation methods for the computation of orbits and families of orbits are provided.
Orbit stability analysis can be run directly on top of this toolbox.

Three tests are provided:
1. Computation of family of Lyapunov orbits starting from L1 linearized orbit
1. Computation of family of Halo orbits from known ICs and propagating it into the Lyapunov orbit bifurcation
1. Computation of maps to identify periodic and quasi-periodic motion

![Family of Lyapunov and Halo Orbits](/results/LyapunovHaloFamilies.png)
![Map for the Earth Moon system with L1 gate open for y=0](/results/EarthMoonL1Map.png)
