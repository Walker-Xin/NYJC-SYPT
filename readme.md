# NYJC Singapore Young Physicist Tournament 2021

## Introduction

This repository hosts the algorithms developed for Singapore Young Physicist Tournament (SYPT) 2021 by me. My team investigated problems Bead Dynamics, Wilberforce Pendulum and Dynamic Hydrophobicity. Numerical simulations are produced for the first two problems and the script for Dynamic Hydrophobicity is only used for visualising data collected. 

## Algorithms

For the mechanics problems, the scripts include solution of the equations of motion by either numerical methods (Bead Dynamics) or analytical solution (Wilberforce Pendulum). Energy profile of the system is analysed to compute the effect of friction. Visualisation and animation using Matplotlib is added whenever possible to provide a better understanding of the system.

For Dynamic Hydrophobicity, the script is largely dedicated to visualising the experimental data we collected. Some phase diagrams are drawn to demonstrate the relationship between local and impact speeds.

The script for Spin Drift is incomplete as I could not find the equations of motion of the system. A wonderful solution and simulation for the flat ground case, however, can be found at https://rotations.berkeley.edu/the-rolling-disk/.

## References for Further Reading

### Wilberforce Pendulum
1. Berg, R. E., & Marshall, T. S. (1991). Wilberforce pendulum oscillations and normal modes. American Journal of Physics, 59(1), 32-38.
2. Zanette, D. H. (2018). Energy exchange between coupled mechanical oscillators: linear regimes. Journal of Physics Communications, 2(9), 095015.
3. Landau, L. D. (1973). Small Oscillations. Classical mechanics. Theoretical Physics, 1.
4. Dolfo, G., & Vigué, J. (2018). Damping of coupled harmonic oscillators. European Journal of Physics, 39(2), 025005.

### Bead Dynamics
1. Dutta, S., & Ray, S. (2011). Bead on a rotating circular hoop: a simple yet feature-rich dynamical system. arXiv preprint arXiv:1112.4697.
2. Raviola, L. A., Véliz, M. E., Salomone, H. D., Olivieri, N. A., & Rodríguez, E. E. (2016). The bead on a rotating hoop revisited: an unexpected resonance. European Journal of Physics, 38(1), 015005.
3. Cross, R. (2016). Coulomb's law for rolling friction. American Journal of Physics, 84(3), 221-230.

### Dynamic Hydrophobicity
1. Gauthier, A., Bouillant, A., Clanet, C., & Quéré, D. (2018). Aerodynamic repellency of impacting liquids. Physical Review Fluids, 3(5), 054002.
2. Landau, L., & Levich, B. (1988). Dragging of a liquid by a moving plate. In Dynamics of curved fronts (pp. 141-153). Academic Press.
3. Siddiqui, M. E., Mukund, V., Scott, J., & Pier, B. (2013). Experimental characterization of transition region in rotating-disk boundary layer. Physics of Fluids, 25(3), 034102.
4. Gauthier, A., Bird, J. C., Clanet, C., & Quéré, D. (2016). Aerodynamic Leidenfrost effect. Physical Review Fluids, 1(8), 084002.
