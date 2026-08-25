# Schwinger Phase Diagram
Phase diagram of the Schwinger model in the presence of a topological term using the [ITensors](https://itensor.org) package for the [Julia Programming Language](https://julialang.org/).

## Overview

This repository accompanies the hands-on tutorials for the lectures on Tensor Networks at the [CERN-MPQ-UIBK School on Quantum Simulation of Fundamental Physics](https://indico.cern.ch/event/1623729/). 

The project provides an example of exploring a simple lattice gauge theory with Matrix Product States. It follows the exercises from the accompanying tutorial. Both cases, integrating the gauge fields out and truncating them to a finite dimension are treated.

## Repository Structure

```
project/
│
├── src/            Source code
├── examples/       Example programs
├── LICENSE
└── README.md
```
The folder `src` contains the codes and different functions that are created by solving the problems on the exercise sheet. The folder `examples` contains a simple example of running a variational optimization for MPS to explore the phase diagram of the Schwinger model


## Installation

To run the code, you will need to have Julia installed. Julia can be downloaded from the [official website](https://julialang.org/downloads/). 

For instructions on getting started, see the [Julia documentation](https://docs.julialang.org/en/v1/manual/getting-started/).

For instructions on how to install the [ITensors](https://docs.itensor.org/ITensors/stable/) and [ITensorMPS](https://docs.itensor.org/ITensorMPS/stable/) packages follow the links to the documentaiton.


## Contributing

Contributions are welcome.

Please
- open an issue for bugs,
- discuss major changes before implementing them,
- follow the existing coding style.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## References

There is a vast amount of literature on MPS and more general tensor network methods available. The following incomplete list of reviews and references therein might provide a useful starting point for learning more about MPS and more general tensor networks.

* S. Coleman, [Ann. Phys. 101, 239–267 (1976)](https://doi.org/10.1016/0003-4916(76)90280-3)
* C. J. Hamer, Z. Weihong, and J. Oitmaa, [Phys. Rev. D56, 55 (1997)](https://doi.org/10.1103/PhysRevD.56.55) 
* T. M. R. Byrnes, P. Sriganesh, R. J. Bursill, and C. J. Hamer, [Phys. Rev. D66, 013002 (2002) ](https://doi.org/10.1103/PhysRevD.66.013002)
* L. Funcke, K. Jansen, and S. Kühn, [Phys. Rev. D101, 054507 (2020)](https://doi.org/10.1103/PhysRevD.101.054507)

## Authors

[Stefan Kühn](https://github.com/kuehnste)
