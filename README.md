# LongRangeIsingVQE
## Author: [Aniruddha Bapat](https://github.com/quristic),  Stephen Jordan
### Date: 07/20/2018

This repository contains code that can be used to simulate the evolution of a quantum state on n qubits under
two possible Hamiltonians:

1. The phase Hamiltonian: This Hamiltonian is a sum of Pauli ZZ and Pauli Z terms on the qubit indices, with
arbitary coefficients.
2. The mixer: This is fixed to be the negative sum of 1-qubit Pauli X operators over all qubit indices

By default, the above two steps are run sequentially in order per iteration. The user can specify the number of 
iterations as a parameter p. The exact schedule (i.e. the evolution times) may be specified by the user. Or,
you can use some in-built functions to find the best schedule.

The "best" schedule usually means one which minimizes the energy expectation of the final state with respect to 
a problem Hamiltonian. The problem Hamiltonian can be specify as a sum of ZZ, Z, and X coefficients terms with
arbitrary coefficients. 

Other figures of merit, such as ground state overlap, may be used instead. 

Finally, there is a function to compute the entanglement entropy with respect to a given cut in the spin chain. Note
that all of the above computations use exact diagonalization, so they're exact but potentially not scalable. 
