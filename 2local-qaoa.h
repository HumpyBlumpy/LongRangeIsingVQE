/*--------------------------------------------------------------------
  Simulator and data analysis tools for Quantum Approximate 
  Optimization Algorithm and Variational Quantum Eigensolvers, written
  by Stephen Jordan in 2018. Microsoft QuARC group.

  Modified (05/03/2018) by Aniruddha Bapat to apply to 2-local Ising 
  Hamiltonians with arbitary ZZ, X, and Z couplings. 
  -------------------------------------------------------------------*/

#define PI 3.14159265359

//Stores the quantum state vector
typedef struct {
  int n;            //the number of qubits
  int N;            //the dimension of the Hilbert space (=2^n)
  double *realcur;  //real parts of amplitudes of current state
  double *imagcur;  //imaginary parts of amplitudes of current state
  double *realbuf;  //real buffer
  double *imagbuf;  //imaginary buffer
}state;

//Stores the Hamiltonian
typedef struct {
  int n;            //the number of qubits
  int N;            //the dimension of the Hilbert space (=2^n)
  double *zzc;  //2D array that stores Jij
  double *xc;  //1D array that stores Bi
  double *zc;  //1D array that stores Ki (1-local Z coefficients)
  double *energy; // 1D array that stores the Z energy for each bistring 
}ham;

//Print out the integer x in binary
void printbits(int x, int n);

//Print out the quantum state vector
void printvec(state psi);

//Return the Euclidean norm of the state vector
double norm(state psi);

//Swap the current state with the buffered state
void swapbuf(state *psi);

//Perform a Hadamard transform on the n qubits in O(NlogN) time.
void Hadamard(state *psi);

//Computes the Z2 inner product between bitstrings a and b of length n.
int Z2inner(int A, int B, int n);

//Convert an array of 1s and 0s to the corresponding number via binary place value.
//The zeroth element of the array contains the 1s place, the 1th element the 2s etc.
int b2i(int *bitstring, int n);

//Convert a string of characters like "1100" to the corresponding number like 12.
//Note that in keeping with standard conventions, this is the reverse order for
//significant bits than in bin2int.
int cb2i(char *bitstring, int n);

//Returns the number of ones in the binary expansion of x
int HammingWeight(int x, int n);

//Returns a skew weight weighted by bi
double SkewWeight(int x, int n, double *b);

//Unitarily evolve for time t according to H = -sum_j X_j
void evolveX(state *psi, double t);

//The energy associated with a bit string according to Hamiltonian H
//H stores hamiltonian coefficients
//x stores the bit string
double energy(int x, ham *H);

//Unitarily evolve for time t according to H = \sum_i<j J_ij Z_i Z_j
void evolveZ(state *psi, ham *H, double t);

//Remove bruised zeros, i.e. things that are exactly zero but appear to
//be 10^-15 or something due to numerical noise.
void debruise(state *psi);

//Allocate memory for the amplitudes and set N = 2^n
//Returns 1 on success, 0 on failure
int allocate(state *psi, int n);

//Deallocate the memory for the amplitudes
void deallocate(state *psi);

//Allocate memory for the coefficients and set N = 2^n
//Returns 1 on success, 0 on failure
int allocateH(ham *H, int n);

//Deallocate the memory for the amplitudes
void deallocateH(ham *psi);

//Initialize psi to the uniform superposition
void uniform(state *psi);

//Initialize psi to the all zeros bit string
void zeros(state *psi);

//Initialize Hamiltonian to a long-range Ising with global X coupling
void lr(ham *H, double zc, double xc, double alpha);

//Return the expectation value of the Hamiltonian in state psi. The Hamiltonian is
//Xcoeff * \sum_i X_i + Zcoeff * \sum_{j>i} Z_i Z_j / |j-i|^alpha.
double expectH(state psi, ham *H);

// Returns the analytically computed energy for p=1 QAOA on angles beta, gamma
double qaoa1energy(ham *H, double beta, double gamma);

// Returns the analytically computed energy for p=1 QAOA on angles beta, gamma
// but with a different sum-of-sigma-X eigenstate
double signedqaoa1energy(ham *H, double beta, double gamma, double *s);

// Entanglement entropy of a pure state, with a cut after the kth qubit
double entent(state * psi, int k);

//The expecation of the global spin-z operator squared
double Jz2(state *psi);

//The expecation of the global spin-z operator squared
double Jz(state *psi);

//The expecation of the global spin-z operator squared
double JzVar(state *psi);
