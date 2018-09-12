/*---------------------------------------------------------------------------
  Title    : Simulator and data analysis tools for Quantum Approximate 
  Optimization Algorithm and Variational Quantum Eigensolvers.
  Authors  : Stephen Jordan (Microsoft QuARC) and Aniruddha Bapat (QuICS) 
  Year     : 2018
  Citation : If you use this code in your research, please cite it.
  -------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include "2local-qaoa.h"

//Print out the integer x in binary
void printbits(int x, int n) {
  int i;
  for(i = n-1; i >= 0; i--) printf("%i", (x>>i)&1);
}

//Print out the quantum state vector
void printvec(state psi) {
  int i;
  double p;
  for(i = 0; i < psi.N; i++) {
    p = psi.realcur[i]*psi.realcur[i] + psi.imagcur[i]*psi.imagcur[i];
    printf("%i(", i);
    printbits(i, psi.n);
    printf("):\t%e\t+i%e\t%e\n", psi.realcur[i], psi.imagcur[i], p);
  }
}

//Return the Euclidean norm of the state vector
double norm(state psi) {
  double val;
  int i;
  val = 0;
  for(i = 0; i < psi.N; i++) {
    val += psi.realcur[i]*psi.realcur[i];
    val += psi.imagcur[i]*psi.imagcur[i];
  }
  return val;
}

//Swap the current state with the buffered state
void swapbuf(state *psi) {
  double *tmp;
  tmp = psi->realcur;
  psi->realcur = psi->realbuf;
  psi->realbuf = tmp;
  tmp = psi->imagcur;
  psi->imagcur = psi->imagbuf;
  psi->imagbuf = tmp;
}  

//Perform a Hadamard transform on the n qubits in O(NlogN) time.
void Hadamard(state *psi) {
  int i,j;
  double root;
  root = 1.0/sqrt(2.0);
  for(i = 0; i < psi->n; i++) {
    //apply a Hadamard gate to the ith qubit and write result into buffer
    for(j = 0; j < psi->N; j++) {
      if((j&(1<<i)) == 0) {
        psi->realbuf[j] = root*(psi->realcur[j] + psi->realcur[j^(1<<i)]);
	psi->imagbuf[j] = root*(psi->imagcur[j] + psi->imagcur[j^(1<<i)]);
      }
      else {
        psi->realbuf[j] = root*(-psi->realcur[j] + psi->realcur[j^(1<<i)]);
	psi->imagbuf[j] = root*(-psi->imagcur[j] + psi->imagcur[j^(1<<i)]);
      }
    }
    swapbuf(psi); //swap cur with buf
  }
}

//Computes the Z2 inner product between bitstrings a and b of length n.
int Z2inner(int A, int B, int n) {
  int i;
  int telescope;
  int AandB;
  AandB = A&B;
  telescope = 0;
  for(i = 0; i < n; i++) telescope ^= (AandB>>i);
  return telescope&1;
}

//Convert an array of 1s and 0s to the corresponding number via binary place value.
//The zeroth element of the array contains the 1s place, the 1th element the 2s etc.
int b2i(int *bitstring, int n) {
  int i;
  int val;
  val = 0;
  for(i = 0; i < n; i++) val += bitstring[n-i-1]*(1<<i);
  return val;
}

//Convert a string of characters like "1100" to the corresponding number like 12.
//Note that in keeping with standard conventions, this is the reverse order for
//significant bits than in bin2int.
int cb2i(char *bitstring, int n) {
  int i;
  int val;
  val = 0;
  //the ascii code for '0' is 48, for '1' is 49, etc.
  for(i = 0; i < n; i++) val += ((int)bitstring[n-i-1]-48)*(1<<i);  
  return val;
}

//Returns the number of ones in the binary expansion of x
int HammingWeight(int x, int n) {
  int weight;
  int i;
  weight = 0;
  for(i = 0; i < n; i++) weight += (x>>i)&1;
  return weight;
}

//Returns b_i*x_i for x_i a bit in the binary expansion of x
double SkewWeight(int x, int n, double *b) {
  double weight;
  int i;
  weight = 0;
  for(i = 0; i < n; i++) weight += b[i]*((x>>i)&1);
  return weight;
}

//Unitarily evolve for time t according to H = - sum_j X_j
void evolveX(state *psi, double t) {
  int i;
  double angle;
  //Hadamard transform into the basis in which H is diagonal
  Hadamard(psi);
  //Do the time evolution in this basis
  for(i = 0; i < psi->N; i++) {
    angle = t*(double)(psi->n - 2*HammingWeight(i,psi->n));
    //multiply the amplitude by e^i(angle)
    psi->realbuf[i] = cos(angle)*psi->realcur[i] - sin(angle)*psi->imagcur[i];
    psi->imagbuf[i] = sin(angle)*psi->realcur[i] + cos(angle)*psi->imagcur[i];
  }
  swapbuf(psi);  //swap cur with buf
  Hadamard(psi); //Hadamard transform back to the computational basis
}

//The energy associated with a bit string according to Hamiltonian
//x stores the bit string
double energy(int x, ham *H) {
  int i;
  int j;
  double val;
  int biti, bitj;
  int n;
  double Ksum=0;
  val = 0;
  n   = H->n;
  for(i = 0; i < n; i++) {
    biti = (x>>i)&1;
    Ksum+=H->zc[i]; 
    for(j = i+1; j < n; j++) {
      bitj = (x>>(j))&1;
      val += H->zzc[n*i+j]*(1-((biti^bitj)<<1));
    }
  }
  val += Ksum - 2*SkewWeight(x, n, H->zc);
  return val;
}

//Unitarily evolve for time t according to H = \sum_i<j J_ij Z_i Z_j
void evolveZ(state *psi, ham *H, double t) {
  int i;
  double angle;
  for(i = 0; i < psi->N; i++) {
    angle = (-t)*energy(i, H);
    psi->realbuf[i] = cos(angle)*psi->realcur[i] - sin(angle)*psi->imagcur[i];
    psi->imagbuf[i] = sin(angle)*psi->realcur[i] + cos(angle)*psi->imagcur[i];
  }
  swapbuf(psi);  //swap cur with buf
}

//Remove bruised zeros, i.e. things that are exactly zero but appear to
//be 10^-15 or something due to numerical noise.
void debruise(state *psi) {
  int i;            //counter variable
  double threshold; //amplitudes smaller than this will be zeroed
  threshold = 1.0E-12;
  for(i = 0; i < psi->N; i++) {
    if(fabs(psi->realcur[i]) < threshold) psi->realcur[i] = 0;
    if(fabs(psi->imagcur[i]) < threshold) psi->imagcur[i] = 0;
  }
}

//Allocate memory for the amplitudes and set N = 2^n
//Returns 1 on success, 0 on failure
int allocate(state *psi, int n) {
  psi->n = n;
  psi->N = 1<<n;
  psi->realcur = (double *)malloc(psi->N*sizeof(double));
  if(psi->realcur == NULL) return 0;
  psi->imagcur = (double *)malloc(psi->N*sizeof(double));
  if(psi->imagcur == NULL) return 0;  
  psi->realbuf = (double *)malloc(psi->N*sizeof(double));
  if(psi->realbuf == NULL) return 0;  
  psi->imagbuf = (double *)malloc(psi->N*sizeof(double));
  if(psi->imagbuf == NULL) return 0;
  return 1;
}

//Deallocate the memory for the amplitudes
void deallocate(state *psi) {
  free(psi->realcur);
  free(psi->imagcur);
  free(psi->realbuf);
  free(psi->imagbuf);
}

int allocateH(ham *H, int n) {
  H->n = n;
  H->N = 1<<n;
  H->zzc = (double *)malloc(n*n*sizeof(double));
  if(H->zzc == NULL) return 0;
  H->xc = (double *)malloc(n*sizeof(double));
  if(H->xc == NULL) return 0;  
  H->zc = (double *)malloc(n*sizeof(double));
  if(H->zc == NULL) return 0;  
  return 1;
}

//Deallocate the memory for the amplitudes
void deallocateH(ham  *H) {
  free(H->zzc);
  free(H->xc);
  free(H->zc);
}

//Initialize psi to the uniform superposition
void uniform(state *psi) {
  int i;
  double amp;
  amp = 1.0/sqrt((double)psi->N);
  for(i = 0; i < psi->N; i++) {
    psi->realcur[i] = amp;
    psi->imagcur[i] = 0;
    psi->realbuf[i] = 0;
    psi->imagbuf[i] = 0;
  }
}  

//Initialize psi to the all zeros bit string
void zeros(state *psi) {
  int i;
  psi->realcur[0] = 1.0;
  psi->imagcur[0] = 0;
  psi->realbuf[0] = 0;
  psi->imagbuf[0] = 0;
  for(i = 1; i < psi->N; i++) {
    psi->realcur[i] = 0;
    psi->imagcur[i] = 0;
    psi->realbuf[i] = 0;
    psi->imagbuf[i] = 0;
  }
}

//Initialize Hamiltonian to a long-range Ising with global X coupling
void lr(ham *H, double Zcoeff, double Xcoeff, double alpha) {
  int i;
  int j;
  int dist2;
  int n = H->n;
  for(i=0; i < n; i++){
    for(j=0; j < n; j++){
      if(i==j) H->zzc[n*i+j] = 0;
      else{
	
	dist2 = (i-j)*(i-j);
	H->zzc[n*i+j] = Zcoeff/((pow((double)dist2,0.5*alpha)));
      }
    }
  }
  for(i=0; i < n; i++){
    H->zc[i]=(double)0.0;
    H->xc[i]=(double)Xcoeff;
  }
}

//Initialize Hamiltonian to a periodic long-range Ising with global X coupling
void lrPer(ham *H, double Zcoeff, double Xcoeff, double alpha) {
  int i;
  int j;
  int dist2, dist2per, distmin;
  int n = H->n;
  for(i=0; i < n; i++){
    for(j=0; j < n; j++){
      if(i==j) H->zzc[n*i+j] = 0;
      else{
	dist2 = abs(i-j);
	dist2per = n-abs(i-j);
	distmin = (dist2<dist2per)?dist2:dist2per;
	H->zzc[n*i+j] = Zcoeff/((pow((double)distmin,alpha)));
      }
    }
  }
  for(i=0; i < n; i++){
    H->zc[i]=(double)0.0;
    H->xc[i]=(double)Xcoeff;
  }
}

// Initialize to H=0
void empty(ham *H) {
  int i;
  int j;
  int n = H->n;
  for(i=0; i < n; i++){
    H->zc[i]=(double)0.0;
    H->xc[i]=(double)0.0;    
    for(j=0; j < n; j++){
      H->zzc[n*i+j] = (double)0.0;
    }
  }
}


//Return the expectation value of the Hamiltonian in state psi. The Hamiltonian is
//Xcoeff * \sum_i b_i X_i +  \sum_{i<j} J_ij Z_i Z_j
double expectH(state psi, ham * H) {
  int x;
  double p;
  double expectZ, expectX;
  state psitmp;
  int success;
  double Bsum;
  for(int i=0; i < H->n; i++){
    Bsum+=H->xc[i];
  }
  //We'll need to do a Hadamard transform, but we don't want to mess up
  //the state in case this is not the end of the protocol. So, we'll copy
  //the state to psitmp. (Note that we are passing psi directly, not by
  //reference, so we are assured not to mess it up.
  success = allocate(&psitmp, psi.n);
  if(!success) {
    printf("Memory allocation failure\n");
    return 0;
  }
  for(x = 0; x < psi.N; x++) {
    psitmp.realcur[x] = psi.realcur[x];
    psitmp.imagcur[x] = psi.imagcur[x];
  }
  expectZ = 0;
  for(x = 0; x < psitmp.N; x++) {
    p = psitmp.realcur[x]*psitmp.realcur[x] + psitmp.imagcur[x]*psitmp.imagcur[x];
    expectZ += p*energy(x,H);
  }
  Hadamard(&psitmp);
  expectX = 0;
  for(x = 0; x < psitmp.N; x++) {
    p = psitmp.realcur[x]*psitmp.realcur[x] + psitmp.imagcur[x]*psitmp.imagcur[x];
    //printf("p(x) = %f\n", p);
    expectX += p*(double)(Bsum - 2*SkewWeight(x,psitmp.n,H->xc));
  }
  deallocate(&psitmp);
  return expectX + expectZ;
}

//// Ani's functions

// Evolve the state by a p-iteration QAOA protocol
void evolveQAOA(state *psi, ham *H, double * beta, double * gamma, int p) {
  int i;
  for(i=0;i<p;i++) {
    evolveZ(psi, H, gamma[i]);
    evolveX(psi,    beta[i]);
  }    
}

// Compute |<psi1|psi2>|
double overlap(state *psi1, state *psi2) {
  int i;
  double real=0;
  double imag=0;
  double mag;
  if(psi1->N!=psi2->N) {
   printf("Error: dimensions must match.\n");
    return 0;
  }
  for(i=0;i<psi1->N;i++) {
    real+=(psi1->realcur[i]*psi2->realcur[i] + psi1->imagcur[i]*psi2->imagcur[i]);
    imag+=(psi1->realcur[i]*psi2->imagcur[i] - psi1->imagcur[i]*psi2->realcur[i]);
  }
  mag = sqrt(real*real + imag*imag);
  return mag;
}

// Returns the analytically computed energy for p=1 QAOA on angles beta, gamma
// but with a different sum-of-sigma-X eigenstate
double qaoa1energy(ham *H, double beta, double gamma){
  double perp=0;
  double par1=0;
  double par2=0;
  double par3=0;
  double prod1=1;
  double prod2=1;
  int n = H->n;
  int i,j,k;
  // X term contributions
  for(i=0; i<n; i++){
    prod1=1.0;
    for(j=0; j<n; j++){
      if(i==j) prod1*=1.0;
      else     prod1*=cos(2*gamma*H->zzc[n*i+j]);
    }
    perp += H->xc[i]*cos(2*gamma*H->zc[i])*prod1;
  }
  // One-local Z term contributions
  for(i=0;i<n;i++){
    prod1=1.0;
    for(j=0; j<n; j++){
      if(i==j) prod1*=1.0;
      else     prod1*=cos(2*gamma*H->zzc[n*i+j]);
    }
    par1 -= H->zc[i]*sin(2*gamma*H->zc[i])*sin(2*beta)*prod1;
  }
  // Two-local ZZ term contributions (Part 1)
  for(i=0; i<n; i++){
    for(j=i+1;j<n; j++){
      prod1=1.0;
      prod2=1.0;
      for(k=0;k<n;k++){
	if(k==i||k==j){
	  prod1*=1.0;
	  prod2*=1.0;
	}
	else{
	  prod1*=cos(2*gamma*H->zzc[n*i+k]);
	  prod2*=cos(2*gamma*H->zzc[n*j+k]);
	}
      }
      par2 += -0.5*sin(4*beta)*H->zzc[n*i+j]*sin(2*gamma*H->zzc[n*i+j])*\
	(cos(2*gamma*H->zc[i])*prod1 + cos(2*gamma*H->zc[j])*prod2);
    }
  }
  // Two-local ZZ term contributions (Part 2)
  for(i=0; i<n; i++){
    for(j=i+1;j<n; j++){
      prod1=1.0;
      prod2=1.0;
      for(k=0;k<n;k++){
	if(k==i||k==j){
	  prod1*=1.0;
	  prod2*=1.0;
	}
	else{
	  prod1*=cos(2*gamma*(H->zzc[n*i+k]+H->zzc[n*j+k]));
	  prod2*=cos(2*gamma*(H->zzc[n*i+k]-H->zzc[n*j+k]));
	}
      }
      par3 += 0.25*(cos(4*beta)-1)*H->zzc[n*i+j]*\
	(cos(2*gamma*(H->zc[i]+H->zc[j]))*prod1 - cos(2*gamma*(H->zc[i]-H->zc[j]))*prod2);
    }
  }
  return perp+par1+par2+par3;
}

// Returns the analytically computed energy for p=1 QAOA on angles beta, gamma
double signedqaoa1energy(ham *H, double beta, double gamma, double *s){
  double perp=0;
  double par1=0;
  double par2=0;
  double par3=0;
  double prod1=1;
  double prod2=1;
  int n = H->n;
  int i,j,k;
  // X term contributions
  for(i=0; i<n; i++){
    prod1=1.0;
    for(j=0; j<n; j++){
      if(i==j) prod1*=1.0;
      else     prod1*=cos(2*gamma*H->zzc[n*i+j]);
    }
    perp += H->xc[i]*cos(2*gamma*H->zc[i])*prod1;
  }
  // One-local Z term contributions
  for(i=0;i<n;i++){
    prod1=1.0;
    for(j=0; j<n; j++){
      if(i==j) prod1*=1.0;
      else     prod1*=cos(2*gamma*H->zzc[n*i+j]);
    }
    par1 -= s[i]*H->zc[i]*sin(2*gamma*H->zc[i])*sin(2*beta)*prod1;
  }
  // Two-local ZZ term contributions (Part 1)
  for(i=0; i<n; i++){
    for(j=i+1;j<n; j++){
      prod1=s[i];
      prod2=s[j];
      for(k=0;k<n;k++){
	if(k==i||k==j){
	  prod1*=1.0;
	  prod2*=1.0;
	}
	else{
	  prod1*=cos(2*gamma*H->zzc[n*i+k]);
	  prod2*=cos(2*gamma*H->zzc[n*j+k]);
	}
      }
      par2 += -0.5*sin(4*beta)*H->zzc[n*i+j]*sin(2*gamma*H->zzc[n*i+j])*\
	(cos(2*gamma*H->zc[i])*prod1 + cos(2*gamma*H->zc[j])*prod2);
    }
  }
  // Two-local ZZ term contributions (Part 2)
  for(i=0; i<n; i++){
    for(j=i+1;j<n; j++){
      prod1=1.0;
      prod2=1.0;
      for(k=0;k<n;k++){
	if(k==i||k==j){
	  prod1*=1.0;
	  prod2*=1.0;
	}
	else{
	  prod1*=cos(2*gamma*(H->zzc[n*i+k]+H->zzc[n*j+k]));
	  prod2*=cos(2*gamma*(H->zzc[n*i+k]-H->zzc[n*j+k]));
	}
      }
      par3 += 0.25*s[i]*s[j]*(cos(4*beta)-1)*H->zzc[n*i+j]*\
	(cos(2*gamma*(H->zc[i]+H->zc[j]))*prod1 - cos(2*gamma*(H->zc[i]-H->zc[j]))*prod2);
    }
  }
  return perp+par1+par2+par3;
}

// Entanglement entropy of a pure state, with a cut after the kth qubit
double entent(state * psi, int k){
  double ent=0;
  int i;
  
  
  return ent;
}

//The eigenvalue of the global spin-z operator squared
//x stores the bit string
double Jz2(state *psi){
  double sum = 0, sq;
  int i,n = psi->n;
  for(i=0;i<psi->N;i++){
    sq = HammingWeight(i,n) - 0.5*n;
    sum += sq*sq*(psi->realcur[i]*psi->realcur[i]+psi->imagcur[i]*psi->imagcur[i]);
  }
  return sum;
}

//The eigenvalue of the global spin-z operator squared
//x stores the bit string
double Jz(state *psi){
  double sum = 0;
  int i,n = psi->n;
  for(i=0;i<psi->N;i++){
    sum += (HammingWeight(i,n) - 0.5*n)*(psi->realcur[i]*psi->realcur[i]+psi->imagcur[i]*psi->imagcur[i]);
  }
  return sum;
}

double JzVar(state *psi){
  return Jz2(psi)-Jz(psi)*Jz(psi);
}
