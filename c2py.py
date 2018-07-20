#from FnLibrary import *
from qaoaLibrary import *
import sys
import time

start_time = time.time()

# Read in parameters
n       = int(sys.argv[1])
N       = 1<<n
alpha   = float(sys.argv[2])
Xcoeff  = float(sys.argv[3])
Zcoeff  = float(sys.argv[4])
p       = int(sys.argv[5])
Tmax    = 10.0

# Initialize state and hamiltonian
psi, H = initialize(n)

# Hamiltonian Input
# OPTION 1: Generate long-range power-law Hamiltonian
qc.lr(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#qc.lrPer(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
H.zc = (c_double * n)(*([-0.00]*n))

# OPTION 2: Read in actual Hamiltonian parameters
#Jij = np.loadtxt('expJ-n%d.dat'%n)
#hamGen(H, Jij, Xcoeff*np.ones(n), np.zeros(n))
#Zcoeff = np.mean([Jij[i,i+1] for i in range(n-1)]) 

Neigs = 2
val, vec = ground(H, Neigs)

######################################################################
# This part focuses on optimizing and comparing with the true ground state of the Hamiltonian

bgOpt, Eopt = optQAOA(psi,H, p, 'BFGS')
printOut(psi, val, vec, bgOpt, Eopt)

trials = 1
EFull  = np.inf
bgFull  = []
#for i in range(trials):
#     bgCur, ECur= optQAOAgreedy(psi, H, p)#100, 2.0, 2000)
#     bgFull = bgFull if ECur>EFull else bgCur
#     EFull = EFull if ECur>EFull else ECur

#printOut(psi, val, vec, bgFull, EFull)

##### ENTANGLEMENT ENTROPY ##########

# True Ground state
#symGS = z2symmetrize(vec[:,0])

#print entent(symGS,1)

#print(entent(symGS,n/2))

# Print execution time
#print time.time() - start_time
