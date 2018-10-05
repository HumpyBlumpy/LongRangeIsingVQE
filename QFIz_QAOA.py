#from FnLibrary import *
from qaoaLibrary import *
import sys
import time
import numpy as np

n=11
# Initialize state and hamiltonian
psi, H = initialize(n)
Zcoeff = 1
Xcoeff = 1 

alpha = 0.5
ppsi = pointer(psi)
qc.uniform(ppsi)

# Hamiltonian Input
# OPTION 1: Generate long-range power-law Hamiltonian
qc.lr(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#qc.lrPer(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
H.zc = (c_double * n)(*([-0.00]*n))


stuff = JzVarOptQAOA(psi, H,3,'BFGS')
print(stuff)
print(qc.Jz(ppsi))
print(qc.Jz2(ppsi))
print('F_Q(J_z)/N='+str(4*qc.JzVar(ppsi)/n))
print(qc.Jx(ppsi))
print(qc.Jx2(ppsi))
print(qc.JxVar(ppsi))
