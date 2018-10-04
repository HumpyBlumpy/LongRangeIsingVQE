#from FnLibrary import *
from qaoaLibrary import *
import sys
import time

n=10
# Initialize state and hamiltonian
psi, H = initialize(n)

# Hamiltonian Input
# OPTION 1: Generate long-range power-law Hamiltonian
#qc.lr(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#qc.lrPer(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#H.zc = (c_double * n)(*([-0.00]*n))

ppsi = pointer(psi)
qc.uniform(ppsi)
print(qc.Jz(ppsi))
print(qc.Jz2(ppsi))
print(qc.JzVar(ppsi))

