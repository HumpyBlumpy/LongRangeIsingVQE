#from FnLibrary import *
from qaoaLibrary import *
import sys
import time

n=10
# Initialize state and hamiltonian
psi, H = initialize(n)
Zcoeff = 1
Xcoeff = 1 
alpha = 1
ppsi = pointer(psi)
qc.uniform(ppsi)

# Hamiltonian Input
# OPTION 1: Generate long-range power-law Hamiltonian
qc.lr(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#qc.lrPer(byref(H), c_double(Zcoeff), c_double(Xcoeff), c_double(alpha))
#H.zc = (c_double * n)(*([-0.00]*n))

pH   = pointer(H)

qc.evolveZ(ppsi, pH, 5)
qc.evolveX(ppsi,    6)
qc.evolveZ(ppsi, pH,8)
qc.evolveX(ppsi,    2)


print(qc.Jz(ppsi))
print(qc.Jz2(ppsi))
print(qc.JzVar(ppsi))
print(qc.Jx(ppsi))
print(qc.Jx2(ppsi))
print(qc.JxVar(ppsi))
