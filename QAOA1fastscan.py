from qaoaLibrary import *
import sys
import time
from matplotlib.pyplot import *

# Perform a fast optimization of QAOA1 angles.
# The code returns optimal beta, gamma and energy
def optQAOA1(H, typeOfOpt='Nelder-mead'):
    pH   = pointer(H)
    fOpt = lambda bg: expectQAOA1(H, bg.tolist())
    bg0  = 0.5*np.ones(2)
    opt    = minimize(fOpt, bg0, method=typeOfOpt)
    bgOpt  = opt.x
    E      = expectQAOA1(H, opt.x.tolist())
    return (bgOpt, E) 


for n in range(10, 100, 10):
    # Hamiltonian
    H = ham()
    qc.allocateH(byref(H),n)

    # OPTION 1: Generate exact long-range Hamiltonian (or any other Hamiltonian of choice)
    J = 1.0
    h = -0.3
    alpha = 0.8

    # Replace qc.lr by qc.lrPer for Periodic boundary conditions
    qc.lr(byref(H),c_double(J), c_double(h), c_double(alpha))
    # Set all single Z terms to 0, or other value of choosing 
    H.zc = (c_double * n)(*([-0.00]*n))

    # OPTION 2: Read Hamiltonian coefficients from file
    #Jij = np.loadtxt('expJ-n%d.dat'%n)
    #hamGen(H, Jij, Xcoeff*np.ones(n), np.zeros(n))
    #Zcoeff = np.mean([Jij[i,i+1] for i in range(n-1)]) 

    # Find the optimal beta, gamma values. This is optional,
    # but if you want to fix one angle at its optimum value, it
    # will help to have these values pre-computed.
    bgOpt, Eopt = optQAOA1(H)
    
    # Now, sweep angles beta, gamma in QAOA1 and compute output energy
    # E(beta, gamma) using the efficient analytical formula.
    # This is fast! You can probably output energy curves for n=100
    # with reasonable computation times
    E = []
    for i in range(100):
        E.append(-expectQAOA1(H,[abs(bgOpt[0]), 0.01*i])/Eopt)
    plot(0.01*np.arange(100), E, label='n=%d'%n)
    ylim(-1.05,0.3)
    ylabel("Energy (normalized so that min. E = -1)")
    xlabel("gamma")
    legend(loc='lower right')
savefig('QAOA1_gamma_n_scan.png')
show()

