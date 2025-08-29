import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh


# To see how the ground state energy evolves with U/J

def curve_ground_state_energy(N,M,alpha):
    g = []
    Hj = BoseHubbard_Hj(N,M)
    Hu = BoseHubbard_Hu(N,M)
    start = time.time()
    for U in alpha:
        H = -Hj + U*Hu
        Eg = eigsh(H,1,which='SA')[0]
        Eg = np.real_if_close(Eg)
        g.append(Eg)
    return g

def plot_ground_state_energy(N,M):
    alpha = np.linspace(0,100,100)
    g = curve_ground_energy(N,M,alpha)
    plt.plot(alpha,g)
    plt.title("Energy of the ground state for N = "+str(N)+", M = "+str(M))
    plt.xlabel(r"$\frac{U}{J}$")
    plt.ylabel(r"$\frac{E_\text{gs}}{J}$")
    plt.show()

# To have a look at the one-body reduced density matrix

def plot_OBDM(N,M,U,J):
    OBDM = OBDM_ground(N,M,U,J)
    plt.imshow(OBDM,cmap="plasma",vmin=0, vmax=N/M)
    plt.colorbar()
    plt.title(f"N = {N}, M = {M}, U = {U} J")
    plt.show()

# To see the evolution of the condensed fraction with U/J

def curve_condensed_fraction(N,M,alpha):
    condensed_fraction = []
    for U in alpha:
        rho1 = OBDM_ground(N,M,U,1)
        c = eigsh(rho1,1,which='LA')[0]/np.trace(rho1)
        condensed_fraction.append(c)
    return condensed_fraction

def plot_condensed_fraction(N,M):
    alpha = np.linspace(0,100,100)
    cf = curve_condensed_fraction(N,M,alpha)
    plt.plot(alpha,cf)
    plt.title("Condensed fraction for N = "+str(N)+", M = "+str(M))
    plt.xlabel(r"$\frac{U}{J}$")
    plt.ylabel(r"$f_\text{cond}$")
    plt.show()















