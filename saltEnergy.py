from math import sqrt, log, log10, exp, gamma, pi, sinh
from scipy.special import exp1, kn
# Note that exp1(x) in scipy agrees with the incomplete Gamma function Gamma(0,x)
import numpy as np

import RNA


# Constants
# Length in angstrom (A)
# Length of one ssRNA backbone bond, noted as l_{ss}
BackboneLen = 6.0
# BackboneLen = 6.4
# Bjerrum length
BjerrumLen = 7
# RNA Helical rise per base pair (in A, noted as l_{ds})
# HelicalRise = 3.4
HelicalRise = 2.8
# Kuhn length
KuhnLen = 19
RodsDist = 20

# Elementary charge (in C)
e0 = 1.602e-19
# Boltzmann constant (in CVK^-1)
kB = 8.617e-5*e0
# in joule (i.e. m^2kgs^-2)
kBJoule = 1.380649e-23
# in kcal mol^-1K^-1
kBcal = 1.987e-3
# Vacuum dielectric constant (in CV^-1A^-1)
epsilon0 = 8.854e-22
# Relative dielectric constant of water
epsilonr = 80
# Avogadro constant (in mol^-1)
NA = 6.022e23
# Eular constant
eular = 0.58

evToJoule = 1.602e-19
jouleTokCal = 2.39e-4

SALT = 1.021

# Default temperature in RNAlib
TEMP = RNA.param().temperature+RNA.K0

# Store linearization for multiloop energy given salt concentration and temp
LinearTable = {}

LoopSalt = {}

PairingSalt = {}

def epsilonr(T):
    return 5321/T + 233.76 - 0.9297*T + 1.417*T*T/1000 - 0.8292*T*T*T/1000000;


def bjerrumlength(T=TEMP, eps=True):
    """Returns the Bjerrum length, the value should roughly be 7 (in A)
    """
    # denom = kB*T*4*pi*epsilon0*epsilonr
    # return e0*e0/denom

    # Precomputation
    if not eps:
        return 2089.6 / T

    return 167092.53/(T*epsilonr(T));

def ionicStrength(rho, Za=1, Zc=1):
    """Returns ionic strength given the salt concentration

    I = 1/2 * (rho_a*Z_a^2 + rho_c*Z_c^2)
    """
    return 1/2 * (rho*Za*Za + rho*Zc*Zc)


# Debye screening length
def kappaInv(rho, T=TEMP):
    """Returns Debye screening length given the salt concentration
    """
    # Exact calculation
    # Change the unit from dm^-3 to A^-3
    # rho = rho * 1e-27
    # nom = epsilon0*epsilonr*kB*T
    # denom = 2*NA*e0*e0*ionicStrength(rho)
    # return sqrt(nom/denom)

    # After precomputation
    # The adjustment of ionic strength unit is taken into account
    return 8.1285/sqrt(bjerrumlength(T)*ionicStrength(rho))



def approxHyper(y):
    a = 1/((y**6)/((2*pi)**6)+1)
    b = y**4/(36*pi**4) - y**3/(24*pi**2) + y*y/(2*pi*pi) - y/2
    c = log(2*pi/y) - 1.96351
    return a*b + (1-a)*c


def tauss(T=TEMP, Zc=1, l=None):
    """Returns line charge density
    """
    if l is None:
        return min(1/BackboneLen, 1/(bjerrumlength(T)*Zc))
    else:
        return softmin(l)(1/BackboneLen, 1/(bjerrumlength(T)*Zc))


def tauds(T=TEMP, Zc=1):
    """Returns line charge density for helix
    """
    return min(1/HelicalRise, 1/(bjerrumlength(T)*Zc))


def softmin(t):
    return lambda p, q: (exp(-t*p)*p+exp(-t*q)*q)/(exp(-t*p)+exp(-t*q))

def loopSaltAux(kmlss, m, T, l=None):
    a = kBcal*T*bjerrumlength(T)*m*BackboneLen*tauss(T, l=l)*tauss(T, l=l)
    b = log(kmlss) - log(pi/2) + eular + approxHyper(kmlss) + 1/kmlss * (1-exp(-kmlss)+kmlss*exp1(kmlss))
    return a*b*100

def loopSaltEnergy(m, rho, T=TEMP, l=None):
    """Returns salt-adjusted free energy of a loop given size m and concentration rho (in M=mol dm^-3)
    The unit of free energy is 10 Cal mol^-1

    """

    before = loopSaltAux(1/kappaInv(rho, T)*m*BackboneLen, m, T, l=l)
    return  before - loopSaltAux(1/kappaInv(SALT,T)*m*BackboneLen, m, T, l=l)


def linearMultiloopSaltEnergy(rho, T=TEMP, interval=range(3,16)):
    """Returns fitted linear function for multillop salt-adjusted free energy
    numpy.linalg.lstsq is used
    """
    x = np.array([i for i in interval])
    y = np.array([loopSaltEnergy(m, rho, T) for m in x])
    A = np.vstack([x, np.ones(len(x))]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]
    return round(m), round(c)


def pairingSaltEnergy(rho, T=TEMP, origin=False):
    """Returns salt-adjusted free energy for the length of helix increasing by 1
    The unit of free energy is 10 Cal mol^-1
    If origin is True, an approximation offered in the origin paper for Bessel function is used, otherwise, scipy.special.kn is used
    """
    # Constants for Bessel function approximation
    # Note: At salt concentration 1 M, the value of approximation is not null

    b1 = 0.0315171
    b2 = 35.0754
    b3 = 1.62292
    b4 = 10
    b5 = 42.6381
    kappa = 1/kappaInv(rho, T)
    a = kBcal*T*bjerrumlength(T)*HelicalRise*tauds(T)*tauds(T)
    if origin:
        return a*(b1 + (b2-b3*log(b4*kappa))/(1+b5*kappa)) * 100
    else:
        return 2*a*(kn(0, RodsDist*kappa)-kn(0, RodsDist/kappaInv(SALT,T)))*100





#####################
##                 ## 
##                 ##
##   Tan & Chen    ##
##                 ##
##                 ##
#####################

# phosphate-phosphate distance (backbone length)
Tan_Chen_d = 6.4

# Note: N is the chain length and the loop size is N – 2

def tan_chen_loop_energy_na(salt, N, x=17, T=TEMP):
    a1 = (0.02*N - 0.026) * log(salt) + 0.54*N + 0.78
    b1 = (-0.01/(N+1) + 0.006) * log(salt) - 7/(N+1)/(N+1) - 0.01
    c1 = 0.07 * log(salt) + 1.8
    d1 = 0.21 * log(salt) + 1.5

    logZloop = a1 * log(N-x/Tan_Chen_d + 1) + b1*(N-x/Tan_Chen_d + 1)**2 - b1
    logZcoil = c1 * N - d1
    return - RNA.GASCONST * T * (logZloop - logZcoil) / 10

def tan_chen_loop_correction_na(salt, N, x=17, T=TEMP):
    return tan_chen_loop_energy_na(salt, N, x, T) - tan_chen_loop_energy_na(SALT, N, x, T)


# doi: 10.1529/biophysj.106.100388
# N is base pair count in helix
def tan_chen_helix_correction_na(N, salt):
    a1 = -0.075 * log(salt) + 0.012 * log(salt) * log(salt)
    b1 = 0.018 * log(salt) * log(salt)
    return (N-1) * (a1 + b1/N) * 100



