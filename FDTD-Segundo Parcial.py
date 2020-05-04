import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import constants

Npts = 200
Ez = np.zeros(Npts)
Hy = np.zeros(Npts)
E0 = 400

eps_0 = constants.epsilon_0
mu_0 = constants.mu_0

eps = np.zeros(Npts)

eps[0:int(Npts/2)] = 10*eps_0
eps[int(Npts/2):Npts] = eps_0

mu = np.zeros(Npts)

mu[0:int(Npts/2)] = mu_0/10
mu[int(Npts/2):Npts] = mu_0

vel = constants.c
dt = 0.18e-9
dx = 2*dt*vel

factor_norm = math.sqrt(eps_0/mu_0)

pc = int(Npts/2)
t0 = 40
spread = 12

Npasos = 450