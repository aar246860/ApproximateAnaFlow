import scipy.special as sc
import numpy as np
from scipy.special import expi
from mpmath import *
mp.dps = 10; mp.prec = 60; mp.pretty = True
import matplotlib.pyplot as plt
import warnings
def app_ext_theis_2d(time, rad, S, TG, var, Len_scale, rate):
    zeta = 1.6
    rw = 0.05
    R = rw + 1.49861*np.sqrt(TG*time/S) 
    lamda1 = -var*(rad*zeta)**2/(2*(Len_scale**2+(rad*zeta)**2))
    lamda2 = var*Len_scale**2/(2*(Len_scale**2+(rad*zeta)**2))
    lamda1R = -var*(R*zeta)**2/(2*(Len_scale**2+(R*zeta)**2))
    lamda2R = var*Len_scale**2/(2*(Len_scale**2+(R*zeta)**2))
    Gamma = np.exp(var/2)*expi(lamda1)-expi(lamda2)
    GammaR = np.exp(var/2)*expi(lamda1R)-expi(lamda2R)
    A1 = (rate/(4*np.pi*TG))*GammaR
    A2 = rate/(2*np.pi*TG)
    sp = (rate/(4*np.pi*TG))*(GammaR-Gamma)
    return sp
