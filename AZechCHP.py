import scipy.special as sc
import numpy as np
from scipy.special import expi
from mpmath import *
mp.dps = 10; mp.prec = 60; mp.pretty = True
import warnings

#for constant-head pumpint#
def app_ext_CHP_2d(time, rad, S, TG, var, Len_scale, s0):
    zeta = 1.6
    rw = 0.05
    R = rw + np.sqrt(np.pi*TG*time/S) 
    lamda1 = -var*(rad*zeta)**2/(2*(Len_scale**2+(rad*zeta)**2))
    lamda2 = var*Len_scale**2/(2*(Len_scale**2+(rad*zeta)**2))
    lamda1R = -var*(R*zeta)**2/(2*(Len_scale**2+(R*zeta)**2))
    lamda2R = var*Len_scale**2/(2*(Len_scale**2+(R*zeta)**2))
    lamda1rw = -var*(rw*zeta)**2/(2*(Len_scale**2+(rw*zeta)**2))
    lamda2rw = var*Len_scale**2/(2*(Len_scale**2+(rw*zeta)**2))
    Gamma = np.exp(var/2)*expi(lamda1)-expi(lamda2)
    GammaR = np.exp(var/2)*expi(lamda1R)-expi(lamda2R)
    Gammarw = np.exp(var/2)*expi(lamda1rw)-expi(lamda2rw)
    sp = s0*(GammaR-Gamma)/(GammaR-Gammarw)
    return sp

#for wellbore flowrate#
def Flow(time, S, TG, var, Len_scale, s0):
    zeta = 1.6
    rw = 0.05
    R = rw + np.sqrt(np.pi*TG*time/S) 
    lamda1 = -var*(rw*zeta)**2/(2*(Len_scale**2+(rw*zeta)**2))
    lamda2 = var*Len_scale**2/(2*(Len_scale**2+(rw*zeta)**2))
    lamda1R = -var*(R*zeta)**2/(2*(Len_scale**2+(R*zeta)**2))
    lamda2R = var*Len_scale**2/(2*(Len_scale**2+(R*zeta)**2))
    lamda1rw = -var*(rw*zeta)**2/(2*(Len_scale**2+(rw*zeta)**2))
    lamda2rw = var*Len_scale**2/(2*(Len_scale**2+(rw*zeta)**2))
    Gamma = np.exp(var/2)*expi(lamda1)-expi(lamda2)
    GammaR = np.exp(var/2)*expi(lamda1R)-expi(lamda2R)
    Gammarw = np.exp(var/2)*expi(lamda1rw)-expi(lamda2rw)
    sp = s0*4*np.pi*TG/(GammaR-Gammarw)
    return sp
