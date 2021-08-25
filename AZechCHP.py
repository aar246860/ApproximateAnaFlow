#!/usr/bin/env python
# coding: utf-8

# In[6]:


from anaflow import theis, ext_theis_2d
import scipy.special as sc
import numpy as np
from scipy.special import expi
from mpmath import *
mp.dps = 10; mp.prec = 60; mp.pretty = True
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import differential_evolution
import warnings


# In[24]:


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


# In[28]:


t1=np.logspace(0,4,num=100,base=10)
#T value
label_ef = "extended Theis"
label_af = "apprimate extended Theis"
TT1=[];
TT1=[];
TT2=[];
TT3=[];
TT4=[];
TT5=[];
for n in range(0, len(t1)):
    TT1.append(app_ext_CHP_2d(t1[n], 0.05, 0.0001, 0.0001, 1, 1, 1))
    TT2.append(app_ext_CHP_2d(t1[n], 0.05, 0.0001, 0.0001, 1, 2, 1))
    TT3.append(app_ext_CHP_2d(t1[n], 0.05, 0.0001, 0.0001, 1, 3, 1))
    TT4.append(app_ext_CHP_2d(t1[n], 0.05, 0.0001, 0.0001, 1, 4, 1))
    TT5.append(app_ext_CHP_2d(t1[n], 0.05, 0.0001, 0.0001, 1, 5, 1))
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
plt.plot(t1, TT4, label=label_ef, color="C"+str(4))
plt.plot(t1, TT5, label=label_ef, color="C"+str(5))


# In[18]:


t1=np.logspace(0,4,num=100,base=10)
#T value
label_ef = "extended Theis"
label_af = "apprimate extended Theis"
TT1=[];
TT1=[];
TT2=[];
TT3=[];
TT4=[];
TT5=[];
for n in range(0, len(t1)):
    TT1.append(Flow(t1[n], 0.0001, 0.0001, 1, 1, 1))
    TT2.append(Flow(t1[n], 0.0001, 0.0001, 1, 2, 1))
    TT3.append(Flow(t1[n], 0.0001, 0.0001, 1, 3, 1))
    TT4.append(Flow(t1[n], 0.0001, 0.0001, 1, 4, 1))
    TT5.append(Flow(t1[n], 0.0001, 0.0001, 1, 5, 1))
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
plt.plot(t1, TT4, label=label_ef, color="C"+str(4))
plt.plot(t1, TT5, label=label_ef, color="C"+str(5))


# In[ ]:




