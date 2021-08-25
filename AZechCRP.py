#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[7]:


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


# In[8]:


#Example1: Temperal drawdown curves predicted by both solutions with varying Len_scale from 1 to 5 observed at r = 1 m#
#r value
t1=np.logspace(0,4,num=100,base=10)
#T value
label_ef = "extended Theis"
label_af = "apprimate extended Theis"
TT1=[];
TT2=[];
TT3=[];
TT4=[];
TT5=[];
TTz1=[];
TTz2=[];
TTz3=[];
TTz4=[];
TTz5=[];
for n in range(0, len(t1)):
    TTz1.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TTz2.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 2, 0.0001))
    TTz3.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 3, 0.0001))
    TTz4.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 4, 0.0001))
    TTz5.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 5, 0.0001))
    TT1.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TT2.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 2, 0.0001))
    TT3.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 3, 0.0001))
    TT4.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 4, 0.0001))
    TT5.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 5, 0.0001))
#in semi-log scale
plt.plot(t1, TTz1, label=label_af, color="C"+str(1), marker=".",markevery=3)
plt.plot(t1, TTz2, label=label_af, color="C"+str(2), marker=".",markevery=3)
plt.plot(t1, TTz3, label=label_af, color="C"+str(3), marker=".",markevery=3)
plt.plot(t1, TTz4, label=label_af, color="C"+str(4), marker=".",markevery=3)
plt.plot(t1, TTz5, label=label_af, color="C"+str(5), marker=".",markevery=3)
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
plt.plot(t1, TT4, label=label_ef, color="C"+str(4))
plt.plot(t1, TT5, label=label_ef, color="C"+str(5))
plt.xlabel("t (sec)")
plt.ylabel("s (m)")
plt.show()


# In[142]:


#Example2: Temperal drawdown curves predicted by both solutions with varying var from 1 to 5 observed at r = 1#
#r value
t1=np.logspace(0,4,num=100,base=10)
var_labels = ["1", "3", "5"]
#T value
TT1=[];
TT2=[];
TT3=[];
TTz1=[];
TTz2=[];
TTz3=[];
for n in range(0, len(t1)):
    TTz1.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TTz2.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 3, 1, 0.0001))
    TTz3.append(ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 5, 1, 0.0001))
    TT1.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 1, 1, 0.0001))
    TT2.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 3, 1, 0.0001))
    TT3.append(app_ext_theis_2d(t1[n], 1, 0.0001, 0.0001, 5, 1, 0.0001))
#in semi-log scale
label_ef = "extended Theis"
label_af = "apprimate extended Theis"

plt.plot(t1, TTz1, label=label_af, marker=".",markevery=2)
plt.plot(t1, TTz2, label=label_af, marker=".",markevery=2)
plt.plot(t1, TTz3, label=label_af, marker=".",markevery=2)
plt.plot(t1, TT1, label=label_ef, color="C"+str(1))
plt.plot(t1, TT2, label=label_ef, color="C"+str(2))
plt.plot(t1, TT3, label=label_ef, color="C"+str(3))
ax2.set_yticklabels(time_labels)
plt.legend()
plt.xlabel("t (sec)")
plt.ylabel("s (m)")
plt.show()


# In[138]:


plt.plot(t1, TTz1, label=label_TG, marker=".",markevery=2)
plt.plot(t1, TTz2, label=label_TG, marker=".",markevery=2)
plt.plot(t1, TTz3, label=label_TG, marker=".",markevery=2)


# In[ ]:




