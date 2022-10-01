#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:44:57 2022

@author: morteza
"""

import numpy as np
import matplotlib.pyplot as plt



vf = np.load("volume_fraction.npy")
vfII = np.load("volume_fractionII.npy")
acceptance_rate = np.load("acceptance_rate.npy")


plt.xlabel("time")
plt.ylabel("volume fraction")
plt.plot(np.arange(0,3000), vf[:3000])


#plt.xlabel("volume fraction")
#plt.ylabel("acceptance rate")
#plt.plot(vfII/1000, acceptance_rate)