#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 16:03:17 2020

@author: athani_m

Generates the historgram of the curvature data
"""

import os
import numpy as np
import matplotlib.pyplot as plt


os.getcwd()
Data = np.loadtxt(fname = "Curvature.txt")

plt.hist(Data[:,1], bins = 100, log = True)
plt.title('Curvature distribution')
plt.xlabel('Curvatures')
plt.ylabel('Counts')
plt.savefig('Histogram.png', dpi=300)