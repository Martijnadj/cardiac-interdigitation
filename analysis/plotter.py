import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy import signal
from numpy import genfromtxt
import pandas as pd
from numpy.fft import fft, ifft
import os
import fnmatch
cbmap=['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
                  
                  
Q_tot = []
Q_thr = []

dt = 1e-5
#dt = 1e-6

Q_tot = genfromtxt("Q_tot.txt", delimiter=',')
Q_thr = genfromtxt("Q_thr.txt", delimiter=',')
time = genfromtxt("time_file.txt", delimiter = ',')


dt = 1e-5*20

plt.scatter(time, Q_tot, c = cbmap[0], label = "Q_tot")
plt.scatter(time, Q_thr, c = cbmap[1], label = "Q_thr")
plt.legend()
plt.show()

SF = Q_tot/Q_thr

plt.scatter(time, SF, c = cbmap[1], label = "SF")

plt.legend()
plt.show()

