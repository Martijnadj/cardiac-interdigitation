import matplotlib.pyplot as plt
import csv
import numpy as np
from scipy import signal
import scipy.optimize
import statsmodels.api as sm
from numpy import genfromtxt
import pandas as pd
from numpy.fft import fft, ifft
import os
import fnmatch
cbmap=['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']
                  
                  


file_list = fnmatch.filter(os.listdir('.'), 'location_*')
file_list.sort()
number_of_files = len(file_list)

files_to_plot = [number_of_files-1,0]

# count number of measurements
results = genfromtxt(file_list[0], delimiter=',')
measurements = len(results)

#find dimensions for the data array
num_rows, num_cols = np.shape(results)

#create data array
data = np.empty([number_of_files, num_rows, num_cols])

#add data to array
file_number = 0
for file in file_list:
    data[file_number] = genfromtxt(file, delimiter=',')
    file_number += 1

number_of_measurements = 10
sizex = 610
measure_loc = [None]*number_of_measurements
for measurement in range(1,number_of_measurements+1):
    measure_loc[measurement-1] = int ((0.5 + int(10 * (number_of_measurements-measurement)/(number_of_measurements-1) + (sizex-11)*(measurement-1)/(number_of_measurements-1))));  


first_0_crossing = [None]*number_of_measurements
for i in range(10):
    first_0_crossing[i]=(np.argmax(data[i,:,1]>0))
    

velocity = [None]*(number_of_measurements-1)
dx = 3e-5
print(measure_loc[9])
print(measure_loc[5])
print(first_0_crossing)
for i in range(9):
    velocity[i] = 100*dx*(measure_loc[i+1] - measure_loc[i])/(data[0,first_0_crossing[i+1],0] - data[0,first_0_crossing[i],0])

for i in range(number_of_measurements-1):
    plt.plot((measure_loc[i],measure_loc[i+1]), (velocity[i],velocity[i]), 'b', alpha = 0.5)
plt.xlabel('Pixel')
plt.ylabel('Velocity (cm/s)')

print('Velocity in atrium = ' + str(100*dx*(measure_loc[9] - measure_loc[6])/(data[0,first_0_crossing[9],0] - data[0,first_0_crossing[6],0])))







plt.show()
