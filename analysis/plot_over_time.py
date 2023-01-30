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
                  
                  


file_list = fnmatch.filter(os.listdir('.'), '*.txt')
file_list.sort()
print(file_list)
number_of_files = len(file_list)

files_to_plot = [number_of_files-1,0]
cutoff_threshold = 0 #difference from initial value from where we consider fft

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



t = data[0,:,0]
				
#cut off first data points to improve frequency result from fft

cutoff_loc = np.zeros(number_of_files)
for file in range(0, number_of_files):
	for i in range (0,measurements):
		if abs(data[file,0,1]-data[file,i,1]) > cutoff_threshold:
			cutoff_loc[file] = i-1
			break
		

FT = [None] * number_of_files
N = measurements - cutoff_loc
T = N*(t[1]-t[0])
for i in range (0,number_of_files):
	FT[i] =  fft(data[i,int(cutoff_loc[i]):,1])

def parabola(x, a, b, c):
	return a*x**2 + b*x + c
	

plt.figure(figsize = (12,6))
plt.subplot(121)
counter = 0
precision = 1e-2
for file in files_to_plot:
	plt.stem(np.arange(N[file])/T[file], np.abs(FT[file]), cbmap[counter], label = 'Measure location ' + str(file), markerfmt=" ")
	counter += 1
f = open("Frequencies.csv", "w")
for file in range(0, number_of_files):
	auto_cor = sm.tsa.acf(np.abs(FT[file]), nlags = 100)
	running_avg = []
	avg_points = 5 #use an odd number
	for x in range(len(auto_cor)-avg_points):
		running_avg.append(np.sum(auto_cor[x:x+avg_points]/avg_points))
	max_loc = np.where(running_avg == np.amax(running_avg[1:]))[0][0]
	print(max_loc)
	if max_loc == 1:
		loc = 0
	elif max_loc == 2: 
		fit_params, pcov = scipy.optimize.curve_fit(parabola, np.arange(int(avg_points/2)*0.05, (len(running_avg)+int(avg_points/2))*0.05-0.001, 0.05)[max_loc-1:max_loc+2], running_avg[max_loc-1:max_loc+2])
		print(fit_params)
		loc = -fit_params[1]/(2*fit_params[0])
	else: 
		fit_params, pcov = scipy.optimize.curve_fit(parabola, np.arange(int(avg_points/2)*0.05, (len(running_avg)+int(avg_points/2))*0.05-0.001, 0.05)[max_loc-2:max_loc+3], running_avg[max_loc-2:max_loc+3])
		print(fit_params)
		loc = -fit_params[1]/(2*fit_params[0])
	f.write(str(loc) + '\n')
	plt.xlim(0,5)
plt.ylim(0, 300000)
plt.xlabel('Freq (Hz)')
plt.ylabel('FFT Amplitude |Stimulus(freq)|')


plt.subplot(122)
counter = 0
for file in files_to_plot:
	print(file)
	plt.plot(t, data[file,:,1], cbmap[counter])
	plt.plot(t, data[file,:,1], cbmap[counter])
	counter += 1
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.tight_layout()
#plt.show()
plt.savefig("Voltage_and_FFT")
