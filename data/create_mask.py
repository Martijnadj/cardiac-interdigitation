import numpy as np

f = open("mask.dat", "w")

for x in range(500):
	for y in range(500):
#		if (y > 150 and y <= 450 and ((x > 40 and x <= 190) or (x > 410 and x <= 560)) or (x > 190 and x <= 410 and y > 289 and y <= 311)):
		if (y > 175 and y <= 325 and ((x > 120 and x <= 195) or (x > 305 and x <= 380)) or (x > 195 and x <= 305 and y >= 245 and y <= 255)):
			f.write(str(x) + "," + str(y) + "\n")
			
f.close()
			
