import numpy as np

f = open("mask2.dat", "w")

for x in range(5000):
	for y in range(5000):
#		if (y > 150 and y <= 450 and ((x > 40 and x <= 190) or (x > 410 and x <= 560)) or (x > 190 and x <= 410 and y > 289 and y <= 311)):
		if (y > 250 and y <= 2250 and ((x > 500 and x <= 1400) or (x > 2600 and x <= 3500)) or (x > 1400 and x <= 2600 and y >= 1190 and y <= 1310)):
			f.write(str(x) + "," + str(y) + "\n")
			
f.close()
			
