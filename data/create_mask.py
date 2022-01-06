import numpy as np

f = open("mask.dat", "w")

for x in range(500):
	for y in range(500):
		if (y > 10 and y <= 160 and ((x > 10 and x <= 85) or (x > 195 and x <= 270)) or (x > 85 and x <= 195 and y > 79 and y <= 90)):
#		if (y > 250 and y <= 2250 and ((x > 500 and x <= 1400) or (x > 2600 and x <= 3500)) or (x > 1400 and x <= 2600 and y >= 1190 and y <= 1310)):
			f.write(str(x) + "," + str(y) + "\n")
			
f.close()
