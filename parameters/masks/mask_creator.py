import numpy as np
import matplotlib.pyplot as plt
import os
import math
from os import path

pixel_size = 0.0025
#in millimeters

#paramaters for mask, see sketch.jpg for clarification

#Left shape
L_shape = "circle"
#Choose either "circle" or "rectangle"
L_radius = 200
L_width = 1200
L_height = 3200

#Isthmus
I_length = 600
I_width = 120

#Right shape
R_angle = 170
#between 0 and 180
R_max_protrusion_left = 300
#smaller than I_length
R_max_height = 3200
R_width = 1000

#Offset
Offset_x = 10
Offset_y = 10


mask_file_name = "files/mask1.dat"

mask_number = 2

global x_max
global y_max
global total_area
total_area = 0

def create_mask_filename():
	global mask_file_name
	global mask_number
	while(path.exists(mask_file_name)):
		mask_file_name = "files/mask" + str(mask_number) + ".dat"
		mask_number += 1
	mask_number -= 1

def write_parameters_to_database():
	f = open("mask_database.txt", "a")
	f.write(mask_file_name + " uses the following parameters:\n\n")
	f.write("Left shape:\n")
	f.write("L_shape = " + L_shape + "\n")
	if (L_shape == "circle"):
		f.write("L_radius = " + str(L_radius) + " pixels or " + str(L_radius*pixel_size) + " mm \n")
	if (L_shape == "rectangle"):
		f.write("L_width = " + str(L_width) + " pixels or " + str(L_width*pixel_size) + " mm \n")
		f.write("L_height = " + str(L_height) + " pixels or " + str(L_height*pixel_size) + " mm \n")
	f.write("\n")

	f.write("Isthmus: \n")
	f.write("I_length = " + str(I_length) + " pixels or " + str(I_length*pixel_size) + " mm \n")
	f.write("I_width = " + str(I_width) + " pixels or " + str(I_width*pixel_size) + " mm \n\n")

	f.write("Right shape: \n")
	f.write("R_angle = " + str(R_angle) + "\n")
	f.write("R_max_protrusion_left = " + str(R_max_protrusion_left) + " pixels or " + str(R_max_protrusion_left*pixel_size) + " mm \n")
	f.write("R_max_height = " + str(R_max_height) + " pixels or " + str(R_max_height*pixel_size) + " mm \n")
	f.write("R_width = " + str(R_width) + " pixels or " + str(R_width*pixel_size) + " mm \n\n")

	f.write("Offset: \n")
	f.write("Offset_x = " + str(Offset_x) + " pixels or " + str(Offset_x*pixel_size) + " mm \n")
	f.write("Offset_y = " + str(Offset_y) + " pixels or " + str(Offset_y*pixel_size) + " mm \n\n")

	f.write("This results in:\n")
	f.write("x_max = " + str(x_max) + " pixels or " + str(x_max*pixel_size) + " mm \n")
	f.write("y_max = " + str(y_max) + " pixels or " + str(y_max*pixel_size) + " mm \n\n")

	f.write("total_area = " + str(total_area) + " pixels or " + str(total_area*pixel_size**2) + " mm^2 \n\n")
	f.write("--------------------------------------------------------------------------- \n\n")
	f.close()

def write_parameters_to_parameter_file():
	f = open("parameters/mask_" + str(mask_number) + ".py", "w")
	if (L_shape == "circle"):
		f.write("L_radius = " + str(L_radius) + "\n")
	if (L_shape == "rectangle"):
		f.write("L_width = " + str(L_width) + "\n")
		f.write("L_height = " + str(L_height) + "\n")
	
	f.write("I_length = " + str(I_length) + "\n")
	f.write("I_width = " + str(I_width) + "\n")

	f.write("R_angle = " + str(R_angle) + "\n")
	f.write("R_max_protrusion_left = " + str(R_max_protrusion_left) + "\n")
	f.write("R_max_height = " + str(R_max_height) + "\n")
	f.write("R_width = " + str(R_width) + "\n")

	f.write("Offset_x = " + str(Offset_x) + "\n")
	f.write("Offset_y = " + str(Offset_y) + "\n")

	f.write("x_max = " + str(x_max) + "\n")
	f.write("y_max = " + str(y_max) + "\n")

	f.write("total_area = " + str(total_area) + "\n")
	f.close()

def find_max_xy():
	global x_max
	global y_max
	if (L_shape == "circle"):
		x_max = Offset_x*2 + L_radius*2 + I_length + R_width
		if R_angle <= 90:
			y_max = int(Offset_y*2 + max(I_width, L_radius*2, min(R_max_height, 2*math.tan(math.radians(R_angle))*R_width+I_width)))
		else:
			y_max = int(Offset_y*2 + max(I_width, L_radius*2, R_max_height))
	if (L_shape == "rectangle"):
		x_max = Offset_x*2 + L_width + I_length + R_width
		if R_angle <= 90:
			y_max = int(Offset_y*2 + max(I_width, L_height, min(R_max_height, 2*math.tan(math.radians(R_angle))*R_width+I_width)))
		else:
			y_max = int(Offset_y*2 + max(I_width, L_height, R_max_height)) 

def construct_shape():
	global total_area
	f = open(mask_file_name, "w")
	global mask
	if (L_shape == "circle"):
		centre_circle_x = Offset_x + L_radius
		centre_circle_y = y_max/2
		if R_angle <= 90:
			centre_wedge_x = Offset_x + centre_circle_x + I_length + L_radius- (I_width/2)/math.tan(math.radians(R_angle))
			for x in range(x_max):
				for y in range(y_max):
					if ((math.sqrt((x-centre_circle_x)**2 + (y-centre_circle_y)**2) < L_radius) or #if within circle
						(((x > centre_circle_x + Offset_x) and x <= (Offset_x + centre_circle_x + I_length + L_radius)) and ((y > y_max/2 - I_width/2) and (y < y_max/2 + I_width/2))) or #if within isthmus
						((x <= x_max - Offset_x) and (y > y_max/2 - R_max_height/2) and (y <= y_max/2 + R_max_height/2) and #if within wedge
						(y < y_max/2 + (x-centre_wedge_x)*math.tan(math.radians(R_angle))) and (y > y_max/2 - (x-centre_wedge_x)*math.tan(math.radians(R_angle))))):
						f.write(str(x) + "," + str(y) + "\n")
						mask[x,y]= 1
						total_area = total_area + 1

		elif R_angle > 90:
			centre_wedge_x = Offset_x + centre_circle_x + I_length + L_radius + (I_width/2)/math.tan(math.radians(180-R_angle))
			for x in range(x_max):
				for y in range(y_max):
					if ((math.sqrt((x-centre_circle_x)**2 + (y-centre_circle_y)**2) < L_radius) or #if within circle
						(((x > centre_circle_x + Offset_x) and x <= (Offset_x + centre_circle_x + I_length + L_radius)) and ((y > y_max/2 - I_width/2) and (y < y_max/2 + I_width/2))) or #if within isthmus
						((x > Offset_x + centre_circle_x + L_radius + I_length) and (x <= x_max - Offset_x) and (y > y_max/2 - R_max_height/2) and (y <= y_max/2 + R_max_height/2)) or #if within rectangle
						((x > (Offset_x + centre_circle_x +  L_radius + I_length - R_max_protrusion_left)) and (x <= (Offset_x + centre_circle_x + I_length + L_radius)) and #if within backward wedge
						((y > y_max/2 + (centre_wedge_x-x)*math.tan(math.radians(180-R_angle))) or (y < y_max/2 - (centre_wedge_x-x)*math.tan(math.radians(180-R_angle)))) and (y > Offset_y) and (y <= y_max - Offset_y))):
						f.write(str(x) + "," + str(y) + "\n")
						mask[x,y]= 1
						total_area = total_area + 1
	elif (L_shape == "rectangle"):
			if R_angle <= 90:
				centre_wedge_x = Offset_x + L_width + I_length - (I_width/2)/math.tan(math.radians(R_angle))
				for x in range(x_max):
					for y in range(y_max):
						if ((x > Offset_x and x <= Offset_x + L_width and y > y_max/2 - L_height/2 and y <= y_max/2 + L_height/2) or #if within left rectangle
							(((x > L_width + Offset_x) and x <= (Offset_x + L_width+ I_length)) and ((y > y_max/2 - I_width/2) and (y < y_max/2 + I_width/2))) or #if within isthmus
							((x <= x_max - Offset_x) and (y > y_max/2 - R_max_height/2) and (y <= y_max/2 + R_max_height/2) and #if within wedge
							(y < y_max/2 + (x-centre_wedge_x)*math.tan(math.radians(R_angle))) and (y > y_max/2 - (x-centre_wedge_x)*math.tan(math.radians(R_angle))))):
							f.write(str(x) + "," + str(y) + "\n")
							mask[x,y]= 1
							total_area = total_area + 1

			elif R_angle > 90:
				centre_wedge_x = Offset_x + L_width + I_length + (I_width/2)/math.tan(math.radians(180-R_angle))
				for x in range(x_max):
					for y in range(y_max):
						if ((x > Offset_x and x <= Offset_x + L_width and y > y_max/2 - L_height/2 and y <= y_max/2 + L_height/2) or #if within left rectangle
							(((x > L_width + Offset_x) and x <= (Offset_x + L_width + I_length)) and ((y > y_max/2 - I_width/2) and (y < y_max/2 + I_width/2))) or #if within isthmus
							((x > Offset_x + L_width + I_length) and (x <= x_max - Offset_x) and (y > y_max/2 - R_max_height/2) and (y <= y_max/2 + R_max_height/2)) or #if within right rectangle
							((x > (Offset_x + L_width + I_length - R_max_protrusion_left)) and (x <=(Offset_x + L_width + I_length)) and #if within backward wedge
							((y > y_max/2 + (centre_wedge_x-x)*math.tan(math.radians(180-R_angle))) or (y < y_max/2 - (centre_wedge_x-x)*math.tan(math.radians(180-R_angle)))) and (y > Offset_y) and (y <= y_max - Offset_y))):
							f.write(str(x) + "," + str(y) + "\n")
							mask[x,y]= 1
							total_area = total_area + 1

	f.close()

def create_image():
	plt.xticks(np.arange(0,x_max,400), np.arange(0,pixel_size*x_max,1))
	plt.yticks(np.arange(0,y_max,400), np.arange(0,pixel_size*y_max,1))
	plt.xlabel('mm')
	plt.ylabel('mm')
	plt.imshow(mask.T, interpolation='nearest')
	plt.savefig("images/mask" + str(mask_number) + ".jpg")



create_mask_filename()
find_max_xy()
mask = np.zeros ((x_max,y_max))
construct_shape()
create_image()
write_parameters_to_database()
write_parameters_to_parameter_file()




			
