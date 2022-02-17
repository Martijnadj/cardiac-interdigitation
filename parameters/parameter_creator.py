#creates a parameter file for a mask file, takes mask_i as input
import sys
import importlib
from os import path

mask_parameter_file = sys.argv[1]

global par
fraction_covered = 0.8
cellsize = 192
def import_mask_file():
    global par
    par = importlib.import_module("masks.parameters."+mask_parameter_file)
    
def write_parameter_file(): 
    filename = mask_parameter_file
    filename = filename+ "_standard.par"
    if(path.exists(filename)):
        print("Filename already exists, abort")
        exit(1)
    f = open(filename, "w")
    f.write("# Cellular Potts parameters \n")
    f.write("T = 50 \n")
    f.write("target_area = " + str(cellsize) + " \n")
    f.write("target_length = 24 \n")
    f.write("lambda = 50 \n")
    f.write("lambda2 = 5 \n")
    f.write("Jtable = ../data/J.dat \n")
    f.write("conn_diss = 2000 \n")
    f.write("vecadherinknockout = true \n")
    f.write("chemotaxis = 0 \n")
    f.write("extensiononly = false \n")
    f.write("border_energy = 0 \n")
    f.write("\n")
    f.write("# note: do not change the following parameters for 'long' cells (lambda2>0) \n")
    f.write("neighbours = 2 \n")
    f.write("periodic_boundaries = false \n")
    f.write("\n")
    f.write("# PDE parameters \n")
    f.write("n_chem = 23 \n")
    f.write("diff_coeff = 1e-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 \n")
    f.write("decay_rate = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 \n")
    f.write("secr_rate = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 \n")
    f.write("saturation = 0 \n")
    f.write("dt = 1e-3 \n")
    f.write("min_stepsize = 1e-6 \n")
    f.write("dx = 2.5e-6 \n")
    f.write("pde_its = 1 \n")
    f.write("eps = 1e-6 \n")
    f.write("\n")
    f.write("\n")
    f.write("# Micropatterns parameters \n")
    f.write("micropatternmask = ../data/mask.dat \n")
    f.write("\n")
    f.write("# Cuda parameters \n")
    f.write("usecuda = true \n")
    f.write("number_of_cores = 80 \n")
    f.write("threads_per_core = 32 \n")
    f.write("\n")
    f.write("# initial conditions (create a 'blob' of cells in the middle) \n")
    f.write("n_init_cells = " + str(int(par.total_area*fraction_covered/cellsize)) + " \n")
    f.write("size_init_cells = 10 \n")
    f.write("sizex = " + str(par.x_max) + "\n")
    f.write("sizey = " + str(par.y_max) + "\n")
    f.write("divisions = 0 \n")
    f.write("mcs = 10 \n")
    f.write("rseed = -1 \n")
    f.write("subfield = 1.0 \n")
    f.write("relaxation = 0 \n")
    f.write("\n")
    f.write("couplingmedium = 0 \n")
    f.write("couplingcell = 5e-10 \n")
    f.write("couplingboundary = 1e-11 \n")
    f.write("couplingoffmask = 0 \n")
    f.write("\n")
    f.write("# output \n")
    f.write("storage_stride = 1 \n")
    f.write("graphics = true \n")
    f.write("store = true \n")
    f.write("datadir = ../images \n")


import_mask_file()
write_parameter_file()




			
