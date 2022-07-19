#creates a parameter file for a mask file, takes mask_i as input
import sys
import importlib
from os import path

mask_parameter_file = sys.argv[1]
noise = [0,0.05,0.1,0.15,0.2,0.25]
interval = [0.45]  

global par
fraction_covered = 1.1
def import_mask_file():
    global par
    par = importlib.import_module("masks.parameters."+mask_parameter_file)
    
def write_parameter_file(): 
    pixel_volume = par.pixel_size*par.pixel_size
    cellsize = int(192/pixel_volume*(0.0025*0.0025))
    target_length = int(24/par.pixel_size*0.0025)
    filename = mask_parameter_file
    filename = filename+ "_standard.par"
    if(path.exists(filename)):
        print("Filename already exists, abort")
        exit(1)
    f = open(filename, "w")
    f.write("# Cellular Potts parameters\n")
    f.write("T = 50\n")
    f.write("target_area = " + str(cellsize) + "\n")
    f.write("target_length = " + str(target_length) + "\n")
    f.write("lambda = 50\n")
    f.write("lambda2 = 5\n")
    f.write("Jtable = ../data/J2.dat\n")
    f.write("conn_diss = 2000\n")
    f.write("vecadherinknockout = true\n")
    f.write("chemotaxis = 0\n")
    f.write("extensiononly = false\n")
    f.write("border_energy = 0\n")
    f.write("\n")
    f.write("# note: do not change the following parameters for 'long' cells (lambda2>0) \n")
    f.write("neighbours = 1\n")
    f.write("periodic_boundaries = false\n")
    f.write("\n")
    f.write("# PDE parameters\n")
    f.write("n_chem = 33\n")
    f.write("diff_coeff = 1e-14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n")
    f.write("decay_rate = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n")
    f.write("secr_rate = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0\n")
    f.write("saturation = 0\n")
    f.write("dt = 2e-5\n")
    f.write("min_stepsize = 1e-18\n")
    f.write("dx = " + str(par.pixel_size*0.001)+"\n")
    f.write("pde_its = 50\n")
    f.write("eps = 1\n")
    f.write("\n")
    f.write("\n")
    f.write("# Micropatterns parameters \n")
    f.write("micropatternmask = ../parameters/masks/files/mask"+mask_parameter_file.lstrip("mask_") + ".dat\n")
    f.write("\n")
    f.write("# Cuda parameters \n")
    f.write("usecuda = true\n")
    f.write("number_of_cores = 272\n")
    f.write("threads_per_core = 32\n")
    f.write("\n")
    f.write("# initial conditions (create a 'blob' of cells in the middle)\n")
    f.write("n_init_cells = " + str(int(par.total_area*fraction_covered/cellsize)) + "\n")
    f.write("size_init_cells = 10\n")
    f.write("sizex = " + str(par.x_max) + "\n")
    f.write("sizey = " + str(par.y_max) + "\n")
    f.write("divisions = 0\n")
    f.write("mcs = 20100\n")
    f.write("rseed = -1\n")
    f.write("subfield = 1.0\n")
    f.write("relaxation = 100\n")
    f.write("\n")
    f.write("couplingmedium = 0\n")
    f.write("couplingAtrialAtrial = 1e-8\n")
    f.write("couplingAtrialPM = 5e-9\n")
    f.write("couplingPMPM = 1e-8\n")
    f.write("couplingcell = 2e-8\n")
    f.write("couplingoffmask = 0\n")
    f.write("\n")
    f.write("beats_per_minute = 60\n")
    f.write("pacing_duration = 5e-3\n")
    f.write("pacing_strength = 5.5e-10\n")
    f.write("\n")
    f.write("#Cell grid parameters \n")
    f.write("celltype1_width = 3 \n")
    f.write("celltype1_length = 20 \n")
    f.write("celltype2_width = 3 \n")
    f.write("celltype2_length = 4 \n")
    #f.write("#FitzHugh-Nagumo\n")
    #f.write("FHN_interval_beats = " + str(interval[0])+"\n")
    #f.write("FHN_pulse_duration = 0.1\n")	
    #f.write("FHN_pulse_strength = 5\n")
    #f.write("FHN_a = 0.3\n")
    #f.write("FHN_b = 0.999\n")
    #f.write("FHN_tau = 10\n")
    #f.write("FHN_a_diff_perc = " + str(noise[0])+"\n")
    #f.write("FHN_b_diff_perc = " + str(noise[0])+"\n")
    #f.write("FHN_tau_diff_perc = " + str(noise[0])+"\n")
    #f.write("FHN_start_0 = -0.575897\n")
    #f.write("FHN_start_1 = -3.84897\n")
    #f.write("\n")
    f.write("\n")
    f.write("# output\n")
    f.write("storage_stride = 1\n")
    f.write("graphics = true\n")
    f.write("store = true\n")
    f.write("datadir = ../images\n")


import_mask_file()
write_parameter_file()




		
