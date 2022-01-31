/* 

Copyright 1996-2006 Roeland Merks

This file is part of Tissue Simulation Toolkit.

Tissue Simulation Toolkit is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

Tissue Simulation Toolkit is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Tissue Simulation Toolkit; if not, write to the Free
Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
02110-1301 USA

*/


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cerrno>
#include <iostream>
#include "parameter.hpp"
#include "output.hpp"
#include "parse.hpp"

Parameter::Parameter() {

  T = 50.;
  target_area = 100;
  ref_adhesive_area = 100;
  target_length = 60;
  lambda = 50;
  lambda2 = 5.0;
  Jtable = strdup("J.dat");
  conn_diss = 2000;
  cluster_connectivity = false;
  vecadherinknockout = false;
  extensiononly = false;
  chemotaxis = 1000;
  border_energy = 100;
  neighbours = 2;
  periodic_boundaries = false;
  gradient = false;
  extended_neighbour_border = false;
  n_chem = 1;
  diff_coeff = new double[1];
  diff_coeff[0] = 1e-13;
  decay_rate = new double[1];
  decay_rate[0] = 1.8e-4;
  secr_rate = new double[1];
  secr_rate[0] = 1.8e-4;
  saturation = 0;
  dt = 2.0;
  min_stepsize = 1e-6;
  dx = 2.0e-6;
  pde_its = 15;
  micropatternmask = strdup("None");
  n_init_cells = 100;
  size_init_cells = 10;
  sizex = 200;
  sizey = 200;
  divisions = 0;
  mcs = 10000;
  rseed = -1;
  subfield = 1.0;
  relaxation = 0;
  storage_stride = 10;
  graphics = true;
  store = false;
  load_mcds = false;
  datadir = strdup("data_film");
  mcds_output = strdup("outstate.xml");
  mcds_input = strdup("false");
  mcds_anneal_steps = 0;
  mcds_denoise_steps = 0; 
  pause_on_start = false;
  useopencl = true;
  usecuda = true;
  number_of_cores = 1;
  threads_per_core = 1;

  opencl_core_path = strdup("../src/reaction_diffusion/pdecore.cl");
  opencl_pref_platform = 0;
  colortable = strdup("../data/default.ctb");

  initial_E_m = -70;
  initial_n = 0.5;
  initial_m = 0.25;
  initial_h = 0.08;
  eps = 1e-8;

  couplingmedium = 1e-5;
  couplingcell = 1e-5;
  couplingboundary = 1e-8;
  couplingoffmask = 0;

}

Parameter::~Parameter() {
  // destruct parameter object
  // free string parameter
  CleanUp();
}

void Parameter::CleanUp(void) {
  if (Jtable) 
     free(Jtable);
  if (diff_coeff) 
     free(diff_coeff);
  if (decay_rate) 
     free(decay_rate);
  if (secr_rate) 
     free(secr_rate);
  if (datadir) 
     free(datadir);
  if (micropatternmask) 
     free(micropatternmask);
}


void Parameter::Read(const char *filename) {
  static bool ReadP=false;
  if (ReadP) {
    //throw "Run Time Error in parameter.cpp: Please Read parameter file only once!!";
    CleanUp();
  } else
    ReadP=true;

  FILE *fp=OpenReadFile(filename);

  T = fgetpar(fp, "T", 50., true);
  target_area = igetpar(fp, "target_area", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  Jtable = sgetpar(fp, "Jtable", "J.dat", true);
  conn_diss = igetpar(fp, "conn_diss", 2000, true);
  ref_adhesive_area = igetpar(fp, "ref_adhesive_area", 100, true);
  target_length = igetpar(fp, "target_length", 60, true);
  lambda = fgetpar(fp, "lambda", 50, true);

  area_constraint_type = igetpar(fp, "area_constraint_type",0, true);
  lambda2 = fgetpar(fp, "lambda2", 5.0, true);
  Jtable = sgetpar(fp, "Jtable", "J.dat", true);
  conn_diss = igetpar(fp, "conn_diss", 2000, true);
  cluster_connectivity = bgetpar(fp, "cluster_connectivity", false, true);
  vecadherinknockout = bgetpar(fp, "vecadherinknockout", false, true);
  extensiononly = bgetpar(fp, "extensiononly", false, true);
  chemotaxis = igetpar(fp, "chemotaxis", 1000, true);
  border_energy = igetpar(fp, "border_energy", 100, true);
  neighbours = igetpar(fp, "neighbours", 2, true);
  periodic_boundaries = bgetpar(fp, "periodic_boundaries", false, true);
  gradient  = bgetpar(fp, "gradient", false, true);
  extended_neighbour_border = bgetpar(fp, "extended_neighbour_border", false, true);
  n_chem = igetpar(fp, "n_chem", 1, true);
  diff_coeff = dgetparlist(fp, "diff_coeff", n_chem, true);
  decay_rate = dgetparlist(fp, "decay_rate", n_chem, true);
  secr_rate = dgetparlist(fp, "secr_rate", n_chem, true);
  saturation = fgetpar(fp, "saturation", 0, true);
  dt = fgetpar(fp, "dt", 2.0, true);
  min_stepsize = fgetpar(fp, "min_stepsize", 1e-6, true);
  dx = fgetpar(fp, "dx", 2.0e-6, true);
  pde_its = igetpar(fp, "pde_its", 15, true);
  micropatternmask = sgetpar(fp, "micropatternmask", "None", true);
  n_init_cells = igetpar(fp, "n_init_cells", 100, true);
  size_init_cells = igetpar(fp, "size_init_cells", 10, true);
  sizex = igetpar(fp, "sizex", 200, true);
  sizey = igetpar(fp, "sizey", 200, true);
  divisions = igetpar(fp, "divisions", 0, true);
  mcs = igetpar(fp, "mcs", 10000, true);
  rseed = igetpar(fp, "rseed", -1, true);
  subfield = fgetpar(fp, "subfield", 1.0, true);
  relaxation = igetpar(fp, "relaxation", 0, true);
  storage_stride = igetpar(fp, "storage_stride", 10, true);
  graphics = bgetpar(fp, "graphics", true, true);
  store = bgetpar(fp, "store", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);
  load_mcds = bgetpar(fp, "load_mcds", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);
  mcds_output = sgetpar(fp, "mcds_output", "outstate.xml", true);
  mcds_input = sgetpar(fp, "mcds_input", "outstate.xml", true);
  mcds_anneal_steps = igetpar(fp, "mcds_anneal_steps", true);
  mcds_denoise_steps = igetpar(fp, "mcds_denoise_steps", true);
  pause_on_start = bgetpar(fp, "pause_on_start", true);
  useopencl = bgetpar(fp, "useopencl", true, true);
  usecuda = bgetpar(fp, "usecuda", true, true);
  number_of_cores = igetpar(fp, "number_of_cores", true);
  threads_per_core = igetpar(fp, "threads_per_core", true);
  opencl_core_path = sgetpar(fp, "opencl_core_path", "../src/reaction_diffusion/pdecore.cl", true);
  opencl_pref_platform = igetpar(fp, "opencl_pref_platform", 0, true);
  adhesion_storage_stride = igetpar(fp, "adhesion_storage_stride", mcs+1, true);
  graphics = bgetpar(fp, "graphics", true, true);
  store = bgetpar(fp, "store", false, true);
  datadir = sgetpar(fp, "datadir", "data_film", true);
  // Act model
  lambda_Act = fgetpar(fp, "lambda_Act", 0, true);
  max_Act = igetpar(fp, "max_Act", 0, true);
  lambda_perimeter = igetpar(fp, "lambda_perimeter",0, true);
  target_perimeter = igetpar(fp, "target_perimeter", 0, true);
  //lymphocyte matrix interaction
  lambda_matrix = fgetpar(fp, "lambda_matrix", 0, true);
  spontaneous_p = fgetpar(fp, "spontaneous_p", 0.001, true);
  decay_p = fgetpar(fp,"decay_p", 0.02, true);
  eden_p = fgetpar(fp,"eden_p", 0.25, true);
  J_pol = igetpar( fp,"J_pol", 0 , true);
  threshold = fgetpar(fp, "threshold", 0 , true);
  start_level = fgetpar(fp, "start_level", 1, true);
  colortable = sgetpar(fp, "colortable", "../data/default.ctb", true);

  initial_E_m = fgetpar(fp, "initial_E_m", -70, true);
  initial_n = fgetpar(fp, "initial_n", 0.5, true);
  initial_m = fgetpar(fp, "initial_m", 0.25, true);
  initial_h = fgetpar(fp, "initial_h", 0.08, true);
  eps = fgetpar(fp, "eps", 1e-8, true);

  couplingmedium = fgetpar(fp, "couplingmedium", 1e-5, true);
  couplingcell = fgetpar(fp, "couplingcell", 1e-5, true);
  couplingboundary = fgetpar(fp, "couplingboundary", 1e-8, true);
  couplingoffmask = fgetpar(fp, "couplingoffmask", 0, true);

  std::cout << std::endl;
}

const char *sbool(const bool &p) {

  const char *true_str="true";
  const char *false_str="false";
  if (p)
    return true_str;
  else
    return false_str;
}

void Parameter::Write(ostream &os) const {
  setlocale(LC_NUMERIC, "C");

  os << " T = " << T << endl;
  os << " target_area = " << target_area << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " lambda2 = " << lambda2 << endl;
  if (Jtable) 
    os << " Jtable = " << Jtable << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " ref_adhesive_area = " << ref_adhesive_area << endl;
  os << " target_length = " << target_length << endl;
  os << " lambda = " << lambda << endl;
  os << " area_constraint_type = " << area_constraint_type << endl;
  os << " lambda2 = " << lambda2 << endl;
  if (Jtable)
    os << " Jtable = " << Jtable << endl;
  os << " conn_diss = " << conn_diss << endl;
  os << " cluster_connectivity = " << cluster_connectivity << endl;
  os << " vecadherinknockout = " << sbool(vecadherinknockout) << endl;
  os << " extensiononly = " << sbool(extensiononly) << endl;
  os << " chemotaxis = " << chemotaxis << endl;
  os << " border_energy = " << border_energy << endl;
  os << " neighbours = " << neighbours << endl;
  os << " periodic_boundaries = " << sbool(periodic_boundaries) << endl;
  os << " gradient = " << sbool(gradient) << endl;
  os << " extended_neighbour_border = " << sbool(extended_neighbour_border) << endl;
  os << " n_chem = " << n_chem << endl;
  os << " diff_coeff = "<< diff_coeff[0] << endl;
  os << " decay_rate = "<< decay_rate[0] << endl;
  os << " secr_rate = "<< secr_rate[0] << endl;
  os << " saturation = " << saturation << endl;
  os << " dt = " << dt << endl;
  os << " min_stepsize = " << min_stepsize << endl;
  os << " dx = " << dx << endl;
  os << " pde_its = " << pde_its << endl;
  os << " micropatternmask = " << micropatternmask << endl;
  os << " n_init_cells = " << n_init_cells << endl;
  os << " size_init_cells = " << size_init_cells << endl;
  os << " sizex = " << sizex << endl;
  os << " sizey = " << sizey << endl;
  os << " divisions = " << divisions << endl;
  os << " mcs = " << mcs << endl;
  os << " rseed = " << rseed << endl;
  os << " subfield = " << subfield << endl;
  os << " relaxation = " << relaxation << endl;
  os << " storage_stride = " << storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;
  os << " load_mcds = " << sbool(load_mcds) << endl;
  os << " mcds_output = " << mcds_output << endl;
  os << " mcds_input = " << mcds_input << endl;
  os << " mcds_anneal_steps = " << mcds_anneal_steps << endl;
  os << " mcds_denoise_steps = " << mcds_denoise_steps << endl;
  os << " pause_on_start = " << pause_on_start << endl; 
  os << " useopencl = " << useopencl << endl;
  os << " usecuda = " << usecuda << endl;
  os << " number_of_cores = " << number_of_cores << endl;
  os << " threads_per_core = " << threads_per_core << endl;
  os << " opencl_pref_platform" << opencl_pref_platform << endl;
  os << " opencl_core_path = " << opencl_core_path << endl;
  if (datadir) 
  os << " adhesion_storage_stride = " << adhesion_storage_stride << endl;
  os << " graphics = " << sbool(graphics) << endl;
  os << " store = " << sbool(store) << endl;
  //Act model
  os << " lambda_Act = " << lambda_Act << endl;
  os << " max_Act = " << max_Act << endl;
  os << " lambda_perimeter = " << lambda_perimeter << endl;
  os << " target_perimeter = " << target_perimeter << endl;
  //lymphocyte matrix interaction
  os << " lambda_matrix = " << lambda_matrix << endl;
  os << " spontaneous_p = " << spontaneous_p << endl;
  os << " decay_p = " << decay_p << endl;
  os << " eden_p = " << eden_p << endl;
  os << " J_pol = " << J_pol << endl;
  if (datadir)
    os << " datadir = " << datadir << endl;
  os << " colortable = " << colortable << endl;

  os << " initial_E_m = " << initial_E_m << endl;
  os << " initial_n = " << initial_n << endl;
  os << " initial_m = " << initial_m << endl;
  os << " initial_h = " << initial_h<< endl;
  os << " eps = " << eps<< endl;

  os << " couplingmedium = " << couplingmedium << endl;
  os << " couplingcell = " << couplingcell << endl;
  os << " couplingboundary = " << couplingboundary << endl;
  os << " couplingoffmask = " << couplingoffmask << endl;
}

ostream &operator<<(ostream &os, Parameter &p) {
  p.Write(os);
  return os;
}

Parameter par;
