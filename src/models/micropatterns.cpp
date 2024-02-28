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
#include <chrono>
#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <thread>
#include <math.h>
#include "dish.hpp"
#include "random.hpp"
#include "cell.hpp"
#include "info.hpp"
#include "parameter.hpp"
#include "sqr.hpp"
#include "graph.hpp"
#include "plotter.hpp"
#include "profiler.hpp"

using namespace std;



INIT {
  try {
    // Define initial distribution of cells
    //CPM->GrowInCellsInMicropattern(par.n_init_cells,par.size_init_cells);
    //CPM->ConstructInitCells(*this);

    //CPM->GrowCellGrid(*this);
    if (par.second_layer)
      CPM->GrowCellGridOnLayers(*this);
    else
      CPM->GrowCellGrid(*this);
    if (par.initial_configuration_file != "None")
      io->ReadConfiguration();
    else 
      CPM->ConstructInitCellGrid(*this);
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file using 
    // parameters n_init_cells, size_init_cells and divisions
    for (int i=0;i<par.divisions;i++) {
      CPM->DivideCells();
    }
    
    
    
    // Assign a random type to each of the cells

    //CPM->SetTypesWithMask();
    //CPM->SetRandomTypes();
    CPM->SetUpTauMatrix(par.sizex,par.sizey);
    CPM->InitializeEdgeList();
    //CPM->InitializeCouplingCoefficient_Gradient();
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);

  
  }
}

TIMESTEP { 
  try {
    static int i=0;
    cout <<"MCS: " << i << endl;
    static Dish *dish;
    if (i == 0 ){
        dish=new Dish();
        dish->CPM->InitializeCouplingCoefficient();
        dish->PDEfield->InitializePDEs(dish->CPM);

    }
    static Info *info=new Info(*dish, *this);
    static Plotter * plotter = new Plotter(dish, this);

    if (!(i%par.t_interval) && par.conv_ext){
		  //refresh links each par.t_interval time steps
		  dish->CPM->RefreshLinks();
	  }

    if (i == par.mcs-1 && par.relaxation > 0){
      string configuration_output_file = par.datadir + "/output_config_" + to_string(par.run_number) + ".json";
      cout << configuration_output_file << endl;
      dish->io->WriteConfiguration(&configuration_output_file[0]);
      char fname[200];
      sprintf(fname,"%s/output_config_%i.png",par.datadir.c_str(),par.run_number);
      PROFILE(plotter_2, plotter->Plot();)
      Write(fname);

    }

     //uncomment for chemotaxis
    if (i>=par.relaxation) {
      if(par.usecuda == true){
          dish->PDEfield->cuPDEsteps(dish->CPM, par.pde_its);
      }
      else if(par.useopencl == true){
        throw std::runtime_error("You are using openCL! This is probably an error.");
        PROFILE(opencl_diff, dish->PDEfield->ODEstepCL(dish->CPM, par.pde_its);)
      }

      else{
        throw std::runtime_error("You are computing secrete / diffuse! This is probably an error.");
        for (int r=0;r<par.pde_its;r++) {
          dish->PDEfield->Secrete(dish->CPM);
          dish->PDEfield->Diffuse(1);
        }
      }
    }


    
    if (par.graphics && !(i%par.storage_stride)) {
      PROFILE(all_plots, plotter->Plot();)
      info->Menu();
      dish->CPM->SetBoundingBox();
    }
    
    if (i == 0 && par.pause_on_start){ 
      info->set_Paused();
    i++;}

    if ((!info->IsPaused() && i < par.relaxation)){ //added second condition for test
      PROFILE(amoebamove, dish->CPM->AmoebaeMove(dish->PDEfield);)
    }


    if ( i == par.mcs){
      dish->ExportMultiCellDS(par.mcds_output);
    }

    if (par.store && !(i%par.storage_stride)) {
      char fname[200];
      sprintf(fname,"%s/image%07d.png",par.datadir.c_str(),i);
      PROFILE(plotter_2, plotter->Plot();)
      Write(fname);
    }

    if (!info->IsPaused()){
      i++;
    }
  } catch(const char* error) {
      cerr << "Caught exception\n";
      std::cerr << error << "\n";
      exit(1);
  }
  PROFILE_PRINT
}

void PDE::Secrete(CellularPotts *cpm) {
  const double dt=par.dt;
  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) {
      // inside cells
      if (cpm->Sigma(x,y)) {
	PDEvars[x*sizey+y]+=par.secr_rate[0]*dt;
      } else {
      // outside cells
	PDEvars[x*sizey+y]-=par.decay_rate[0]*dt*PDEvars[x*sizey+y];
      }
    }
  }
  PROFILE_PRINT
}

void Plotter::Plot()  {
  graphics->BeginScene();
  graphics->ClearImage(); 
  
  //Somewhere here show mask
  plotPDEDensity();
  plotCPMCellTypes();
  plotCPMLines(); 
  plotPDEContourLines();
  graphics->EndScene();
}

/* //Paci
int PDE::MapColour(double val) {
  if ((int(val*1000 + 100)/2) +155>255)
    return 255;
  else if ((int(val*1000 + 100)/2) +155<156)
    return 156;
  else  
    return (int(val*1000 + 100)/2) +155;
} */


 //FHN
int PDE::MapColour(double val) {
  double max_value = 50;
  double min_value = -80;
  if (val > max_value)
    return 255;
  else if (val < min_value)
    return 156;
  else{ 
    return (int)((val-min_value)/(max_value-min_value)*100) + 155;
  }
}


int main(int argc, char *argv[]) {
  
  extern Parameter par;
  try {  
    par.Read(argv[1]);
    long seed = Seed(par.rseed);
    ofstream metadatafile;
    metadatafile.open("metadata.txt", std::ofstream::out | std::ofstream::app);
    metadatafile << "Random seed = " << seed << "\n";
    metadatafile.close();
    start_graphics(argc, argv);
  } catch(const char* error) {
    std::cerr << error << std::endl;
    return 1;
  } /*catch(...) {
    std::cerr << "An unknown exception was caught" << std::endl;
    return 1;
  }*/

  return 0;
}
