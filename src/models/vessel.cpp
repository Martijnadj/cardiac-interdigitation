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
#include <stdio.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <thread>
#include "dish.hpp"
#include "random.hpp"
#include "cell.hpp"
#include "info.hpp"
#include "parameter.hpp"
#include "plotter.hpp"
#include "profiler.hpp"
#include "graph.hpp"

using namespace std;

INIT {
  try {
    // Define initial distribution of cells
     CPM->GrowInCells(par.n_init_cells,par.size_init_cells,par.subfield);
     CPM->ConstructInitCells(*this);
    
    // If we have only one big cell and divide it a few times
    // we start with a nice initial clump of cells. 
    // 
    // The behavior can be changed in the parameter file.
    for (int i=0;i<par.divisions;i++) {
      CPM->DivideCells();
    }

    CPM->InitializeEdgeList();
        
  } catch(const char* error) {
    cerr << "Caught exception\n";
    std::cerr << error << "\n";
    exit(1);
  }

}

TIMESTEP { 
  try {
    static int i=0;
    static Dish *dish=new Dish();
    static Info *info=new Info(*dish, *this);
    static Plotter plotter = Plotter(dish, this);
    if (i>=par.relaxation) {
      if (par.useopencl){
        PROFILE(opencl_diff, dish->PDEfield->ODEstepCL(dish->CPM, par.pde_its);)
      }
      else{
        for (int r=0;r<par.pde_its;r++) {
	  dish->PDEfield->Secrete(dish->CPM);
	  dish->PDEfield->Diffuse(1);
        }
      }
    }
    PROFILE(amoebamove, dish->CPM->AmoebaeMove(dish->PDEfield);)
    
    if (par.graphics && !(i%par.storage_stride)) {
      PROFILE(all_plots, plotter.Plot();) 
      char title[400];
      snprintf(title,399,"CellularPotts: %.2f hr",dish->PDEfield->TheTime()/3600);
      info->Menu();
    }
    if (par.store && !(i%par.storage_stride)) {
      char fname[200],fname_mcds[200];
      snprintf(fname,199,"%s/extend%05d.png",par.datadir,i);
      snprintf(fname_mcds,199,"%s/extend%05d.xml",par.datadir,i);
      Write(fname);
      if (!(i%(par.storage_stride*10)))
        dish->ExportMultiCellDS(fname_mcds);
    }
    i++;
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
	sigma[0][x][y]+=par.secr_rate[0]*dt;
      } else {
      // outside cells
	sigma[0][x][y]-=par.decay_rate[0]*dt*sigma[0][x][y];
      }
    }
  }
  PROFILE_PRINT
}

int PDE::MapColour(double val) {
  return (((int)((val/((val)+1.))*100))%100)+155;
}

void Plotter::Plot()  {
  graphics->BeginScene();
  graphics->ClearImage();
  
  plotPDEDensity();
  plotCPMCellTypes();
  plotCPMLines(); 
  plotPDEContourLines();
  
  graphics->EndScene();
}

int main(int argc, char *argv[]) {
  extern Parameter par;
  try {  
    par.Read(argv[1]);
    Seed(par.rseed);
    start_graphics(argc, argv);
  } catch(const char* error) {
    std::cerr << error << std::endl;
    return 1;
  }
  return 0;
}
