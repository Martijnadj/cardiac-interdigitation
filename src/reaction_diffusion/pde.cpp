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
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <chrono>

#include "crash.hpp"
#include "parameter.hpp"
#include "ca.hpp"
#include "pde.hpp"
#include "conrec.hpp"
#include "graph.hpp"

/* STATIC DATA MEMBER INITIALISATION */
const int PDE::nx[9] = {0, 1, 1, 1, 0,-1,-1,-1, 0 };
const int PDE::ny[9] = {0, 1, 0,-1,-1,-1, 0, 1, 1 };

extern Parameter par;

/** PRIVATE **/

PDE::PDE(const int l, const int sx, const int sy) {
  PDEvars=0;
  thetime=0;
  PDEsteps=0;
  sizex=sx;
  sizey=sy;
  layers=l;
  PDEvars=AllocatePDEvars(l,sx,sy);
  alt_PDEvars=AllocatePDEvars(l,sx,sy);
  btype=1; //boundary type
}


PDE::PDE(void) {

  PDEvars=0;
  alt_PDEvars=0;
  sizex=0; sizey=0; layers=0;
  thetime=0;
  if (par.useopencl){this->SetupOpenCL();}
}

// destructor (virtual)
PDE::~PDE(void) {
  if (PDEvars) {
    free(PDEvars[0][0]);
    free(PDEvars[0]);
    free(PDEvars);
    PDEvars=0;
  }
  if (alt_PDEvars) {
    free(alt_PDEvars[0][0]);
    free(alt_PDEvars[0]);
    free(alt_PDEvars);
    alt_PDEvars=0;
  }
}

PDEFIELD_TYPE ***PDE::AllocatePDEvars(const int layers, const int sx, const int sy) {
  PDEFIELD_TYPE ***mem;
  sizex=sx; sizey=sy;
  mem=(PDEFIELD_TYPE ***)malloc(layers*sizeof(PDEFIELD_TYPE **));
  if (mem==NULL) {
    MemoryWarning();
  }
  mem[0]=(PDEFIELD_TYPE **)malloc(layers*sizex*sizeof(PDEFIELD_TYPE *));
  if (mem[0]==NULL) { 
    MemoryWarning();
  }
  for (int i=1;i<layers;i++) {
    mem[i]=mem[i-1]+sizex;
  }  
  mem[0][0]=(PDEFIELD_TYPE *)malloc(layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (mem[0][0]==NULL) {
    MemoryWarning();
  }
  for (int i=1;i<layers*sizex;i++) {
    mem[0][i]=mem[0][i-1]+sizey;
  }
  /* Clear PDE plane */
  for (int i=0;i<layers*sizex*sizey;i++) {
    mem[0][0][i]=0.; 
  }
  return mem;
}

void PDE::InitializePDEvars(void){
    for (int i = 0; i< sizex* sizey;i++){

      PDEvars[0][0][i] = -0.050;
      PDEvars[1][0][i] = 0.32;
      PDEvars[2][0][i] = 0.0002;
      PDEvars[3][0][i] = 0;
      PDEvars[4][0][i] = 0;
      PDEvars[5][0][i] = 1;
      PDEvars[6][0][i] = 1;
      PDEvars[7][0][i] = 1;
      PDEvars[8][0][i] = 0;
      PDEvars[9][0][i] = 1;
      PDEvars[10][0][i] = 0;
      PDEvars[11][0][i] = 0.75;
      PDEvars[12][0][i] = 0.75; 
      PDEvars[13][0][i] = 0;
      PDEvars[14][0][i] = 0.1;
      PDEvars[15][0][i] = 1;
      PDEvars[16][0][i] = 0;
      PDEvars[17][0][i] = 9.2;
      PDEvars[18][0][i] = 0;
      PDEvars[19][0][i] = 0.75;
      PDEvars[20][0][i] = 0.3;
      PDEvars[21][0][i] = 0.9;
      PDEvars[22][0][i] = 0.1;
    }


}




void PDE::Plot(Graphics *g,const int l) {
  // l=layer: default layer is 0
  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g->Point(MapColour(PDEvars[l][x][y]),x,y);
      g->Point(MapColour(PDEvars[l][x][y]),x+1,y);
      g->Point(MapColour(PDEvars[l][x][y]),x,y+1);
      g->Point(MapColour(PDEvars[l][x][y]),x+1,y+1);
    }
  }
}

// Plot the value of the PDE only in the medium of the CPM
void PDE::Plot(Graphics *g, CellularPotts *cpm, const int l) {
  // suspend=true suspends calling of DrawScene
  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) { 
      if (cpm->Sigma(x,y)==0) {
	// Make the pixel four times as large
	// to fit with the CPM plane
        g->Point(MapColour(PDEvars[l][x][y]),x,y);
	g->Point(MapColour(PDEvars[l][x][y]),x+1,y);
	g->Point(MapColour(PDEvars[l][x][y]),x,y+1);
	g->Point(MapColour(PDEvars[l][x][y]),x+1,y+1);
      }
    }
  }
}

void PDE::ContourPlot(Graphics *g, int l, int colour) {
  // calls "conrec" routine by Paul Bourke, as downloaded from
  // http://astronomy.swin.edu.au/~pbourke/projection/conrec

  // number of contouring levels
  int nc = 10;

  // A one dimensional array z(0:nc-1) that saves as a list of the contour levels in increasing order.   
  double *z=(double *)malloc(nc*sizeof(double));
  double min=Min(l), max=Max(l);
  double step=(max-min)/nc;
  if (min == 0 && max == 0) return;

  for (int i=0;i<nc;i++) {
    z[i]=(i+1)*step;
  }
  double *x=(double *)malloc(sizex*sizeof(double));
  for (int i=0;i<sizex;i++) {
    x[i]=i;
  }
  double *y=(double *)malloc(sizey*sizeof(double));
  for (int i=0;i<sizey;i++) {
    y[i]=i;
  }

  conrec(PDEvars[l],0,sizex-1,0,sizey-1,x,y,nc,z,g,colour);
  
  free(x);
  free(y);
  free(z);
}



int MapColour3(double val, int l) {
	int step=0;
	if (l==2){
	step = (240)/par.max_Act;
  return (int)(256-val*step-1);}
	else if (l==3 && val>1){
	step = (256)/2;
  return (int)(500+val*step);}
	else
	return 0;
}


void PDE::PlotInCells (Graphics *g, CellularPotts *cpm, const int l) {
  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) {
      if( cpm->Sigma(x,y)>0) {
        if (par.lambda_Act>0) {
	   g->Rectangle(MapColour3(cpm->actPixels[{x,y}],l), x, y);
        } else {
	  g->Rectangle(255, x, y);
        }
        if (par.lambda_matrix>0) {
          if (cpm->matrix[x][y]>0) {
            g->Rectangle(256, x, y);
          }
        }
      } else if (cpm->Sigma(x,y)==-2) {
        g->Rectangle(10, x, y);
      }
      if (cpm->Sigma(x,y)==-3) {
	g->Rectangle(256, x, y);
      }
    }
  }
}



void PDE::SetupOpenCL(){
  extern CLManager clm;

  program = clm.make_program(par.opencl_core_path, OPENCL_PDE_TYPE);

  //Secretion and diffusion variables
  PDEFIELD_TYPE dt = (PDEFIELD_TYPE) par.dt;
  PDEFIELD_TYPE dx2 = (PDEFIELD_TYPE) par.dx*par.dx;
  // This ignores all but the first value?!
  PDEFIELD_TYPE decay_rate = (PDEFIELD_TYPE) * par.decay_rate.data();
  PDEFIELD_TYPE secr_rate = (PDEFIELD_TYPE) * par.secr_rate.data();
  
  
  int btype = 1;
  if (par.periodic_boundaries) btype=2;

  //Allocate memory on the GPU
  clm.cpm = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(int)*sizex*sizey); 
  clm.numberofedges = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(int)*sizex*sizey); 
  clm.couplingcoefficient = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(PDEFIELD_TYPE)*sizex*sizey); 
  clm.pdeA = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(PDEFIELD_TYPE)*sizex*sizey*layers);
  clm.pdeB = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(PDEFIELD_TYPE)*sizex*sizey*layers); 
  clm.diffco = cl::Buffer(clm.context, CL_MEM_READ_WRITE, sizeof(PDEFIELD_TYPE)*layers);


  //Making kernel and setting arguments
  kernel_ODEstep = cl::Kernel(program,"ODEstep");      

  kernel_ODEstep.setArg(0, clm.cpm);
  kernel_ODEstep.setArg(1, clm.pdeA);
  kernel_ODEstep.setArg(2, clm.pdeB);
  kernel_ODEstep.setArg(3, sizeof(int), &sizex);
  kernel_ODEstep.setArg(4, sizeof(int), &sizey);
  kernel_ODEstep.setArg(5, sizeof(int), &layers);
  kernel_ODEstep.setArg(6, sizeof(PDEFIELD_TYPE), &decay_rate);
  kernel_ODEstep.setArg(7, sizeof(PDEFIELD_TYPE), &dt);
  kernel_ODEstep.setArg(8, sizeof(PDEFIELD_TYPE), &dx2);
  kernel_ODEstep.setArg(9, clm.diffco);
  kernel_ODEstep.setArg(10,sizeof(PDEFIELD_TYPE), &secr_rate);
  kernel_ODEstep.setArg(11, sizeof(int),  &btype);
  kernel_ODEstep.setArg(12, clm.numberofedges);
  kernel_ODEstep.setArg(13, clm.couplingcoefficient);


  PDEFIELD_TYPE diff_coeff[layers];

  for (int index = 0; index < layers; index++){
    diff_coeff[index] = (PDEFIELD_TYPE) par.diff_coeff[index];
  }

  clm.queue.enqueueWriteBuffer(clm.diffco,
    CL_TRUE, 0, sizeof(PDEFIELD_TYPE)*layers, diff_coeff);

  openclsetup = true;
} 


void PDE::ODEstepCL(CellularPotts *cpm, int repeat){
    extern CLManager clm; 
    if (!openclsetup ){this->SetupOpenCL();}
    //A B scheme used to keep arrays on GPU
    int errorcode = 0;

    
    
    //Write the cellSigma array to GPU for secretion
    clm.queue.enqueueWriteBuffer(clm.cpm,
    CL_TRUE, 0, sizeof(int)*sizex*sizey, cpm->getSigma()[0]);
    clm.queue.enqueueWriteBuffer(clm.numberofedges,
    CL_TRUE, 0, sizeof(int)*sizex*sizey, cpm->getNumberofedges()[0]);
    clm.queue.enqueueWriteBuffer(clm.couplingcoefficient,
    CL_TRUE, 0, sizeof(PDEFIELD_TYPE)*sizex*sizey, cpm->getCouplingCoefficient()[0]);
    //Writing pdefield PDEvars is only necessary if modified outside of clm.pdeA)kernel
    if (first_round) {
      clm.queue.enqueueWriteBuffer(clm.pdeA,  CL_TRUE, 0, sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars[0][0]);
      first_round = false;
    }
    //Main loop executing kernel and switching between A and B arrays
    for (int index = 0; index < repeat; index ++){
      for (int innerloop = 0; innerloop < 4; innerloop++){

          using std::chrono::high_resolution_clock;
          using std::chrono::duration_cast;
          using std::chrono::duration;
          using std::chrono::milliseconds;

          auto t1 = high_resolution_clock::now();

        kernel_ODEstep.setArg(14, sizeof(int), &PDEsteps);
        if(innerloop == 0 || innerloop == 1){
          kernel_ODEstep.setArg(1, clm.pdeA);
          kernel_ODEstep.setArg(2, clm.pdeB);
        }
        else{
          kernel_ODEstep.setArg(1, clm.pdeB);
          kernel_ODEstep.setArg(2, clm.pdeA);
        }
        if(innerloop == 0 || innerloop == 2){
          errorcode = clm.queue.enqueueNDRangeKernel(kernel_ODEstep,
                      cl::NullRange, cl::NDRange(sizex*sizey), cl::NullRange);
          errorcode = clm.queue.finish();
        }
        else{
          errorcode = clm.queue.enqueueNDRangeKernel(kernel_ODEstep,
                      cl::NullRange, cl::NDRange(1), cl::NullRange);
          errorcode = clm.queue.finish();
        }
        if (errorcode != 0){
          printf("Error during OpenCL secretion and diffusion");
          exit(0);
        }

        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1;
        cout << "For PDEsteps = " << PDEsteps << ", " << ms_double.count() << " ms has elapsed" << endl;
        PDEsteps += 1;

      }
    }
    //Reading from correct array containing the output
    if (clm.pde_AB == 0) {
      clm.queue.enqueueReadBuffer(clm.pdeB,CL_TRUE,0,
                            sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars[0][0]);
    }
    else {
      clm.queue.enqueueReadBuffer(clm.pdeA,CL_TRUE,0,
                            sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars[0][0]);
    }
    if (errorcode != CL_SUCCESS) cout << "error:" << errorcode << endl;
    
}


// public
void PDE::Diffuse(int repeat) {
  
  // Just diffuse everywhere (cells are transparent), using finite difference
  // (We're ignoring the problem of how to cope with moving cell
  // boundaries right now)
  
  const PDEFIELD_TYPE dt=par.dt;
  const PDEFIELD_TYPE dx2=par.dx*par.dx;
  cout << "a" << endl;
  for (int r=0;r<repeat;r++) {
    //NoFluxBoundaries();
    cout << "b" << endl;
    if (par.periodic_boundaries) {
      PeriodicBoundaries();
    } else {
      AbsorbingBoundaries();
      //NoFluxBoundaries();
    }
    cout << "c" << endl;
    for (int l=0;l<layers;l++) {
      for (int x=1;x<sizex-1;x++)
        for (int y=1;y<sizey-1;y++) {
          
          PDEFIELD_TYPE sum=0.;
          sum+=PDEvars[l][x+1][y];
          sum+=PDEvars[l][x-1][y];
          sum+=PDEvars[l][x][y+1];
          sum+=PDEvars[l][x][y-1];
          sum-=4*PDEvars[l][x][y];
          alt_PDEvars[l][x][y]=PDEvars[l][x][y]+sum*dt*par.diff_coeff[l]/dx2;
        }
    }
    cout << "d" << endl;
    PDEFIELD_TYPE ***tmp;
    tmp=PDEvars;
    PDEvars=alt_PDEvars;
    alt_PDEvars=tmp;
  
    thetime+=dt;
  }
}

double PDE::GetChemAmount(const int layer) {
  // Sum the total amount of chemical in the lattice
  // in layer l
  // (This is useful to check particle conservation)
  double sum=0.;
  if (layer==-1) { // default argument: sum all chemical species
    for (int l=0;l<layers;l++) {
      for (int x=1;x<sizex-1;x++) {
	for (int y=1;y<sizey-1;y++) {
	  sum+=PDEvars[l][x][y];
	}
      }
    }
  } else {
    for (int x=1;x<sizex-1;x++)
      for (int y=1;y<sizey-1;y++) {
	sum+=PDEvars[layer][x][y];
      }
  } 
  return sum;
}

// private
void PDE::NoFluxBoundaries(void) {
  // all gradients at the edges become zero, 
  // so nothing flows out
  // Note that four corners points are not defined (0.)
  // but they aren't used in the calculations
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l][x][0]=PDEvars[l][x][1];
      PDEvars[l][x][sizey-1]=PDEvars[l][x][sizey-2];
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l][0][y]=PDEvars[l][1][y];
      PDEvars[l][sizex-1][y]=PDEvars[l][sizex-2][y];
    }
  }
}


// private
void PDE::AbsorbingBoundaries(void) {
  // all boundaries are sinks, 
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l][x][0]=0.;
      PDEvars[l][x][sizey-1]=0.;
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l][0][y]=0.;
      PDEvars[l][sizex-1][y]=0.;
    }
  }
}

// private
void PDE::PeriodicBoundaries(void) {
  // periodic...
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l][x][0]=PDEvars[l][x][sizey-2];
      PDEvars[l][x][sizey-1]=PDEvars[l][x][1];
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l][0][y]=PDEvars[l][sizex-2][y];
      PDEvars[l][sizex-1][y]=PDEvars[l][1][y];
    }
  }
}

void PDE::GradC(int layer, int first_grad_layer) {
  // calculate the first and second order gradients and put
  // them in the next chemical fields
  if (par.n_chem<5) {
    throw("PDE::GradC: Not enough chemical fields");
  }

  // GradX
  for (int y=0;y<sizey;y++) {
    for (int x=1;x<sizex-1;x++) {
      PDEvars[first_grad_layer][x][y]=(PDEvars[layer][x+1][y]-PDEvars[layer][x-1][y])/2.;
    } 
  }
  // GradY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      PDEvars[first_grad_layer+1][x][y]=(PDEvars[layer][x][y+1]-PDEvars[layer][x][y-1])/2.;
    } 
  }
  // GradXX
  for (int y=0;y<sizey;y++) {
    for (int x=1;x<sizex-1;x++) {
      PDEvars[first_grad_layer+2][x][y]=PDEvars[layer][x+1][y]-PDEvars[layer][x-1][y]-2*PDEvars[layer][x][y];
    } 
  }

  // GradYY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      PDEvars[first_grad_layer+3][x][y]=PDEvars[layer][x][y-1]-PDEvars[layer][x][y+1]-2*PDEvars[layer][x][y];
    } 
  }
}

void PDE::PlotVectorField(Graphics &g, int stride, int linelength, int first_grad_layer) {
  // Plot vector field assuming it's in layer 1 and 2
  for (int x=1;x<sizex-1;x+=stride) {
    for (int y=1;y<sizey-1;y+=stride) {
      
      // calculate line
      int x1,y1,x2,y2;
      
      x1=(int)(x-linelength*PDEvars[first_grad_layer][x][y]);
      y1=(int)(y-linelength*PDEvars[first_grad_layer+1][x][y]);
      x2=(int)(x+linelength*PDEvars[first_grad_layer][x][y]);
      y2=(int)(y+linelength*PDEvars[first_grad_layer+1][x][y]);
      if (x1<0) x1=0;
      if (x1>sizex-1) x1=sizex-1;
      if (y1<0) y1=0;
      if (y1>sizey-1) y1=sizey-1;
      if (x2<0) x2=0;
      if (x2>sizex-1) x2=sizex-1;
      if (y2<0) y2=0;
      if (y2>sizey-1) y2=sizey-1;

      // And draw it :-)
      // perhaps I can add arrowheads later to make it even nicer :-)
      g.Line(x1,y1,x2,y2,1);
    }
  }
}

bool PDE::plotPos(int x, int y, Graphics * graphics, int layer) {
  layer = 0;
  double val = PDEvars[layer][x][y];
  if (val > -100){
    graphics->Rectangle(MapColour(val), x, y); 
    return false;
  }
  return true;
}


void PDE::SetSpeciesName(int l, const char *name) {
    species_names[l]=string(name);
}


void PDE::InitLinearYGradient(int spec, double conc_top, double conc_bottom) {
    for (int y=0;y<sizey;y++) {
      double val=(double)conc_top+y*((double)(conc_bottom-conc_top)/(double)sizey);
    for (int x=0;x<sizex;x++) {
      PDEvars[spec][x][y]=val;
    }
    cerr << y << " " << val << endl;
  }
}
