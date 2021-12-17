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
#include <cusparse.h>

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
  btype = 1;
  dx2 = par.dx*par.dx;
  dt = par.dt;
  usePDEorAltPDE = true;
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
  free(couplingcoefficient);
  cudaFree(PDEvars);
  cudaFree(alt_PDEvars);
  cudaFree(couplingcoefficient);
  cudaFree(upperH);
  cudaFree(diagH);
  cudaFree(lowerH);
  cudaFree(BH);
  cudaFree(XH);
  cudaFree(upperV);
  cudaFree(diagV);
  cudaFree(lowerV);
  cudaFree(BV);
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

void PDE::AllocateTridiagonalvars(int sx, int sy){
  
  sizex=sx; sizey=sy;


  //Allocate lowerH
  lowerH=(PDEFIELD_TYPE **)malloc(sizey*sizeof(PDEFIELD_TYPE*));
  if (lowerH==NULL)
    MemoryWarning();

  lowerH[0]=(PDEFIELD_TYPE *)malloc(sizey*sizex*sizeof(PDEFIELD_TYPE));
  if (lowerH[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizey;i++) 
    lowerH[i]=lowerH[i-1]+sizex;}
  /* clear matrix */

  {for (int i=0;i<sizey*sizex;i++) 
    lowerH[0][i]=0; }

  //Allocate upperH
  upperH=(PDEFIELD_TYPE **)malloc(sizey*sizeof(PDEFIELD_TYPE*));
  if (upperH==NULL)
    MemoryWarning();

  upperH[0]=(PDEFIELD_TYPE *)malloc(sizey*sizex*sizeof(PDEFIELD_TYPE));
  if (upperH[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizey;i++) 
    upperH[i]=upperH[i-1]+sizex;}
  /* clear matrix */

  {for (int i=0;i<sizey*sizex;i++) 
    upperH[0][i]=0; }


  //Allocate diagH
  diagH=(PDEFIELD_TYPE **)malloc(sizey*sizeof(PDEFIELD_TYPE*));
  if (diagH==NULL)
    MemoryWarning();

  diagH[0]=(PDEFIELD_TYPE *)malloc(sizey*sizex*sizeof(PDEFIELD_TYPE));
  if (diagH[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizey;i++) 
    diagH[i]=diagH[i-1]+sizex;}
  /* clear matrix */

  {for (int i=0;i<sizey*sizex;i++) 
    diagH[0][i]=0; }


  //Allocate BH
  BH=(PDEFIELD_TYPE **)malloc(sizey*sizeof(PDEFIELD_TYPE*));
  if (BH==NULL)
    MemoryWarning();

  BH[0]=(PDEFIELD_TYPE *)malloc(sizey*sizex*sizeof(PDEFIELD_TYPE));
  if (BH[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizey;i++) 
    BH[i]=BH[i-1]+sizex;}

  /* clear matrix */
  {for (int i=0;i<sizey*sizex;i++) 
    BH[0][i]=0; }

    //Allocate XH
  XH=(PDEFIELD_TYPE **)malloc(sizey*sizeof(PDEFIELD_TYPE*));
  if (BH==NULL)
    MemoryWarning();

  XH[0]=(PDEFIELD_TYPE *)malloc(sizey*sizex*sizeof(PDEFIELD_TYPE));
  if (XH[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizey;i++) 
    XH[i]=XH[i-1]+sizex;}

  /* clear matrix */
  {for (int i=0;i<sizey*sizex;i++) 
    XH[0][i]=0; }


  //Allocate lowerV
  lowerV=(PDEFIELD_TYPE **)malloc(sizex*sizeof(PDEFIELD_TYPE*));
  if (lowerV==NULL)
    MemoryWarning();

  lowerV[0]=(PDEFIELD_TYPE *)malloc(sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (lowerV[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizex;i++) 
    lowerV[i]=lowerV[i-1]+sizey;}
  /* clear matrix */

  {for (int i=0;i<sizex*sizey;i++) 
    lowerV[0][i]=0; }


  //Allocate upperV
  upperV=(PDEFIELD_TYPE **)malloc(sizex*sizeof(PDEFIELD_TYPE*));
  if (upperV==NULL)
    MemoryWarning();

  upperV[0]=(PDEFIELD_TYPE *)malloc(sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (upperV[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizex;i++) 
    upperV[i]=upperV[i-1]+sizey;}
  /* clear matrix */

  {for (int i=0;i<sizex*sizey;i++) 
    upperV[0][i]=0; }


  //Allocate diagV
  diagV=(PDEFIELD_TYPE **)malloc(sizex*sizeof(PDEFIELD_TYPE*));
  if (diagV==NULL)
    MemoryWarning();

  diagV[0]=(PDEFIELD_TYPE *)malloc(sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (diagV[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizex;i++) 
    diagV[i]=diagV[i-1]+sizey;}
  /* clear matrix */

  {for (int i=0;i<sizex*sizey;i++) 
    diagV[0][i]=0; }


  //Allocate BV
  BV=(PDEFIELD_TYPE **)malloc(sizex*sizeof(PDEFIELD_TYPE*));
  if (BV==NULL)
    MemoryWarning();

  BV[0]=(PDEFIELD_TYPE *)malloc(sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (BV[0]==NULL)  
    MemoryWarning();

    {for (int i=1;i<sizex;i++) 
    BV[i]=BV[i-1]+sizey;}
  /* clear matrix */

  {for (int i=0;i<sizex*sizey;i++) 
    BV[0][i]=0; }

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
  PDEFIELD_TYPE decay_rate = (PDEFIELD_TYPE) * par.decay_rate;
  PDEFIELD_TYPE secr_rate = (PDEFIELD_TYPE) * par.secr_rate;
  
  
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

void PDE::InitializeCuda(CellularPotts *cpm){
  AllocateTridiagonalvars(sizex, sizey);
  cout << "A" << endl;
  couplingcoefficient = cpm->getCouplingCoefficient();
  cudaMallocManaged(&PDEvars[0][0], layers*sizex*sizey*sizeof(PDEFIELD_TYPE));

  cudaMallocManaged(&alt_PDEvars[0][0], layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&couplingcoefficient[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  cout << "B" << endl;

  //Needed for ADI steps
  cudaMallocManaged(&upperH[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&diagH[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&lowerH[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  /*cout << "BH[0] = " << BH[0] << endl;
  cout << "BH[0]+sizex =" << BH[0]+sizex << endl;
  cout << "BH[1] = " << BH[1] << endl; */
  cudaMallocManaged(&BH[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&XH[0], sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&upperV[0], sizey*sizex*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&diagV[0], sizey*sizex*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&lowerV[0], sizey*sizex*sizeof(PDEFIELD_TYPE));
  cudaMallocManaged(&BV[0], sizey*sizex*sizeof(PDEFIELD_TYPE));

  handleH = 0;
  pbuffersizeH = 0;
  pbufferH = NULL;
  statusH=cusparseCreate(&handleH);
  cusparseSgtsvInterleavedBatch_bufferSizeExt(handleH, 0, sizex, lowerH[0], diagH[0], upperH[0], BH[0], sizey, &pbuffersizeH); //Compute required buffersize for horizontal sweep

  cudaMalloc( &pbufferH, sizeof(char)* pbuffersizeH);

  handleV = 0;
  pbuffersizeV = 0;
  pbufferV = NULL;
  statusV=cusparseCreate(&handleV);
  cusparseSgtsvInterleavedBatch_bufferSizeExt(handleV, 0, sizey, lowerV[0], diagV[0], upperV[0], BV[0], sizex, &pbuffersizeV); //Compute required buffersize for vertical sweep
  cudaMalloc( &pbufferV, sizeof(char)* pbuffersizeV);
}

void PDE::InitializePDEs(CellularPotts *cpm){
  InitializePDEvars();
  InitializeCuda(cpm);
}


__global__ void InitializeDiagonals(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* lowerH, PDEFIELD_TYPE* upperH, PDEFIELD_TYPE* diagH, PDEFIELD_TYPE* lowerV, PDEFIELD_TYPE* upperV, PDEFIELD_TYPE* diagV, PDEFIELD_TYPE* couplingcoefficient){
  //This function could in theory be parellelized further, split into 6 (each part only assigning 1 value.), but this is probably slower
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int xloc; //position we currently want to assign to
  int yloc;
  int idcc; //id corresponding to the couplingcoefficient (+sizey to get the value right, +1 to get the value above)
  for (int id = index; id < sizex*sizey; id += stride){
    xloc = id%sizex;
    yloc = id/sizex;
    idcc = xloc*sizey + yloc;
    if(xloc == 0){
      lowerH[id] = 0;
      diagH[id] = couplingcoefficient[idcc+sizey]/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+sizey]/dx2;
    }
    else if(xloc == sizex -1){
      lowerH[id] = -couplingcoefficient[idcc-sizey]/dx2;;
      diagH[id] = couplingcoefficient[idcc-sizey]/dx2 + twooverdt;
      upperH[id] = 0;
    }
    else{
      lowerH[id] = -couplingcoefficient[idcc-sizey]/dx2;
      diagH[id] = (couplingcoefficient[idcc+sizey]+couplingcoefficient[idcc-sizey])/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+sizey]/dx2;
    }
    xloc = id/sizey;
    yloc = id%sizey;
    idcc = xloc*sizey + yloc;

    if(yloc == 0){
      lowerH[id] = 0;
      diagH[id] = couplingcoefficient[idcc+1]/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+1]/dx2;
    }
    else if(yloc == sizey -1){
      lowerH[id] = -couplingcoefficient[idcc-1]/dx2;;
      diagH[id] = couplingcoefficient[idcc-1]/dx2 + twooverdt;
      upperH[id] = 0;
    }
    else{
      lowerH[id] = -couplingcoefficient[idcc-1]/dx2;
      diagH[id] = (couplingcoefficient[idcc+1]+couplingcoefficient[idcc-1])/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+1]/dx2;
      
    }  

  
  }

}

__global__ void InitializeHorizontalVectors(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* couplingcoefficient, PDEFIELD_TYPE* alt_PDEvars){
  //This function could in theory be parellelized further, split into 6 (each part only assigning 1 value.), but this is probably slower
  
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int xloc;
  int yloc;
  int idcc; //id corresponding to the couplingcoefficient and alt_PDEvars(+sizey to get the value right, +1 to get the value above)
  for (int id = index; id < sizex*sizey; id += stride){
    xloc = id%sizex;
    yloc = id/sizex;
    idcc = xloc*sizey + yloc;

    if (yloc == 0)
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+1]*(alt_PDEvars[idcc+1] - alt_PDEvars[idcc]))/dx2; 
    else if (yloc == sizey-1)
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc-1]*(alt_PDEvars[idcc-1] - alt_PDEvars[idcc]))/dx2;  
    else 
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+1]*(alt_PDEvars[idcc+1] - alt_PDEvars[idcc]) + couplingcoefficient[idcc-1]*(alt_PDEvars[idcc-1] - alt_PDEvars[idcc]))/dx2;
  
  }

}

__global__ void InitializeVerticalVectors(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* couplingcoefficient, PDEFIELD_TYPE* alt_PDEvars){

  //This function could in theory be parellelized further, split into 6 (each part only assigning 1 value.), but this is probably slower
  
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int xloc;
  int yloc;
  int idcc; //id corresponding to the couplingcoefficient and alt_PDEvars(+sizey to get the value right, +1 to get the value above)
  for (int id = index; id < sizex*sizey; id += stride){
    xloc = id/sizey;
    yloc = id%sizey;
    idcc = xloc*sizey + yloc;

    if (xloc == 0)
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+sizey]*(alt_PDEvars[idcc+sizey] - alt_PDEvars[idcc]))/dx2; 
    else if (xloc == sizey-1)
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc-sizey]*(alt_PDEvars[idcc-sizey] - alt_PDEvars[idcc]))/dx2;  
    else 
      BH[idcc] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+sizey]*(alt_PDEvars[idcc+sizey] - alt_PDEvars[idcc]) + couplingcoefficient[idcc-sizey]*(alt_PDEvars[idcc-sizey] - alt_PDEvars[idcc]))/dx2;
  
  }

}

__global__ void NewPDEfieldH(int sizex, int sizey, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    PDEvars[(id%sizex)*sizey + id/sizex] = BH[id]; //Conversion is needed because PDEvars iterates over rowas first and then columns, while BH does the opposite
  }
}

__global__ void NewPDEfieldV(int sizex, int sizey, PDEFIELD_TYPE* BV, PDEFIELD_TYPE* PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    PDEvars[id] = BV[id];
  }
}

__device__ void ComputeDerivs(PDEFIELD_TYPE t, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, int id){ // computes the derivatives at some time point with forward Euler
    //equations for Paci2020
  
    // Software implementation of the Paci2020 model of the action potential 
    // of human induced pluripotent stem cell-derived cardiomyocytes, 
    // used in 10.1016/j.bpj.2020.03.018
    //
    // This software is provided for NON-COMMERCIAL USE ONLY 
    // (read the license included in the zip file).
  
  
    //-------------------------------------------------------------------------------
    // State variables
    //-------------------------------------------------------------------------------
  
    // 0: Vm (volt) (in Membrane)
    // 1: Ca_SR (millimolar) (in calcium_dynamics)
    // 2: Cai (millimolar) (in calcium_dynamics)
    // NOT USED 3: g (dimensionless) (in calcium_dynamics)
    // 4: d (dimensionless) (in i_CaL_d_gate)
    // 5: f1 (dimensionless) (in i_CaL_f1_gate)
    // 6: f2 (dimensionless) (in i_CaL_f2_gate)
    // 7: fCa (dimensionless) (in i_CaL_fCa_gate)
    // 8: Xr1 (dimensionless) (in i_Kr_Xr1_gate)
    // 9: Xr2 (dimensionless) (in i_Kr_Xr2_gate)
    // 10: Xs (dimensionless) (in i_Ks_Xs_gate)
    // 11: h (dimensionless) (in i_Na_h_gate)
    // 12: j (dimensionless) (in i_Na_j_gate)
    // 13: m (dimensionless) (in i_Na_m_gate)
    // 14: Xf (dimensionless) (in i_f_Xf_gate)
    // 15: q (dimensionless) (in i_to_q_gate)
    // 16: r (dimensionless) (in i_to_r_gate)
    // 17: Nai (millimolar) (in sodium_dynamics)
    // 18: m_L (dimensionless) (in i_NaL_m_gate)
    // 19: h_L (dimensionless) (in i_NaL_h_gate)
    // 20: RyRa (dimensionless) (in calcium_dynamics)
    // 21: RyRo (dimensionless) (in calcium_dynamics)
    // 22: RyRc (dimensionless) (in calcium_dynamics)
  
    //// Constants
    PDEFIELD_TYPE F = 96485.3415;     // coulomb_per_mole (in model_parameters)
    PDEFIELD_TYPE R = 8.314472;       // joule_per_mole_kelvin (in model_parameters)
    PDEFIELD_TYPE T = 310.0;          // kelvin (in model_parameters) //37°C
  
    //// Cell geometry
    PDEFIELD_TYPE V_SR = 583.73;        // micrometre_cube (in model_parameters)
    PDEFIELD_TYPE Vc   = 8800.0;        // micrometre_cube (in model_parameters)
    PDEFIELD_TYPE Cm   = 9.87109e-11;   // farad (in model_parameters)
  
    //// Extracellular concentrations
    PDEFIELD_TYPE Nao = 151.0; // millimolar (in model_parameters)
    PDEFIELD_TYPE Ko  = 5.4;   // millimolar (in model_parameters)
    PDEFIELD_TYPE Cao = 1.8;   // millimolar (in model_parameters)
  
    //// Intracellular concentrations
    // Naio = 10 mM y[17]
    PDEFIELD_TYPE Ki = 150.0;   // millimolar (in model_parameters)
    // Cai  = 0.0002 mM y[2]
    // caSR = 0.3 mM y[1]
  
    // time (second)
  
    //// Nernst potential
    PDEFIELD_TYPE E_Na = R*T/F*log(Nao/y[17]);
    PDEFIELD_TYPE E_Ca = 0.5*R*T/F*log(Cao/y[2]);
    PDEFIELD_TYPE E_K  = R*T/F*log(Ko/Ki);
    PDEFIELD_TYPE PkNa = 0.03;   // dimensionless (in electric_potentials)
    PDEFIELD_TYPE E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*y[17]));
  
    //// INa adapted from DOI:10.3389/fphys.2018.00080
    PDEFIELD_TYPE g_Na        = 3671.2302; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*6447.1896;
    PDEFIELD_TYPE i_Na        =  g_Na*pow((float)y[13],3.0f)*y[11]*y[12]*(y[0] - E_Na);
  
    PDEFIELD_TYPE m_inf       = 1 / (1 + exp((y[0]*1000 + 39)/-11.2));
    PDEFIELD_TYPE tau_m       = (0.00001 + 0.00013*exp(-pow((float)((y[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((y[0]*1000 + 42)/-5)));
    dydt[13]   = (m_inf-y[13])/tau_m;
  
    PDEFIELD_TYPE h_inf       = 1 / (1 + exp((y[0]*1000 + 66.5)/6.8));
    PDEFIELD_TYPE tau_h       = (0.00007 + 0.034 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(y[0]*1000 + 79)/14)));
    dydt[11]   = (h_inf-y[11])/tau_h;
  
    PDEFIELD_TYPE j_inf       = h_inf;
    PDEFIELD_TYPE tau_j       = 10*(0.0007 + 0.15 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(y[0]*1000 + 79)/14)));
    dydt[12]   = (j_inf-y[12])/tau_j;
  
  
    //// INaL
    PDEFIELD_TYPE myCoefTauM  = 1;
    PDEFIELD_TYPE tauINaL     = 200; //ms
    PDEFIELD_TYPE GNaLmax     = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
    PDEFIELD_TYPE Vh_hLate    = 87.61;
    PDEFIELD_TYPE i_NaL       = GNaLmax* pow((float)y[18],3.0f)*y[19]*(y[0]-E_Na);
  
    PDEFIELD_TYPE m_inf_L     = 1/(1+exp(-(y[0]*1000+42.85)/(5.264)));
    PDEFIELD_TYPE alpha_m_L   = 1/(1+exp((-60-y[0]*1000)/5));
    PDEFIELD_TYPE beta_m_L    = 0.1/(1+exp((y[0]*1000+35)/5))+0.1/(1+exp((y[0]*1000-50)/200));
    PDEFIELD_TYPE tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
    dydt[18]   = (m_inf_L-y[18])/tau_m_L*1000;
  
    PDEFIELD_TYPE h_inf_L     = 1/(1+exp((y[0]*1000+Vh_hLate)/(7.488)));
    PDEFIELD_TYPE tau_h_L     = 1 * tauINaL;
    dydt[19]   = (h_inf_L-y[19])/tau_h_L*1000;
  
    //// If adapted from DOI:10.3389/fphys.2018.00080
    PDEFIELD_TYPE g_f         = 1; //((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*22.2763088;
    PDEFIELD_TYPE fNa         = 0.37;
    PDEFIELD_TYPE fK          = 1 - fNa;
    PDEFIELD_TYPE i_fK        = fK*g_f*y[14]*(y[0] - E_K);
    PDEFIELD_TYPE i_fNa       = fNa*g_f*y[14]*(y[0] - E_Na);
    PDEFIELD_TYPE i_f         = i_fK + i_fNa;
  
    PDEFIELD_TYPE Xf_infinity = 1.0/(1.0 + exp((y[0]*1000 + 69)/8));
    PDEFIELD_TYPE tau_Xf      = 5600 / (1 + exp((y[0]*1000 + 65)/7) + exp(-(y[0]*1000 + 65)/19));
    dydt[14]   = 1000*(Xf_infinity-y[14])/tau_Xf;
    
  
    //// ICaL
    PDEFIELD_TYPE g_CaL       = 8.635702e-5;   // metre_cube_per_F_per_s (in i_CaL)
    PDEFIELD_TYPE i_CaL;  
    PDEFIELD_TYPE precision = 0.0001;     
    if(y[0]< precision && y[0] > -precision) //hopital
      i_CaL =  g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T) *y[4]*y[5]*y[6]*y[7] / (2.0*F/(R*T)) * (y[2] - 0.341*Cao);
    else
      i_CaL = g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*y[0]*pow(F,2.0)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7];
  
    PDEFIELD_TYPE d_infinity  = 1.0/(1.0+exp(-(y[0]*1000.0+9.1)/7.0));
    PDEFIELD_TYPE alpha_d     = 0.25+1.4/(1.0+exp((-y[0]*1000.0-35.0)/13.0));
    PDEFIELD_TYPE beta_d      = 1.4/(1.0+exp((y[0]*1000.0+5.0)/5.0));
    PDEFIELD_TYPE gamma_d     = 1.0/(1.0+exp((-y[0]*1000.0+50.0)/20.0));
    PDEFIELD_TYPE tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
    dydt[4]    = (d_infinity-y[4])/tau_d;
  
    PDEFIELD_TYPE f1_inf      = 1.0/(1.0+exp((y[0]*1000.0+26.0)/3.0));
    PDEFIELD_TYPE constf1;
    if (f1_inf-y[5] > 0.0)
        constf1 = 1.0+1433.0*(y[2]-50.0*1.0e-6);
    else
        constf1 = 1.0;
  
    PDEFIELD_TYPE tau_f1      = (20.0+1102.5*exp(-pow((float)((y[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-y[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf1/1000.0;
    dydt[5]    = (f1_inf-y[5])/tau_f1;
  
    PDEFIELD_TYPE f2_inf      = 0.33+0.67/(1.0+exp((y[0]*1000.0+32.0)/4.0));
    PDEFIELD_TYPE constf2     = 1.0;
    PDEFIELD_TYPE tau_f2      = (600.0*exp(-pow((float)(y[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-y[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf2/1000.0;
    dydt[6]    = (f2_inf-y[6])/tau_f2;
  
    PDEFIELD_TYPE alpha_fCa   = 1.0/(1.0+pow((float)(y[2]/0.0006),8.0f));
    PDEFIELD_TYPE beta_fCa    = 0.1/(1.0+exp((y[2]-0.0009)/0.0001));
    PDEFIELD_TYPE gamma_fCa   = 0.3/(1.0+exp((y[2]-0.00075)/0.0008));
    PDEFIELD_TYPE fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
  
    PDEFIELD_TYPE constfCa;
    if ((y[0] > -0.06) && (fCa_inf > y[7]))
        constfCa = 0.0;
    else
        constfCa = 1.0;
  
    PDEFIELD_TYPE tau_fCa     = 0.002;   // second (in i_CaL_fCa_gate)
    dydt[7]    = constfCa*(fCa_inf-y[7])/tau_fCa;
  
    //// Ito
    PDEFIELD_TYPE g_to        = 29.9038;   // S_per_F (in i_to)  
    PDEFIELD_TYPE i_to        = g_to*(y[0]-E_K)*y[15]*y[16];
  
    PDEFIELD_TYPE q_inf       = 1.0/(1.0+exp((y[0]*1000.0+53.0)/13.0));
    PDEFIELD_TYPE tau_q       = (6.06+39.102/(0.57*exp(-0.08*(y[0]*1000.0+44.0))+0.065*exp(0.1*(y[0]*1000.0+45.93))))/1000.0;
    dydt[15]   = (q_inf-y[15])/tau_q;
  
  
    PDEFIELD_TYPE r_inf       = 1.0/(1.0+exp(-(y[0]*1000.0-22.3)/18.75));
    PDEFIELD_TYPE tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(y[0]*1000.0+30.61))+0.369*exp(-0.12*(y[0]*1000.0+23.84))))/1000.0;
    dydt[16]   = (r_inf-y[16])/tau_r;
  
    //// IKs
    PDEFIELD_TYPE g_Ks        = 2.041;   // S_per_F (in i_Ks)
    PDEFIELD_TYPE i_Ks        = g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
  
    PDEFIELD_TYPE Xs_infinity = 1.0/(1.0+exp((-y[0]*1000.0-20.0)/16.0));
    PDEFIELD_TYPE alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-y[0]*1000.0)/6.0));
    PDEFIELD_TYPE beta_Xs     = 1.0/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
    PDEFIELD_TYPE tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
    dydt[10]   = (Xs_infinity-y[10])/tau_Xs;
  
    //// IKr
    PDEFIELD_TYPE L0           = 0.025;   // dimensionless (in i_Kr_Xr1_gate)
    PDEFIELD_TYPE Q            = 2.3;     // dimensionless (in i_Kr_Xr1_gate)
    PDEFIELD_TYPE g_Kr         = 29.8667;   // S_per_F (in i_Kr)
    PDEFIELD_TYPE i_Kr         = g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4);
  
    PDEFIELD_TYPE V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(L0*pow((float)(1.0+Cao/0.58),4.0f)))-0.019);
  
    PDEFIELD_TYPE Xr1_inf      = 1.0/(1.0+exp((V_half-y[0]*1000.0)/4.9));
    PDEFIELD_TYPE alpha_Xr1    = 450.0/(1.0+exp((-45.0-y[0]*1000.0)/10.0));
    PDEFIELD_TYPE beta_Xr1     = 6.0/(1.0+exp((30.0+y[0]*1000.0)/11.5));
    PDEFIELD_TYPE tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
    dydt[8]     = (Xr1_inf-y[8])/tau_Xr1;
  
    PDEFIELD_TYPE Xr2_infinity = 1.0/(1.0+exp((y[0]*1000.0+88.0)/50.0));
    PDEFIELD_TYPE alpha_Xr2    = 3.0/(1.0+exp((-60.0-y[0]*1000.0)/20.0));
    PDEFIELD_TYPE beta_Xr2     = 1.12/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
    PDEFIELD_TYPE tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
    dydt[9]    = (Xr2_infinity-y[9])/tau_Xr2;
  
    //// IK1
    PDEFIELD_TYPE alpha_K1    = 3.91/(1.0+exp(0.5942*(y[0]*1000.0-E_K*1000.0-200.0)));
    PDEFIELD_TYPE beta_K1     = (-1.509*exp(0.0002*(y[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(y[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(y[0]*1000.0-E_K*1000.0)));
    PDEFIELD_TYPE XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
    PDEFIELD_TYPE g_K1        = 28.1492;   // S_per_F (in i_K1)
    PDEFIELD_TYPE i_K1        = g_K1*XK1_inf*(y[0]-E_K)*sqrt(Ko/5.4);
  
    //// INaCa
    PDEFIELD_TYPE KmCa        = 1.38;   // millimolar (in i_NaCa)
    PDEFIELD_TYPE KmNai       = 87.5;   // millimolar (in i_NaCa)
    PDEFIELD_TYPE Ksat        = 0.1;    // dimensionless (in i_NaCa)
    PDEFIELD_TYPE gamma       = 0.35;   // dimensionless (in i_NaCa)
    PDEFIELD_TYPE alpha       = 2.16659;
    PDEFIELD_TYPE kNaCa      = 3917.0463; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaCaRedMed) * 6514.47574;   // A_per_F (in i_NaCa)
    PDEFIELD_TYPE i_NaCa      = kNaCa*(exp(gamma*y[0]*F/(R*T))*pow((float)y[17],3.0f)*Cao-exp((gamma-1.0)*y[0]*F/(R*T))*pow((float)Nao,3.0f)*y[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*y[0]*F/(R*T))));
  
    //// INaK
    PDEFIELD_TYPE Km_K        = 1.0;    // millimolar (in i_NaK)
    PDEFIELD_TYPE Km_Na       = 40.0;   // millimolar (in i_NaK)
    PDEFIELD_TYPE PNaK        = 2.74240;// A_per_F (in i_NaK)
    PDEFIELD_TYPE i_NaK       = PNaK*Ko/(Ko+Km_K)*y[17]/(y[17]+Km_Na)/(1.0+0.1245*exp(-0.1*y[0]*F/(R*T))+0.0353*exp(-y[0]*F/(R*T)));
  
    //// IpCa
    PDEFIELD_TYPE KPCa        = 0.0005;   // millimolar (in i_PCa)
    PDEFIELD_TYPE g_PCa       = 0.4125;   // A_per_F (in i_PCa)
    PDEFIELD_TYPE i_PCa       = g_PCa*y[2]/(y[2]+KPCa);
  
    //// Background currents
    PDEFIELD_TYPE g_b_Na      = 1.14;         // S_per_F (in i_b_Na)
    PDEFIELD_TYPE i_b_Na      = g_b_Na*(y[0]-E_Na);
  
    PDEFIELD_TYPE g_b_Ca      = 0.8727264;    // S_per_F (in i_b_Ca)
    PDEFIELD_TYPE i_b_Ca      = g_b_Ca*(y[0]-E_Ca);
  
    //// Sarcoplasmic reticulum
    PDEFIELD_TYPE VmaxUp		= 0.82205;
    PDEFIELD_TYPE Kup			=  4.40435e-4;
    PDEFIELD_TYPE i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)y[2],2.0f));
  
    PDEFIELD_TYPE V_leak		= 4.48209e-4;
    PDEFIELD_TYPE i_leak      = (y[1]-y[2])*V_leak;
  
    dydt[3]    = 0;
  
    // RyR
    PDEFIELD_TYPE g_irel_max	= 55.808061;
    PDEFIELD_TYPE RyRa1       = 0.05169;
    PDEFIELD_TYPE RyRa2       = 0.050001;
    PDEFIELD_TYPE RyRahalf    = 0.02632;
    PDEFIELD_TYPE RyRohalf    = 0.00944;
    PDEFIELD_TYPE RyRchalf    = 0.00167;
  
    PDEFIELD_TYPE RyRSRCass   = (1 - 1/(1 +  exp((y[1]-0.3)/0.1)));
    PDEFIELD_TYPE i_rel       = g_irel_max*RyRSRCass*y[21]*y[22]*(y[1]-y[2]);
  
    PDEFIELD_TYPE RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*y[2]-(RyRahalf))/0.0082));
    PDEFIELD_TYPE RyRtauadapt = 1; //s
    dydt[20]   = (RyRainfss- y[20])/RyRtauadapt;
  
    PDEFIELD_TYPE RyRoinfss   = (1 - 1/(1 +  exp((1000*y[2]-(y[20]+ RyRohalf))/0.003)));
    PDEFIELD_TYPE RyRtauact;
    if (RyRoinfss>= y[21])
      RyRtauact = 18.75e-3;       //s
    else
      RyRtauact = 0.1*18.75e-3;   //s
  
    dydt[21]    = (RyRoinfss- y[21])/(RyRtauact);
  
    PDEFIELD_TYPE RyRcinfss   = (1/(1 + exp((1000*y[2]-(y[20]+RyRchalf))/0.001)));
    PDEFIELD_TYPE RyRtauinact;
    if (RyRcinfss>= y[22])
      RyRtauinact = 2*87.5e-3;    //s
    else
      RyRtauinact = 87.5e-3;      //s
  
    dydt[22]    = (RyRcinfss- y[22])/(RyRtauinact);
  
  
  
  
    //// Ca2+ buffering
    PDEFIELD_TYPE Buf_C       = 0.25;   // millimolar (in calcium_dynamics)
    PDEFIELD_TYPE Buf_SR      = 10.0;   // millimolar (in calcium_dynamics)
    PDEFIELD_TYPE Kbuf_C      = 0.001;   // millimolar (in calcium_dynamics)
    PDEFIELD_TYPE Kbuf_SR     = 0.3;   // millimolar (in calcium_dynamics)
    PDEFIELD_TYPE Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(y[2]+Kbuf_C), 2.0f));
    PDEFIELD_TYPE Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(y[1]+Kbuf_SR), 2.0f));
  
    //// Ionic concentrations
    //Nai
    dydt[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
    //Cai
    dydt[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
     //caSR
    dydt[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
  
    //// Stimulation
  //  PDEFIELD_TYPE i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
  //  PDEFIELD_TYPE i_stim_End 				= 1000.0;   // second (in stim_mode)
  //  PDEFIELD_TYPE i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
  //  PDEFIELD_TYPE i_stim_Start 			= 0.0;   // second (in stim_mode)
  //  PDEFIELD_TYPE i_stim_frequency        = 60.0;   // per_second (in stim_mode)
    //PDEFIELD_TYPE stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
  //  PDEFIELD_TYPE i_stim_Period 			= 60.0/i_stim_frequency;
  
    //if stim_flag~=0 && stim_flag~=1
    //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
    //end
  
    /*
    if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
        i_stim = stim_flag*i_stim_Amplitude/Cm;
    else
        i_stim = 0.0;
    */
    PDEFIELD_TYPE i_stim = 0;
  
    //// Membrane potential
    dydt[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);
}
__global__ void RungeKuttaStep(PDEFIELD_TYPE dt, PDEFIELD_TYPE thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int i;

  PDEFIELD_TYPE dydt[23];
  PDEFIELD_TYPE ak2[23];
  PDEFIELD_TYPE ak3[23];
  PDEFIELD_TYPE ak4[23];
  PDEFIELD_TYPE ak5[23];
  PDEFIELD_TYPE ak6[23];
  PDEFIELD_TYPE ytemp[23]; 

  PDEFIELD_TYPE y[23];

  PDEFIELD_TYPE a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0;
  for (int id = index; id < sizex*sizey; id += stride){
    for (int l = 0; l < layers; l++) //fill with current PDE values
      y[l] = PDEvars[l*sizex*sizey + l];
  


    //Paci2020
 
    ComputeDerivs(thetime,y,dydt,id);
    for (i=0;i<layers;i++) //First step.
      ytemp[i]=y[i]+b21*dt*dydt[i];
    ComputeDerivs(thetime+a2*dt,ytemp,ak2,id);// Second step.

    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b31*dydt[i]+b32*ak2[i]);
    ComputeDerivs(thetime+a3*dt,ytemp,ak3,id); //Third step.
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b41*dydt[i]+b42*ak2[i]+b43*ak3[i]);
    ComputeDerivs(thetime+a4*dt,ytemp,ak4,id); //Fourth step.
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b51*dydt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    ComputeDerivs(thetime+a5*dt,ytemp,ak5,id); //Fifth step.
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b61*dydt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    ComputeDerivs(thetime+a6*dt,ytemp,ak6,id); //Sixth step.
    for (i=0;i<layers;i++) //Accumulate increments with proper weights.
      alt_PDEvars[i*sizex*sizey + id]=PDEvars[i*sizex*sizey + id]+(c1*dydt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i])*dt;
  }
}

void PDE::cuPDEsteps(CellularPotts * cpm, int repeat){
  //setup matrices for upperdiagonal, diagonal and lower diagonal for both the horizontal and vertical direction, since these remain the same during once MCS
  InitializeDiagonals<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, lowerH[0], upperH[0], diagH[0], lowerV[0], upperV[0], diagV[0], couplingcoefficient[0]);
  for (int iteration = 0; iteration < repeat; iteration++){
    //Do an ODE step of size dt/2
    RungeKuttaStep<<<par.number_of_cores, par.threads_per_core>>>(dt, thetime, layers, sizex, sizey, PDEvars[0][0], alt_PDEvars[0][0]);
    
    //Do a horizontal ADI sweep of size dt/2
    InitializeHorizontalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BH[0], couplingcoefficient[0], alt_PDEvars[0][0]);
    statusH = cusparseSgtsvInterleavedBatch(handleH, 0, sizex, lowerH[0], diagH[0], upperH[0], BH[0], sizey, &pbuffersizeH);
    if (statusH != CUSPARSE_STATUS_SUCCESS)
    {
      cout << statusH << endl;
    }
    NewPDEfieldH<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, BH[0], PDEvars[0][0]);

    //Do an ODE step of size dt/2
    RungeKuttaStep<<<par.number_of_cores, par.threads_per_core>>>(dt, thetime, layers, sizex, sizey, PDEvars[0][0], alt_PDEvars[0][0]);

    //Do a vertical ADI sweep of size dt/2
    InitializeVerticalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BV[0], couplingcoefficient[0], alt_PDEvars[0][0]);
    statusV = cusparseSgtsvInterleavedBatch(handleV, 0, sizey, lowerV[0], diagV[0], upperV[0], BV[0], sizex, &pbuffersizeV);
    if (statusV != CUSPARSE_STATUS_SUCCESS)
    {
      cout << statusV << endl;
    }
    NewPDEfieldV<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, BV[0], PDEvars[0][0]); //////

    thetime = thetime + dt;
  }
}

// public
void PDE::Diffuse(int repeat) {
  
  // Just diffuse everywhere (cells are transparent), using finite difference
  // (We're ignoring the problem of how to cope with moving cell
  // boundaries right now)
  
  const PDEFIELD_TYPE dt=par.dt;
  const PDEFIELD_TYPE dx2=par.dx*par.dx;

  for (int r=0;r<repeat;r++) {
    //NoFluxBoundaries();
    if (par.periodic_boundaries) {
      PeriodicBoundaries();
    } else {
      AbsorbingBoundaries();
      //NoFluxBoundaries();
    }
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
