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

    //equations for Paci2020
  
    // Software implementation of the Paci2020 model of the action potential 
    // of human induced pluripotent stem cell-derived cardiomyocytes, 
    // used in 10.1016/j.bpj.2020.03.018
    //
    // This software is provided for NON-COMMERCIAL USE ONLY 
    // License included:

/*
ACADEMIC PUBLIC LICENSE (software implementation of the Paci2020 model, v1.0)

Copyright (c) 2021 Michelangelo Paci, Tampere University Foundation sr, All rights reserved

The following license governs the use of the software implementation of the Paci2020 model of the action potential of human induced pluripotent stem cell-derived cardiomyocytes, used in the paper "All-Optical Electrophysiology Refines Populations of In Silico Human iPSC-CMs for Drug Evaluation" (10.1016/j.bpj.2020.03.018), in non-commercial academic environments. In case of need for license extensions, please contact inventions@tuni.fi

------------------------------------------------------------------------------------------------------
Preamble

This license contains the terms and conditions of using the software implementation of the Paci2020 model in non-commercial settings: at academic institutions for teaching and research use, and at not-for-profit research organizations. You will find that this license provides non-commercial users of the software implementation of the Paci2020 model with rights that are similar to the well-known GNU General Public License 2.0, yet it retains the possibility for the software implementation of the Paci2020 model authors to financially support the development by selling commercial licenses. In fact, if you intend to use the software implementation of the Paci2020 model in a "for-profit" environment, where the software implementation of the Paci2020 model simulations are conducted to develop or enhance a product (including commercial or industry-sponsored research at academic institutions), or to use the software implementation of the Paci2020 model in a commercial service offering, then you need to obtain a license extension for the software implementation of the Paci2020 model. In that case, please contact inventions@tuni.fi.

What are the rights given to non-commercial users? Similarly, to GPL 2.0, you have the right to use the software, to distribute copies, to receive source code, to change the software and distribute your modifications or the modified software. Also, similarly to the GPL 2.0, if you distribute verbatim or modified copies of this software, they must be distributed under this license.

By modeling the GPL 2.0, this license guarantees that you’re safe when using the software implementation of the Paci2020 model in your work, for teaching, and research. This license guarantees that the software implementation of the Paci2020 model will remain available free of charge for non-profit use. You can modify the software implementation of the Paci2020 model to your purposes, and you can also share your modifications. Even in case of the authors abandoning the software implementation of the Paci2020 model entirely, this license permits anyone to continue developing it from the last release, and to create further releases under this license.

The precise terms and conditions for using, copying, distribution and modification follow.
------------------------------------------------------------------------------------------------------

Terms and Conditions for Use, Copying, Distribution and Modification
Definitions

•	"Program" means a copy of the software implementation of the Paci2020 model and all the files included in this archive, which are said to be distributed under this Academic Public License.
•	"Work based on the Program" means either the Program or any derivative work under copyright law: that is to say, a work containing the Program or a portion of it, either verbatim or with modifications and/or translated into another language. (Hereinafter, translation is included without limitation in the term "modification".)
•	"Using the Program" means any act of creating executables that contain or directly use libraries that are part of the Program, running any of the tools that are part of the Program, or creating works based on the Program.
•	Each licensee is addressed as "you".

§1. Permission is hereby granted to use the Program free of charge for any non-commercial purpose, including teaching and research at universities, colleges and other educational institutions, non-commercial research at organizations that are either not-for-profit or reinvest all profits in their scientific research, and personal not-for-profit purposes. For using the Program for commercial purposes, including but not restricted to commercial research at academic institutions, industrially sponsored research at academic institutions, consulting activities, and design of commercial hardware or software products or services, you have to contact inventions@tuni.fi for an appropriate license.

§2. You may copy and distribute verbatim copies of the source code of the program via any medium, provided that you add a conspicuous and appropriate copyright notice and a warranty disclaimer to each copy. Retain all notices relating to this license and the lack of any warranty. Forward a copy of this license to all other recipients of the program.

§3. You are entitled to change your copies of the Program or a part thereof and thus create a work based on the Program. You may copy and distribute changes or work in accordance with the provisions of Section 2 provided you also meet all of the following conditions: a) You must ensure that the changed files are provided with noticeable comments stating the author of the change and when this change was made. b) You must ensure that all work that you distribute or publish, that contains or is derived from the Program or parts thereof, as a whole, is licensed under the conditions of this license.

These requirements apply to the changed work as a whole. If identifiable sections of this work do not come from the Program and can be considered separate, this license and its terms do not apply to those sections if you distribute them as separate work. However, if you distribute the same sections as part of a whole that is based on the Program, the distribution of the whole must be done in accordance with the terms of this license as outlined in §2, independently of who wrote it.

The mere merging of another work that is not based on the Program with the Program (or a work based on the Program) does not bring the other work into the scope of this license.

§4. You may copy and distribute the Program (or a work based on it, in accordance with §3, in object code or executable form in accordance with the provisions of above Sections 2 and 3, provided that you also add the complete corresponding machine-readable source code. For an executable program, complete source code means the entire source code for all modules contained therein, as well as all associated interface definition files and scripts, with which the compilation and installation of the executable file is controlled.

§5. Any attempt to copy, modify, sublicense or distribute the Program in any other way than specified in this license is void, and will automatically terminate your rights under this license. However, parties who have received copies or rights from you under this license will not lose their license as long as these parties fully comply with the terms.

§6. You do not have to accept this license because you have not signed it. However, if you want to change or distribute the Program (or a work based on the Program), you automatically consent to this license and all its terms for copying, distributing or changing the program or the works based upon it.

§7. Each time you redistribute the Program (or any work based on the Program), the recipient automatically acquires a license from the initial licensor to copy, distribute or modify the Program in accordance with these terms and conditions. You may not impose any further restrictions on the recipient's exercise of the rights granted here. You are not responsible for ensuring that this license is enforced by third parties.

§8. If, as a result of a court decision or violation of a patent right, or for any other reason (not limited to patent issues), conditions are imposed that conflict with the terms of this license, you will not be released from the terms of this license. If you cannot distribute the Program because you would have to meet obligations under this license and other obligations at the same time, you may not distribute the Program at all.

§9. If the distribution and/or use of the Program in certain countries is restricted either by patents or by copyrighted interfaces, the original copyright holder who puts the Program under this license may add an explicit geographic distribution restriction that excludes these countries. In this case, this license contains the restriction as if it was written in the body of this license.

NO WARRANTY

§10. SINCE THE PROGRAM IS LICENSED FOR FREE, THERE IS NO GUARANTEE FOR THE PROGRAM. THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT ANY EXPRESSED OR IMPLIED GUARANTEE, INCLUDING BUT NOT LIMITED TO THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS YOURS. IF THE PROGRAM TURNS OUT TO BE FAULTY, YOU ARE RESPONSIBLE FOR THE COSTS FOR ALL NECESSARY MAINTENANCE, REPAIR OR CORRECTION WORK.

§11. UNDER NO CIRCUMSTANCES WILL A COPYRIGHT HOLDER OR ANY OTHER PARTY WHO CAN MODIFY AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE BE LIABLE FOR DAMAGE, INCLUDING GENERAL, SPECIAL, ACCIDENTAL OR OTHER DAMAGE. THE DISCLAIMER ALSO INCLUDES CONSEQUENTIAL DAMAGES THAT RESULT FROM USING THE PROGRAM ALONE OR IN CONJUNCTION WITH OTHER PROGRAMS, INCLUDING BUT NOT LIMITED TO THE LOSS OR CORRUPTION OF DATA.

IF YOU DO NOT AGREE WITH THIS LICENSE TERMS, DO NOT USE, COPY, CHANGE OR DISTRIBUTE THE PROGRAM (OR A WORK BASED ON THE PROGRAM).

This license was
- initially written by Andras Varga (public domain) for OMNeT++ https://omnetpp.org/intro/license, 
- then adapted by the openCARP project https://opencarp.org/download/license,
- now adapted for the software implementation of the Paci2020 model.
The adaptation is licensed under CC0 1.0 (Public Domain Dedication).
*/
  
  
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


/*ODE solver is an adaptation from Press, W. H. (2007). 
Numerical recipes : the art of scientific computing (3rd ed.). 
/New York, N.Y., [etc.]: Cambridge University Press.

*/




#include <stdio.h>
#include <unistd.h>
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
#include <cuda_profiler_api.h>
#define ARRAY_SIZE 2

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

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
  bpm = par.beats_per_minute;
  pacing_interval = 1/(bpm/60);
  PDEvars = new PDEFIELD_TYPE[layers*sizex*sizey];
  alt_PDEvars = new PDEFIELD_TYPE[layers*sizex*sizey];
  min_stepsize = par.min_stepsize;
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
    cudaFree(PDEvars);
  }
  if (alt_PDEvars) {
    cudaFree(alt_PDEvars);
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

/*
PDEFIELD_TYPE ***PDE::AllocatePDEvars(const int layers, const int sx, const int sy) { //Omschrijven naar [xcoordinaat][ycoordinaat][layercoordinaat] (eerst compileren zo)
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
    mem[i]=mem[i-1]+sizex*sizey;
  }
  mem[0][0]=(PDEFIELD_TYPE *)malloc(layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  if (mem[0][0]==NULL) {
    MemoryWarning();
  }
  for (int i=1;i<layers*sizex;i++) {
    mem[0][i]=mem[0][i-1]+sizey;
  }
  
  //Clear PDE plane 
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
  // clear matrix 

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
  // clear matrix 

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
  // clear matrix 

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

  // clear matrix 
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

  // clear matrix
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
  // clear matrix

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
  // clear matrix 

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
  // clear matrix

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
  // clear matrix

  {for (int i=0;i<sizex*sizey;i++) 
    BV[0][i]=0; }

}*/

void PDE::InitializePDEvars(CellularPotts *cpm){
  PDEFIELD_TYPE PDEinit[ARRAY_SIZE];
  bool* mask = cpm->getMask()[0];
  /* For Paci2018
  PDEinit[0] = -0.070;
  PDEinit[1] = 0.32;
  PDEinit[2] = 0.0002;
  PDEinit[3] = 0;
  PDEinit[4] = 0;
  PDEinit[5] = 1;
  PDEinit[6] = 1;
  PDEinit[7] = 1;
  PDEinit[8] = 0;
  PDEinit[9] = 1;
  PDEinit[10] = 0;
  PDEinit[11] = 0.75;
  PDEinit[12] = 0.75; 
  PDEinit[13] = 0;
  PDEinit[14] = 0.1;
  PDEinit[15] = 1;
  PDEinit[16] = 0;
  PDEinit[17] = 9.2;
  PDEinit[18] = 0;
  PDEinit[19] = 0.75;
  PDEinit[20] = 0.3;
  PDEinit[21] = 0.9;
  PDEinit[22] = 0.1;*/

  PDEinit[0] = -1.2275879383;
  PDEinit[1] = -0.6109462976;
  for (int layer = 0; layer<layers; layer++){
    for (int i = layer*sizex*sizey; i<(layer+1)*sizex*sizey; i++){
      PDEvars[i] = PDEinit[layer];
      //if (layer == 0 && (i%(sizey+1) == 130))
      //  PDEvars[i] = 100;
    }
  }

  for (int i = 0; i < sizex*sizey; i++)
    if(!mask[i])
          PDEvars[i] = -10;
  //PDEvars[246183] = 10000;

}




void PDE::Plot(Graphics *g,const int l) {
  // l=layer: default layer is 0
  for (int x=0;x<sizex;x++) {
    for (int y=0;y<sizey;y++) {
      // Make the pixel four times as large
      // to fit with the CPM plane
      g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x,y);
      g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x+1,y);
      g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x,y+1);
      g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x+1,y+1);
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
        g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x,y);
        g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x+1,y);
        g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x,y+1);
        g->Point(MapColour(PDEvars[l*sizex*sizey+x*sizey+y]),x+1,y+1);
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

  conrec(&PDEvars[l*sizex*sizey],0,sizex-1,0,sizey-1,x,y,nc,z,g,colour);
  
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
      clm.queue.enqueueWriteBuffer(clm.pdeA,  CL_TRUE, 0, sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars);
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
                            sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars);
    }
    else {
      clm.queue.enqueueReadBuffer(clm.pdeA,CL_TRUE,0,
                            sizeof(PDEFIELD_TYPE)*sizex*sizey*layers, PDEvars);
    }
    if (errorcode != CL_SUCCESS) cout << "error:" << errorcode << endl;
    
}

void PDE::InitializeCuda(CellularPotts *cpm){
  //AllocateTridiagonalvars(sizex, sizey);

  cudaMalloc((void**) &d_couplingcoefficient, sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMalloc((void**) &d_celltype, sizex*sizey*sizeof(int));

  cudaMalloc((void**) &d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMemcpy(d_PDEvars, PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  cudaMalloc((void**) &d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMemcpy(d_alt_PDEvars, alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);


  //Needed for ADI steps
  gpuErrchk(cudaMallocManaged(&upperH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&diagH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&lowerH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&BH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&XH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&upperV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&diagV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&lowerV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&BV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&next_stepsize, sizey*sizex*sizeof(PDEFIELD_TYPE)));

  handleH = 0;
  pbuffersizeH = 0;
  pbufferH = NULL;
  statusH=cusparseCreate(&handleH);
  cusparseSgtsvInterleavedBatch_bufferSizeExt(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, &pbuffersizeH); //Compute required buffersize for horizontal sweep

  gpuErrchk(cudaMalloc( &pbufferH, sizeof(char)* pbuffersizeH));
  

  handleV = 0;
  pbuffersizeV = 0;
  pbufferV = NULL;
  statusV=cusparseCreate(&handleV);
  cusparseSgtsvInterleavedBatch_bufferSizeExt(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, &pbuffersizeV); //Compute required buffersize for vertical sweep
  gpuErrchk(cudaMalloc( &pbufferV, sizeof(char)* pbuffersizeV));

  

}

__global__ void InitializeLastStepsize(PDEFIELD_TYPE min_stepsize, PDEFIELD_TYPE* next_stepsize, int sizex, int sizey){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    next_stepsize[id] = min_stepsize;
  }
}


void PDE::InitializePDEs(CellularPotts *cpm){
  InitializePDEvars(cpm);
  InitializeCuda(cpm);
  InitializeLastStepsize<<<par.number_of_cores, par.threads_per_core>>>(min_stepsize, next_stepsize, sizex, sizey);
  cudaDeviceSynchronize();
}




__global__ void InitializeDiagonals(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* lowerH, PDEFIELD_TYPE* upperH, PDEFIELD_TYPE* diagH, PDEFIELD_TYPE* lowerV, PDEFIELD_TYPE* upperV, PDEFIELD_TYPE* diagV, PDEFIELD_TYPE* couplingcoefficient){
  //This function could in theory be parellelized further, split into 6 (each part only assigning 1 value.), but this is probably slower
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int xloc; //position we currently want to assign to
  int yloc;
  int idcc; //id corresponding to the couplingcoefficient (+sizey to get the value right, +1 to get the value above)
  for (int id = index; id < sizex*sizey; id += stride){
    xloc = id/sizey; //needed to obtain interleaved format
    yloc = id%sizey;
    idcc = xloc*sizey + yloc;
    if(xloc == 0){
      lowerH[id] = 0;
      diagH[id] = couplingcoefficient[idcc+sizey]/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+sizey]/dx2;  
    }
    else if(xloc == sizex -1){
      lowerH[id] = -couplingcoefficient[idcc-sizey]/dx2;
      diagH[id] = couplingcoefficient[idcc-sizey]/dx2 + twooverdt;
      upperH[id] = 0;
    }
    else{
      lowerH[id] = -couplingcoefficient[idcc-sizey]/dx2;
      diagH[id] = (couplingcoefficient[idcc+sizey]+couplingcoefficient[idcc-sizey])/dx2 + twooverdt;
      upperH[id] = -couplingcoefficient[idcc+sizey]/dx2;
    }

    xloc = id%sizex; //needed to obtain interleaved format
    yloc = id/sizex;
    idcc = xloc*sizey + yloc;
    if(yloc == 0){
      lowerV[id] = 0;
      diagV[id] = couplingcoefficient[idcc+1]/dx2 + twooverdt;
      upperV[id] = -couplingcoefficient[idcc+1]/dx2;
    }
    else if(yloc == sizey -1){
      lowerV[id] = -couplingcoefficient[idcc-1]/dx2;
      diagV[id] = couplingcoefficient[idcc-1]/dx2 + twooverdt;
      upperV[id] = 0;
    }
    else{
      lowerV[id] = -couplingcoefficient[idcc-1]/dx2;
      diagV[id] = (couplingcoefficient[idcc+1]+couplingcoefficient[idcc-1])/dx2 + twooverdt;
      upperV[id] = -couplingcoefficient[idcc+1]/dx2;
      
    }
  }
}

__global__ void InitializeHorizontalVectors(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* couplingcoefficient, PDEFIELD_TYPE* alt_PDEvars){
  
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int xloc;
  int yloc;
  int idcc; //id corresponding to the couplingcoefficient and alt_PDEvars(+sizey to get the value right, +1 to get the value above)
  for (int id = index; id < sizex*sizey; id += stride){
    xloc = id/sizey; //needed to obtain interleaved format
    yloc = id%sizey;
    idcc = xloc*sizey + yloc;

    if (yloc == 0)
      BH[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+1]*(alt_PDEvars[idcc+1] - alt_PDEvars[idcc]))/dx2; 
    else if (yloc == sizey-1)
      BH[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc-1]*(alt_PDEvars[idcc-1] - alt_PDEvars[idcc]))/dx2;  
    else 
      BH[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+1]*(alt_PDEvars[idcc+1] - alt_PDEvars[idcc]) + couplingcoefficient[idcc-1]*(alt_PDEvars[idcc-1] - alt_PDEvars[idcc]))/dx2;
  }
}

__global__ void InitializeVerticalVectors(int sizex, int sizey, PDEFIELD_TYPE twooverdt, PDEFIELD_TYPE dx2, PDEFIELD_TYPE* BV, PDEFIELD_TYPE* couplingcoefficient, PDEFIELD_TYPE* alt_PDEvars){

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

    if (xloc == 0)
      BV[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+sizey]*(alt_PDEvars[idcc+sizey] - alt_PDEvars[idcc]))/dx2; 
    else if (xloc == sizex-1)
      BV[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc-sizey]*(alt_PDEvars[idcc-sizey] - alt_PDEvars[idcc]))/dx2;  
    else 
      BV[id] = twooverdt*alt_PDEvars[idcc] + (couplingcoefficient[idcc+sizey]*(alt_PDEvars[idcc+sizey] - alt_PDEvars[idcc]) + couplingcoefficient[idcc-sizey]*(alt_PDEvars[idcc-sizey] - alt_PDEvars[idcc]))/dx2;
  }

}



__global__ void NewPDEfieldH0(int sizex, int sizey, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* PDEvars){ //Take the values from BH and assign the new values of the first layers of PDEvars
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride)
    PDEvars[id] = BH[id];      
}


__global__ void NewPDEfieldV0(int sizex, int sizey, PDEFIELD_TYPE* BV, PDEFIELD_TYPE* PDEvars){ //Take the values from BV and assign the new values of the first layers of PDEvars
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    PDEvars[sizey*(id%sizex)+id/sizex] = BV[id]; //Conversion is needed because PDEvars iterates over columns first and then rows, while BV does the opposite 
  }
}


__global__ void NewPDEfieldOthers(int sizex, int sizey, int layers, PDEFIELD_TYPE* BH, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars){ //copy the other values from alt_PDEvars to PDEvars
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index+sizex*sizey; id < layers*sizex*sizey; id += stride)
    PDEvars[id] = alt_PDEvars[id]; 
}





#if 0
__global__ void RungeKuttaStepOld(PDEFIELD_TYPE dt, PDEFIELD_TYPE thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype){
  
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int i;



  
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE ak2[ARRAY_SIZE];
  PDEFIELD_TYPE ak3[ARRAY_SIZE];
  PDEFIELD_TYPE ak4[ARRAY_SIZE];
  PDEFIELD_TYPE ak5[ARRAY_SIZE];
  PDEFIELD_TYPE ak6[ARRAY_SIZE];
  PDEFIELD_TYPE ytemp[ARRAY_SIZE]; 
  PDEFIELD_TYPE y[ARRAY_SIZE];
  PDEFIELD_TYPE yerr[ARRAY_SIZE];
  PDEFIELD_TYPE yout[ARRAY_SIZE];

  static PDEFIELD_TYPE a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
  dc5 = -277.00/14336.0;
  PDEFIELD_TYPE dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;

  //Declare variables needed for Paci2020 model and assign the constants

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
    
      //// Nernst potential
  PDEFIELD_TYPE E_Na;
  PDEFIELD_TYPE E_Ca;
  PDEFIELD_TYPE E_K;
  PDEFIELD_TYPE PkNa = 0.03;   // dimensionless (in electric_potentials)
  PDEFIELD_TYPE E_Ks;
    
  //// INa adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_Na = 3671.2302; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*6447.1896;
  PDEFIELD_TYPE i_Na;
    
  PDEFIELD_TYPE m_inf;
  PDEFIELD_TYPE tau_m;
    
  PDEFIELD_TYPE h_inf;
  PDEFIELD_TYPE tau_h;
    
  PDEFIELD_TYPE j_inf;
  PDEFIELD_TYPE tau_j;
    
    
  //// INaL
  PDEFIELD_TYPE myCoefTauM  = 1;
  PDEFIELD_TYPE tauINaL = 200; //ms
  PDEFIELD_TYPE GNaLmax = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
  PDEFIELD_TYPE Vh_hLate = 87.61;
  PDEFIELD_TYPE i_NaL;
    
  PDEFIELD_TYPE m_inf_L;
  PDEFIELD_TYPE alpha_m_L;
  PDEFIELD_TYPE beta_m_L;
  PDEFIELD_TYPE tau_m_L;
    
  PDEFIELD_TYPE h_inf_L;
  PDEFIELD_TYPE tau_h_L = 1 * tauINaL;
    
  //// If adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_f = 1; //((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*22.2763088;
  PDEFIELD_TYPE fNa = 0.37;
  PDEFIELD_TYPE fK = 1 - fNa;
  PDEFIELD_TYPE i_fK;
  PDEFIELD_TYPE i_fNa;
  PDEFIELD_TYPE i_f;
    
  PDEFIELD_TYPE Xf_infinity;
  PDEFIELD_TYPE tau_Xf; 
    
      //// ICaL
  PDEFIELD_TYPE g_CaL = 8.635702e-5;   // metre_cube_per_F_per_s (in i_CaL)
  PDEFIELD_TYPE i_CaL;  
  PDEFIELD_TYPE precision = 0.0001;     
    
  PDEFIELD_TYPE d_infinity;
  PDEFIELD_TYPE alpha_d;
  PDEFIELD_TYPE beta_d;
  PDEFIELD_TYPE gamma_d;
  PDEFIELD_TYPE tau_d;
    
  PDEFIELD_TYPE f1_inf;
  PDEFIELD_TYPE constf1;
    
  PDEFIELD_TYPE tau_f1;
    
  PDEFIELD_TYPE f2_inf;
  PDEFIELD_TYPE constf2 = 1.0;
  PDEFIELD_TYPE tau_f2;
    
  PDEFIELD_TYPE alpha_fCa;
  PDEFIELD_TYPE beta_fCa;
  PDEFIELD_TYPE gamma_fCa;
  PDEFIELD_TYPE fCa_inf;
    
  PDEFIELD_TYPE constfCa;
    
  PDEFIELD_TYPE tau_fCa     = 0.002;   // second (in i_CaL_fCa_gate)
    
  //// Ito
  PDEFIELD_TYPE g_to = 29.9038;   // S_per_F (in i_to)  
  PDEFIELD_TYPE i_to;
    
  PDEFIELD_TYPE q_inf;
  PDEFIELD_TYPE tau_q;
    
    
  PDEFIELD_TYPE r_inf;
  PDEFIELD_TYPE tau_r;
    
  //// IKs
  PDEFIELD_TYPE g_Ks = 2.041;   // S_per_F (in i_Ks)
  PDEFIELD_TYPE i_Ks; // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
  PDEFIELD_TYPE Xs_infinity;
  PDEFIELD_TYPE alpha_Xs;
  PDEFIELD_TYPE beta_Xs;
  PDEFIELD_TYPE tau_Xs;
    
  //// IKr
  PDEFIELD_TYPE L0 = 0.025;   // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE Q = 2.3;     // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE g_Kr = 29.8667;   // S_per_F (in i_Kr)
  PDEFIELD_TYPE i_Kr;
    
  PDEFIELD_TYPE V_half;
    
  PDEFIELD_TYPE Xr1_inf;
  PDEFIELD_TYPE alpha_Xr1;
  PDEFIELD_TYPE beta_Xr1;
  PDEFIELD_TYPE tau_Xr1;
    
  PDEFIELD_TYPE Xr2_infinity;
  PDEFIELD_TYPE alpha_Xr2;
  PDEFIELD_TYPE beta_Xr2;
  PDEFIELD_TYPE tau_Xr2;
    
  //// IK1
  PDEFIELD_TYPE alpha_K1;
  PDEFIELD_TYPE beta_K1;
  PDEFIELD_TYPE XK1_inf;
  PDEFIELD_TYPE g_K1 = 28.1492;   // S_per_F (in i_K1)
  PDEFIELD_TYPE i_K1;
    
  //// INaCa
  PDEFIELD_TYPE KmCa = 1.38;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE KmNai = 87.5;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE Ksat = 0.1;    // dimensionless (in i_NaCa)
  PDEFIELD_TYPE gamma = 0.35;   // dimensionless (in i_NaCa)
  PDEFIELD_TYPE alpha = 2.16659;
  PDEFIELD_TYPE kNaCa = 3917.0463; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaCaRedMed) * 6514.47574;   // A_per_F (in i_NaCa)
  PDEFIELD_TYPE i_NaCa;

  //// INaK
  PDEFIELD_TYPE Km_K = 1.0;    // millimolar (in i_NaK)
  PDEFIELD_TYPE Km_Na = 40.0;   // millimolar (in i_NaK)
  PDEFIELD_TYPE PNaK = 2.74240;// A_per_F (in i_NaK)
  PDEFIELD_TYPE i_NaK;
    
  //// IpCa
  PDEFIELD_TYPE KPCa = 0.0005;   // millimolar (in i_PCa)
  PDEFIELD_TYPE g_PCa = 0.4125;   // A_per_F (in i_PCa)
  PDEFIELD_TYPE i_PCa;
    
  //// Background currents
  PDEFIELD_TYPE g_b_Na = 1.14;         // S_per_F (in i_b_Na)
  PDEFIELD_TYPE i_b_Na;
    
  PDEFIELD_TYPE g_b_Ca = 0.8727264;    // S_per_F (in i_b_Ca)
  PDEFIELD_TYPE i_b_Ca;

  PDEFIELD_TYPE i_up;
  PDEFIELD_TYPE i_leak;
    
  //// Sarcoplasmic reticulum
  PDEFIELD_TYPE VmaxUp = 0.82205;
  PDEFIELD_TYPE Kup	= 4.40435e-4;
    
  PDEFIELD_TYPE V_leak = 4.48209e-4;
    
  // RyR
  PDEFIELD_TYPE g_irel_max = 55.808061;
  PDEFIELD_TYPE RyRa1 = 0.05169;
  PDEFIELD_TYPE RyRa2 = 0.050001;
  PDEFIELD_TYPE RyRahalf = 0.02632;
  PDEFIELD_TYPE RyRohalf = 0.00944;
  PDEFIELD_TYPE RyRchalf = 0.00167;
    
  PDEFIELD_TYPE RyRSRCass;
  PDEFIELD_TYPE i_rel;
    
  PDEFIELD_TYPE RyRainfss;
  PDEFIELD_TYPE RyRtauadapt = 1; //s
    
  PDEFIELD_TYPE RyRoinfss;
  PDEFIELD_TYPE RyRtauact;
    
  PDEFIELD_TYPE RyRcinfss;
  PDEFIELD_TYPE RyRtauinact;

  //// Ca2+ buffering
  PDEFIELD_TYPE Buf_C = 0.25;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Buf_SR = 10.0;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_C = 0.001;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_SR = 0.3;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Cai_bufc;
  PDEFIELD_TYPE Ca_SR_bufSR;
    
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
  
  for (int id = index; id < sizex*sizey; id += stride){
    if (celltype[id] < 1)
      for (int l = 0; l < layers; l++) //fill with current PDE values
        alt_PDEvars[l*sizex*sizey + id]= PDEvars[l*sizex*sizey + id];
    else{    
      for (int l = 0; l < layers; l++) //fill with current PDE values
        y[l] = PDEvars[l*sizex*sizey + id];

      //-----FIRST STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
      
        //// Nernst potential
      E_Na = R*T/F*log(Nao/y[17]);
      E_Ca = 0.5*R*T/F*log(Cao/y[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*y[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)y[13],3.0f)*y[11]*y[12]*(y[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((y[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((y[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((y[0]*1000 + 42)/-5)));
      dydt[13]   = (m_inf-y[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((y[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(y[0]*1000 + 79)/14)));
      dydt[11]   = (h_inf-y[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(y[0]*1000 + 79)/14)));
      dydt[12]   = (j_inf-y[12])/tau_j;
    
    
      //// INaL
      tauINaL     = 200; //ms
      GNaLmax     = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
      Vh_hLate    = 87.61;
      i_NaL       = GNaLmax* pow((float)y[18],3.0f)*y[19]*(y[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(y[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-y[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((y[0]*1000+35)/5))+0.1/(1+exp((y[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      dydt[18]   = (m_inf_L-y[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((y[0]*1000+Vh_hLate)/(7.488)));
      tau_h_L     = 1 * tauINaL;
      dydt[19]   = (h_inf_L-y[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*y[14]*(y[0] - E_K);
      i_fNa       = fNa*g_f*y[14]*(y[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((y[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((y[0]*1000 + 65)/7) + exp(-(y[0]*1000 + 65)/19));
      dydt[14]   = 1000*(Xf_infinity-y[14])/tau_Xf;
      
    
      //// ICaL
      //Prevent division by 0
      if(y[0]< precision && y[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T) *y[4]*y[5]*y[6]*y[7] / (2.0*F/(R*T)) * (y[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*y[0]*pow(F,2.0)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7];
    
      d_infinity  = 1.0/(1.0+exp(-(y[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-y[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((y[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-y[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      dydt[4]    = (d_infinity-y[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((y[0]*1000.0+26.0)/3.0));
      if (f1_inf-y[5] > 0.0)
          constf1 = 1.0+1433.0*(y[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((y[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-y[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf1/1000.0;
      dydt[5]    = (f1_inf-y[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((y[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(y[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-y[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf2/1000.0;
      dydt[6]    = (f2_inf-y[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(y[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((y[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((y[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((y[0] > -0.06) && (fCa_inf > y[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      dydt[7]    = constfCa*(fCa_inf-y[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(y[0]-E_K)*y[15]*y[16];
    
      q_inf       = 1.0/(1.0+exp((y[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(y[0]*1000.0+44.0))+0.065*exp(0.1*(y[0]*1000.0+45.93))))/1000.0;
      dydt[15]   = (q_inf-y[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(y[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(y[0]*1000.0+30.61))+0.369*exp(-0.12*(y[0]*1000.0+23.84))))/1000.0;
      dydt[16]   = (r_inf-y[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-y[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-y[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      dydt[10]   = (Xs_infinity-y[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-y[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-y[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+y[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      dydt[8]     = (Xr1_inf-y[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((y[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-y[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      dydt[9]    = (Xr2_infinity-y[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(y[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(y[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(y[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(y[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(y[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*y[0]*F/(R*T))*pow((float)y[17],3.0f)*Cao-exp((gamma-1.0)*y[0]*F/(R*T))*pow((float)Nao,3.0f)*y[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*y[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*y[17]/(y[17]+Km_Na)/(1.0+0.1245*exp(-0.1*y[0]*F/(R*T))+0.0353*exp(-y[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*y[2]/(y[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(y[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(y[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)y[2],2.0f));
    
      i_leak      = (y[1]-y[2])*V_leak;
    
      dydt[3]    = 0;
    
      // RyR
    
      RyRSRCass   = (1 - 1/(1 +  exp((y[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*y[21]*y[22]*(y[1]-y[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*y[2]-(RyRahalf))/0.0082));
      dydt[20]   = (RyRainfss- y[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*y[2]-(y[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= y[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      dydt[21]    = (RyRoinfss- y[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*y[2]-(y[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= y[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      dydt[22]    = (RyRcinfss- y[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(y[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(y[1]+Kbuf_SR), 2.0f));
    
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
    
      //// Membrane potential
      dydt[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);


    //-----SECOND STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+b21*dt*dydt[i];

      E_Na = R*T/F*log(Nao/ytemp[17]);
      E_Ca = 0.5*R*T/F*log(Cao/ytemp[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*ytemp[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)ytemp[13],3.0f)*ytemp[11]*ytemp[12]*(ytemp[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((ytemp[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((ytemp[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((ytemp[0]*1000 + 42)/-5)));
      ak2[13]   = (m_inf-ytemp[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((ytemp[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak2[11]   = (h_inf-ytemp[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak2[12]   = (j_inf-ytemp[12])/tau_j;
    
    
      //// INaL
      i_NaL       = GNaLmax* pow((float)ytemp[18],3.0f)*ytemp[19]*(ytemp[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(ytemp[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-ytemp[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((ytemp[0]*1000+35)/5))+0.1/(1+exp((ytemp[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      ak2[18]   = (m_inf_L-ytemp[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((ytemp[0]*1000+Vh_hLate)/(7.488)));
      ak2[19]   = (h_inf_L-ytemp[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*ytemp[14]*(ytemp[0] - E_K);
      i_fNa       = fNa*g_f*ytemp[14]*(ytemp[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((ytemp[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((ytemp[0]*1000 + 65)/7) + exp(-(ytemp[0]*1000 + 65)/19));
      ak2[14]   = 1000*(Xf_infinity-ytemp[14])/tau_Xf;
      
    
      //// ICaL
      //Prevent division by 0
      if(ytemp[0]< precision && ytemp[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T) *ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7] / (2.0*F/(R*T)) * (ytemp[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*ytemp[0]*pow(F,2.0)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7];
    
      d_infinity  = 1.0/(1.0+exp(-(ytemp[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-ytemp[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((ytemp[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-ytemp[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      ak2[4]    = (d_infinity-ytemp[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((ytemp[0]*1000.0+26.0)/3.0));
      if (f1_inf-ytemp[5] > 0.0)
          constf1 = 1.0+1433.0*(ytemp[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((ytemp[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-ytemp[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf1/1000.0;
      ak2[5]    = (f1_inf-ytemp[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((ytemp[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(ytemp[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-ytemp[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf2/1000.0;
      ak2[6]    = (f2_inf-ytemp[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(ytemp[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((ytemp[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((ytemp[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((ytemp[0] > -0.06) && (fCa_inf > ytemp[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      ak2[7]    = constfCa*(fCa_inf-ytemp[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(ytemp[0]-E_K)*ytemp[15]*ytemp[16];
    
      q_inf       = 1.0/(1.0+exp((ytemp[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(ytemp[0]*1000.0+44.0))+0.065*exp(0.1*(ytemp[0]*1000.0+45.93))))/1000.0;
      ak2[15]   = (q_inf-ytemp[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(ytemp[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(ytemp[0]*1000.0+30.61))+0.369*exp(-0.12*(ytemp[0]*1000.0+23.84))))/1000.0;
      ak2[16]   = (r_inf-ytemp[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-ytemp[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-ytemp[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      ak2[10]   = (Xs_infinity-ytemp[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-ytemp[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-ytemp[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+ytemp[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      ak2[8]     = (Xr1_inf-ytemp[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((ytemp[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-ytemp[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      ak2[9]    = (Xr2_infinity-ytemp[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(ytemp[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(ytemp[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(ytemp[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(ytemp[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(ytemp[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*ytemp[0]*F/(R*T))*pow((float)ytemp[17],3.0f)*Cao-exp((gamma-1.0)*ytemp[0]*F/(R*T))*pow((float)Nao,3.0f)*ytemp[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*ytemp[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*ytemp[17]/(ytemp[17]+Km_Na)/(1.0+0.1245*exp(-0.1*ytemp[0]*F/(R*T))+0.0353*exp(-ytemp[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*ytemp[2]/(ytemp[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(ytemp[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(ytemp[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)ytemp[2],2.0f));
    
      i_leak      = (ytemp[1]-ytemp[2])*V_leak;
    
      ak2[3]    = 0;
    
    
      RyRSRCass   = (1 - 1/(1 +  exp((ytemp[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*ytemp[21]*ytemp[22]*(ytemp[1]-ytemp[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*ytemp[2]-(RyRahalf))/0.0082));
      ak2[20]   = (RyRainfss- ytemp[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*ytemp[2]-(ytemp[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= ytemp[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      ak2[21]    = (RyRoinfss- ytemp[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*ytemp[2]-(ytemp[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= ytemp[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      ak2[22]    = (RyRcinfss- ytemp[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(ytemp[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(ytemp[1]+Kbuf_SR), 2.0f));
    
      //// Ionic concentrations
      //Nai
      ak2[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
      //Cai
      ak2[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
       //caSR
      ak2[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
    
      //// Stimulation
    //  i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
    //  i_stim_End 				= 1000.0;   // second (in stim_mode)
    //  i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
    //  i_stim_Start 			= 0.0;   // second (in stim_mode)
    //  i_stim_frequency        = 60.0;   // per_second (in stim_mode)
      //stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
    //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1
      //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
      //end
    
      /*
      if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
          i_stim = stim_flag*i_stim_Amplitude/Cm;
      else
          i_stim = 0.0;
      */
    
      //// Membrane potential
      ak2[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);


    //-----THIRD STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b31*dydt[i]+b32*ak2[i]);

      E_Na = R*T/F*log(Nao/ytemp[17]);
      E_Ca = 0.5*R*T/F*log(Cao/ytemp[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*ytemp[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)ytemp[13],3.0f)*ytemp[11]*ytemp[12]*(ytemp[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((ytemp[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((ytemp[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((ytemp[0]*1000 + 42)/-5)));
      ak3[13]   = (m_inf-ytemp[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((ytemp[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak3[11]   = (h_inf-ytemp[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak3[12]   = (j_inf-ytemp[12])/tau_j;
    
    
      //// INaL
      i_NaL       = GNaLmax* pow((float)ytemp[18],3.0f)*ytemp[19]*(ytemp[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(ytemp[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-ytemp[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((ytemp[0]*1000+35)/5))+0.1/(1+exp((ytemp[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      ak3[18]   = (m_inf_L-ytemp[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((ytemp[0]*1000+Vh_hLate)/(7.488)));
      ak3[19]   = (h_inf_L-ytemp[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*ytemp[14]*(ytemp[0] - E_K);
      i_fNa       = fNa*g_f*ytemp[14]*(ytemp[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((ytemp[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((ytemp[0]*1000 + 65)/7) + exp(-(ytemp[0]*1000 + 65)/19));
      ak3[14]   = 1000*(Xf_infinity-ytemp[14])/tau_Xf;
      
    
      //// ICaL  
      //Prevent division by 0
      if(ytemp[0]< precision && ytemp[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T) *ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7] / (2.0*F/(R*T)) * (ytemp[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*ytemp[0]*pow(F,2.0)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7];
    
      d_infinity  = 1.0/(1.0+exp(-(ytemp[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-ytemp[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((ytemp[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-ytemp[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      ak3[4]    = (d_infinity-ytemp[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((ytemp[0]*1000.0+26.0)/3.0));
      if (f1_inf-ytemp[5] > 0.0)
          constf1 = 1.0+1433.0*(ytemp[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((ytemp[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-ytemp[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf1/1000.0;
      ak3[5]    = (f1_inf-ytemp[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((ytemp[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(ytemp[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-ytemp[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf2/1000.0;
      ak3[6]    = (f2_inf-ytemp[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(ytemp[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((ytemp[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((ytemp[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((ytemp[0] > -0.06) && (fCa_inf > ytemp[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;

      ak3[7]    = constfCa*(fCa_inf-ytemp[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(ytemp[0]-E_K)*ytemp[15]*ytemp[16];
    
      q_inf       = 1.0/(1.0+exp((ytemp[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(ytemp[0]*1000.0+44.0))+0.065*exp(0.1*(ytemp[0]*1000.0+45.93))))/1000.0;
      ak3[15]   = (q_inf-ytemp[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(ytemp[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(ytemp[0]*1000.0+30.61))+0.369*exp(-0.12*(ytemp[0]*1000.0+23.84))))/1000.0;
      ak3[16]   = (r_inf-ytemp[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-ytemp[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-ytemp[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      ak3[10]   = (Xs_infinity-ytemp[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-ytemp[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-ytemp[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+ytemp[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      ak3[8]     = (Xr1_inf-ytemp[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((ytemp[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-ytemp[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      ak3[9]    = (Xr2_infinity-ytemp[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(ytemp[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(ytemp[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(ytemp[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(ytemp[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(ytemp[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*ytemp[0]*F/(R*T))*pow((float)ytemp[17],3.0f)*Cao-exp((gamma-1.0)*ytemp[0]*F/(R*T))*pow((float)Nao,3.0f)*ytemp[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*ytemp[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*ytemp[17]/(ytemp[17]+Km_Na)/(1.0+0.1245*exp(-0.1*ytemp[0]*F/(R*T))+0.0353*exp(-ytemp[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*ytemp[2]/(ytemp[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(ytemp[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(ytemp[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)ytemp[2],2.0f));
    
      i_leak      = (ytemp[1]-ytemp[2])*V_leak;
    
      ak3[3]    = 0;
  
    
      RyRSRCass   = (1 - 1/(1 +  exp((ytemp[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*ytemp[21]*ytemp[22]*(ytemp[1]-ytemp[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*ytemp[2]-(RyRahalf))/0.0082));
      ak3[20]   = (RyRainfss- ytemp[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*ytemp[2]-(ytemp[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= ytemp[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      ak3[21]    = (RyRoinfss- ytemp[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*ytemp[2]-(ytemp[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= ytemp[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      ak3[22]    = (RyRcinfss- ytemp[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(ytemp[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(ytemp[1]+Kbuf_SR), 2.0f));
    
      //// Ionic concentrations
      //Nai
      ak3[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
      //Cai
      ak3[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
       //caSR
      ak3[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
    
      //// Stimulation
    //  i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
    //  i_stim_End 				= 1000.0;   // second (in stim_mode)
    //  i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
    //  i_stim_Start 			= 0.0;   // second (in stim_mode)
    //  i_stim_frequency        = 60.0;   // per_second (in stim_mode)
      //stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
    //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1
      //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
      //end
    
      /*
      if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
          i_stim = stim_flag*i_stim_Amplitude/Cm;
      else
          i_stim = 0.0;
      */
    
      //// Membrane potential
      ak3[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);


    //-----FOURTH STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b41*dydt[i]+b42*ak2[i]+b43*ak3[i]);
    
      E_Na = R*T/F*log(Nao/ytemp[17]);
      E_Ca = 0.5*R*T/F*log(Cao/ytemp[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*ytemp[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)ytemp[13],3.0f)*ytemp[11]*ytemp[12]*(ytemp[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((ytemp[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((ytemp[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((ytemp[0]*1000 + 42)/-5)));
      ak4[13]   = (m_inf-ytemp[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((ytemp[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak4[11]   = (h_inf-ytemp[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak4[12]   = (j_inf-ytemp[12])/tau_j;
    
    
      //// INaL
      i_NaL       = GNaLmax* pow((float)ytemp[18],3.0f)*ytemp[19]*(ytemp[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(ytemp[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-ytemp[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((ytemp[0]*1000+35)/5))+0.1/(1+exp((ytemp[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      ak4[18]   = (m_inf_L-ytemp[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((ytemp[0]*1000+Vh_hLate)/(7.488)));
      ak4[19]   = (h_inf_L-ytemp[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*ytemp[14]*(ytemp[0] - E_K);
      i_fNa       = fNa*g_f*ytemp[14]*(ytemp[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((ytemp[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((ytemp[0]*1000 + 65)/7) + exp(-(ytemp[0]*1000 + 65)/19));
      ak4[14]   = 1000*(Xf_infinity-ytemp[14])/tau_Xf;
      
    
      //// ICaL
      //Prevent division by 0
      if(ytemp[0]< precision && ytemp[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T) *ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7] / (2.0*F/(R*T)) * (ytemp[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*ytemp[0]*pow(F,2.0)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7];
    
      d_infinity  = 1.0/(1.0+exp(-(ytemp[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-ytemp[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((ytemp[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-ytemp[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      ak4[4]    = (d_infinity-ytemp[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((ytemp[0]*1000.0+26.0)/3.0));
      if (f1_inf-ytemp[5] > 0.0)
          constf1 = 1.0+1433.0*(ytemp[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((ytemp[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-ytemp[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf1/1000.0;
      ak4[5]    = (f1_inf-ytemp[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((ytemp[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(ytemp[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-ytemp[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf2/1000.0;
      ak4[6]    = (f2_inf-ytemp[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(ytemp[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((ytemp[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((ytemp[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((ytemp[0] > -0.06) && (fCa_inf > ytemp[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      ak4[7]    = constfCa*(fCa_inf-ytemp[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(ytemp[0]-E_K)*ytemp[15]*ytemp[16];
    
      q_inf       = 1.0/(1.0+exp((ytemp[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(ytemp[0]*1000.0+44.0))+0.065*exp(0.1*(ytemp[0]*1000.0+45.93))))/1000.0;
      ak4[15]   = (q_inf-ytemp[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(ytemp[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(ytemp[0]*1000.0+30.61))+0.369*exp(-0.12*(ytemp[0]*1000.0+23.84))))/1000.0;
      ak4[16]   = (r_inf-ytemp[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-ytemp[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-ytemp[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      ak4[10]   = (Xs_infinity-ytemp[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-ytemp[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-ytemp[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+ytemp[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      ak4[8]     = (Xr1_inf-ytemp[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((ytemp[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-ytemp[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      ak4[9]    = (Xr2_infinity-ytemp[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(ytemp[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(ytemp[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(ytemp[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(ytemp[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(ytemp[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*ytemp[0]*F/(R*T))*pow((float)ytemp[17],3.0f)*Cao-exp((gamma-1.0)*ytemp[0]*F/(R*T))*pow((float)Nao,3.0f)*ytemp[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*ytemp[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*ytemp[17]/(ytemp[17]+Km_Na)/(1.0+0.1245*exp(-0.1*ytemp[0]*F/(R*T))+0.0353*exp(-ytemp[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*ytemp[2]/(ytemp[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(ytemp[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(ytemp[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)ytemp[2],2.0f));
    
      i_leak      = (ytemp[1]-ytemp[2])*V_leak;
    
      ak4[3]    = 0;
    
    
      RyRSRCass   = (1 - 1/(1 +  exp((ytemp[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*ytemp[21]*ytemp[22]*(ytemp[1]-ytemp[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*ytemp[2]-(RyRahalf))/0.0082));
      ak4[20]   = (RyRainfss- ytemp[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*ytemp[2]-(ytemp[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= ytemp[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      ak4[21]    = (RyRoinfss- ytemp[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*ytemp[2]-(ytemp[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= ytemp[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      ak4[22]    = (RyRcinfss- ytemp[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(ytemp[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(ytemp[1]+Kbuf_SR), 2.0f));
    
      //// Ionic concentrations
      //Nai
      ak4[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
      //Cai
      ak4[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
       //caSR
      ak4[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
    
      //// Stimulation
    //  i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
    //  i_stim_End 				= 1000.0;   // second (in stim_mode)
    //  i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
    //  i_stim_Start 			= 0.0;   // second (in stim_mode)
    //  i_stim_frequency        = 60.0;   // per_second (in stim_mode)
      //stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
    //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1
      //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
      //end
    
      /*
      if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
          i_stim = stim_flag*i_stim_Amplitude/Cm;
      else
          i_stim = 0.0;
      */
    
      //// Membrane potential
      ak4[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);
    //-----FIFTH STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b51*dydt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);

      E_Na = R*T/F*log(Nao/ytemp[17]);
      E_Ca = 0.5*R*T/F*log(Cao/ytemp[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*ytemp[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)ytemp[13],3.0f)*ytemp[11]*ytemp[12]*(ytemp[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((ytemp[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((ytemp[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((ytemp[0]*1000 + 42)/-5)));
      ak5[13]   = (m_inf-ytemp[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((ytemp[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak5[11]   = (h_inf-ytemp[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak5[12]   = (j_inf-ytemp[12])/tau_j;
    
    
      //// INaL
      i_NaL       = GNaLmax* pow((float)ytemp[18],3.0f)*ytemp[19]*(ytemp[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(ytemp[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-ytemp[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((ytemp[0]*1000+35)/5))+0.1/(1+exp((ytemp[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      ak5[18]   = (m_inf_L-ytemp[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((ytemp[0]*1000+Vh_hLate)/(7.488)));
      ak5[19]   = (h_inf_L-ytemp[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*ytemp[14]*(ytemp[0] - E_K);
      i_fNa       = fNa*g_f*ytemp[14]*(ytemp[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((ytemp[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((ytemp[0]*1000 + 65)/7) + exp(-(ytemp[0]*1000 + 65)/19));
      ak5[14]   = 1000*(Xf_infinity-ytemp[14])/tau_Xf;
      
    
      //// ICaL 
      //Prevent division by 0
      if(ytemp[0]< precision && ytemp[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T) *ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7] / (2.0*F/(R*T)) * (ytemp[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*ytemp[0]*pow(F,2.0)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7];
    
      d_infinity  = 1.0/(1.0+exp(-(ytemp[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-ytemp[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((ytemp[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-ytemp[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      ak5[4]    = (d_infinity-ytemp[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((ytemp[0]*1000.0+26.0)/3.0));
      if (f1_inf-ytemp[5] > 0.0)
          constf1 = 1.0+1433.0*(ytemp[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((ytemp[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-ytemp[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf1/1000.0;
      ak5[5]    = (f1_inf-ytemp[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((ytemp[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(ytemp[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-ytemp[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf2/1000.0;
      ak5[6]    = (f2_inf-ytemp[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(ytemp[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((ytemp[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((ytemp[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((ytemp[0] > -0.06) && (fCa_inf > ytemp[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      ak5[7]    = constfCa*(fCa_inf-ytemp[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(ytemp[0]-E_K)*ytemp[15]*ytemp[16];
    
      q_inf       = 1.0/(1.0+exp((ytemp[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(ytemp[0]*1000.0+44.0))+0.065*exp(0.1*(ytemp[0]*1000.0+45.93))))/1000.0;
      ak5[15]   = (q_inf-ytemp[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(ytemp[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(ytemp[0]*1000.0+30.61))+0.369*exp(-0.12*(ytemp[0]*1000.0+23.84))))/1000.0;
      ak5[16]   = (r_inf-ytemp[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-ytemp[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-ytemp[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      ak5[10]   = (Xs_infinity-ytemp[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-ytemp[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-ytemp[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+ytemp[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      ak5[8]     = (Xr1_inf-ytemp[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((ytemp[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-ytemp[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      ak5[9]    = (Xr2_infinity-ytemp[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(ytemp[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(ytemp[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(ytemp[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(ytemp[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(ytemp[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*ytemp[0]*F/(R*T))*pow((float)ytemp[17],3.0f)*Cao-exp((gamma-1.0)*ytemp[0]*F/(R*T))*pow((float)Nao,3.0f)*ytemp[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*ytemp[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*ytemp[17]/(ytemp[17]+Km_Na)/(1.0+0.1245*exp(-0.1*ytemp[0]*F/(R*T))+0.0353*exp(-ytemp[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*ytemp[2]/(ytemp[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(ytemp[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(ytemp[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)ytemp[2],2.0f));
    
      i_leak      = (ytemp[1]-ytemp[2])*V_leak;
    
      ak5[3]    = 0;
  
    
      RyRSRCass   = (1 - 1/(1 +  exp((ytemp[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*ytemp[21]*ytemp[22]*(ytemp[1]-ytemp[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*ytemp[2]-(RyRahalf))/0.0082));
      ak5[20]   = (RyRainfss- ytemp[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*ytemp[2]-(ytemp[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= ytemp[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      ak5[21]    = (RyRoinfss- ytemp[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*ytemp[2]-(ytemp[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= ytemp[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      ak5[22]    = (RyRcinfss- ytemp[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(ytemp[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(ytemp[1]+Kbuf_SR), 2.0f));
    
      //// Ionic concentrations
      //Nai
      ak5[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
      //Cai
      ak5[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
       //caSR
      ak5[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
    
      //// Stimulation
    //  i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
    //  i_stim_End 				= 1000.0;   // second (in stim_mode)
    //  i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
    //  i_stim_Start 			= 0.0;   // second (in stim_mode)
    //  i_stim_frequency        = 60.0;   // per_second (in stim_mode)
      //stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
    //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1
      //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
      //end
    
      /*
      if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
          i_stim = stim_flag*i_stim_Amplitude/Cm;
      else
          i_stim = 0.0;
      */
    
      //// Membrane potential
      ak5[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);      

    //-----SIXTH STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------      
    for (i=0;i<layers;i++)
      ytemp[i]=y[i]+dt*(b61*dydt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);


      E_Na = R*T/F*log(Nao/ytemp[17]);
      E_Ca = 0.5*R*T/F*log(Cao/ytemp[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*ytemp[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)ytemp[13],3.0f)*ytemp[11]*ytemp[12]*(ytemp[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((ytemp[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((ytemp[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((ytemp[0]*1000 + 42)/-5)));
      ak6[13]   = (m_inf-ytemp[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((ytemp[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak6[11]   = (h_inf-ytemp[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((ytemp[0]*1000 + 41)/5.5) + exp(-(ytemp[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(ytemp[0]*1000 + 79)/14)));
      ak6[12]   = (j_inf-ytemp[12])/tau_j;
    
    
      //// INaL
      i_NaL       = GNaLmax* pow((float)ytemp[18],3.0f)*ytemp[19]*(ytemp[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(ytemp[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-ytemp[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((ytemp[0]*1000+35)/5))+0.1/(1+exp((ytemp[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      ak6[18]   = (m_inf_L-ytemp[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((ytemp[0]*1000+Vh_hLate)/(7.488)));
      ak6[19]   = (h_inf_L-ytemp[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*ytemp[14]*(ytemp[0] - E_K);
      i_fNa       = fNa*g_f*ytemp[14]*(ytemp[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((ytemp[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((ytemp[0]*1000 + 65)/7) + exp(-(ytemp[0]*1000 + 65)/19));
      ak6[14]   = 1000*(Xf_infinity-ytemp[14])/tau_Xf;
      
    
      //// ICaL
      //Prevent division by 0
      if(ytemp[0]< precision && ytemp[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T) *ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7] / (2.0*F/(R*T)) * (ytemp[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*ytemp[0]*pow((float)F,2.0f)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*ytemp[0]*pow(F,2.0)/(R*T)*(ytemp[2]*exp(2.0*ytemp[0]*F/(R*T))-0.341*Cao)/(exp(2.0*ytemp[0]*F/(R*T))-1.0)*ytemp[4]*ytemp[5]*ytemp[6]*ytemp[7];
    
      d_infinity  = 1.0/(1.0+exp(-(ytemp[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-ytemp[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((ytemp[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-ytemp[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      ak6[4]    = (d_infinity-ytemp[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((ytemp[0]*1000.0+26.0)/3.0));
      if (f1_inf-ytemp[5] > 0.0)
          constf1 = 1.0+1433.0*(ytemp[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((ytemp[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-ytemp[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf1/1000.0;
      ak6[5]    = (f1_inf-ytemp[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((ytemp[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(ytemp[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-ytemp[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+ytemp[0]*1000.0)/10.0)))*constf2/1000.0;
      ak6[6]    = (f2_inf-ytemp[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(ytemp[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((ytemp[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((ytemp[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((ytemp[0] > -0.06) && (fCa_inf > ytemp[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      ak6[7]    = constfCa*(fCa_inf-ytemp[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(ytemp[0]-E_K)*ytemp[15]*ytemp[16];
    
      q_inf       = 1.0/(1.0+exp((ytemp[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(ytemp[0]*1000.0+44.0))+0.065*exp(0.1*(ytemp[0]*1000.0+45.93))))/1000.0;
      ak6[15]   = (q_inf-ytemp[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(ytemp[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(ytemp[0]*1000.0+30.61))+0.369*exp(-0.12*(ytemp[0]*1000.0+23.84))))/1000.0;
      ak6[16]   = (r_inf-ytemp[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(ytemp[0]-E_Ks)*pow((float)ytemp[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/ytemp[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-ytemp[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-ytemp[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      ak6[10]   = (Xs_infinity-ytemp[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(ytemp[0]-E_K)*ytemp[8]*ytemp[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-ytemp[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-ytemp[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+ytemp[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      ak6[8]     = (Xr1_inf-ytemp[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((ytemp[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-ytemp[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+ytemp[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      ak6[9]    = (Xr2_infinity-ytemp[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(ytemp[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(ytemp[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(ytemp[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(ytemp[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(ytemp[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*ytemp[0]*F/(R*T))*pow((float)ytemp[17],3.0f)*Cao-exp((gamma-1.0)*ytemp[0]*F/(R*T))*pow((float)Nao,3.0f)*ytemp[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*ytemp[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*ytemp[17]/(ytemp[17]+Km_Na)/(1.0+0.1245*exp(-0.1*ytemp[0]*F/(R*T))+0.0353*exp(-ytemp[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*ytemp[2]/(ytemp[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(ytemp[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(ytemp[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)ytemp[2],2.0f));
    
      i_leak      = (ytemp[1]-ytemp[2])*V_leak;
    
      ak6[3]    = 0;
    
    
      RyRSRCass   = (1 - 1/(1 +  exp((ytemp[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*ytemp[21]*ytemp[22]*(ytemp[1]-ytemp[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*ytemp[2]-(RyRahalf))/0.0082));
      ak6[20]   = (RyRainfss- ytemp[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*ytemp[2]-(ytemp[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= ytemp[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      ak6[21]    = (RyRoinfss- ytemp[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*ytemp[2]-(ytemp[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= ytemp[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      ak6[22]    = (RyRcinfss- ytemp[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(ytemp[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(ytemp[1]+Kbuf_SR), 2.0f));
    
      //// Ionic concentrations
      //Nai
      ak6[17]   = -Cm*(i_Na+i_NaL+i_b_Na+3.0*i_NaK+3.0*i_NaCa+i_fNa)/(F*Vc*1.0e-18);
      //Cai
      ak6[2]    = Cai_bufc*(i_leak-i_up+i_rel-(i_CaL+i_b_Ca+i_PCa-2.0*i_NaCa)*Cm/(2.0*Vc*F*1.0e-18));
       //caSR
      ak6[1]    = Ca_SR_bufSR*Vc/V_SR*(i_up-(i_rel+i_leak));
    
      //// Stimulation
      //  i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
      //  i_stim_End 				= 1000.0;   // second (in stim_mode)
      //  i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
      //  i_stim_Start 			= 0.0;   // second (in stim_mode)
      //  i_stim_frequency        = 60.0;   // per_second (in stim_mode)
      //stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
      //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1_mode)
      //  i_stim_Period 			= 60.0/i_stim_frequency;
    
      //if stim_flag~=0 && stim_flag~=1
      //error('Paci2020: wrong pacing! stimFlag can be only 0 (spontaneous) or 1 (paced)');
      //end
    
      /*
      if ((time >= i_stim_Start) && (time <= i_stim_End) && (time-i_stim_Start-floor((time-i_stim_Start)/i_stim_Period)*i_stim_Period <= i_stim_PulseDuration))
          i_stim = stim_flag*i_stim_Amplitude/Cm;
      else
          i_stim = 0.0;
      */
    
      //// Membrane potential
      ak6[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);
  
      //-----WRITE NEW VALUES ------------------------------------------------------------------------------------------------------------------------------------------------------------------      
      for (i=0;i<layers;i++){ //Accumulate increments with proper weights.
        yout[i]=y[i]+dt*(c1*dydt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
        yerr[i]=dt*(dc1*dydt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
      }

      #pragma unroll
      for (i=0;i<layers;i++) //Accumulate increments with proper weights.
        alt_PDEvars[i*sizex*sizey + id]=y[i]+(c1*dydt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i])*dt;
    }  
  
  }
}
#endif

__device__ void derivsPaci(PDEFIELD_TYPE current_time, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength){

  //Declare variables needed for Paci2020 model and assign the constants

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
    
      //// Nernst potential
  PDEFIELD_TYPE E_Na;
  PDEFIELD_TYPE E_Ca;
  PDEFIELD_TYPE E_K;
  PDEFIELD_TYPE PkNa = 0.03;   // dimensionless (in electric_potentials)
  PDEFIELD_TYPE E_Ks;
    
  //// INa adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_Na = 3671.2302; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*6447.1896;
  PDEFIELD_TYPE i_Na;
    
  PDEFIELD_TYPE m_inf;
  PDEFIELD_TYPE tau_m;
    
  PDEFIELD_TYPE h_inf;
  PDEFIELD_TYPE tau_h;
    
  PDEFIELD_TYPE j_inf;
  PDEFIELD_TYPE tau_j;
    
    
  //// INaL
  PDEFIELD_TYPE myCoefTauM  = 1;
  PDEFIELD_TYPE tauINaL = 200; //ms
  PDEFIELD_TYPE GNaLmax = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
  PDEFIELD_TYPE Vh_hLate = 87.61;
  PDEFIELD_TYPE i_NaL;
    
  PDEFIELD_TYPE m_inf_L;
  PDEFIELD_TYPE alpha_m_L;
  PDEFIELD_TYPE beta_m_L;
  PDEFIELD_TYPE tau_m_L;
    
  PDEFIELD_TYPE h_inf_L;
  PDEFIELD_TYPE tau_h_L = 1 * tauINaL;
    
  //// If adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_f = 1; //((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*22.2763088;
  PDEFIELD_TYPE fNa = 0.37;
  PDEFIELD_TYPE fK = 1 - fNa;
  PDEFIELD_TYPE i_fK;
  PDEFIELD_TYPE i_fNa;
  PDEFIELD_TYPE i_f;
    
  PDEFIELD_TYPE Xf_infinity;
  PDEFIELD_TYPE tau_Xf; 
    
      //// ICaL
  PDEFIELD_TYPE g_CaL = 8.635702e-5;   // metre_cube_per_F_per_s (in i_CaL)
  PDEFIELD_TYPE i_CaL;  
  PDEFIELD_TYPE precision = 0.0001;     
    
  PDEFIELD_TYPE d_infinity;
  PDEFIELD_TYPE alpha_d;
  PDEFIELD_TYPE beta_d;
  PDEFIELD_TYPE gamma_d;
  PDEFIELD_TYPE tau_d;
    
  PDEFIELD_TYPE f1_inf;
  PDEFIELD_TYPE constf1;
    
  PDEFIELD_TYPE tau_f1;
    
  PDEFIELD_TYPE f2_inf;
  PDEFIELD_TYPE constf2 = 1.0;
  PDEFIELD_TYPE tau_f2;
    
  PDEFIELD_TYPE alpha_fCa;
  PDEFIELD_TYPE beta_fCa;
  PDEFIELD_TYPE gamma_fCa;
  PDEFIELD_TYPE fCa_inf;
    
  PDEFIELD_TYPE constfCa;
    
  PDEFIELD_TYPE tau_fCa     = 0.002;   // second (in i_CaL_fCa_gate)
    
  //// Ito
  PDEFIELD_TYPE g_to = 29.9038;   // S_per_F (in i_to)  
  PDEFIELD_TYPE i_to;
    
  PDEFIELD_TYPE q_inf;
  PDEFIELD_TYPE tau_q;
    
    
  PDEFIELD_TYPE r_inf;
  PDEFIELD_TYPE tau_r;
    
  //// IKs
  PDEFIELD_TYPE g_Ks = 2.041;   // S_per_F (in i_Ks)
  PDEFIELD_TYPE i_Ks; // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
  PDEFIELD_TYPE Xs_infinity;
  PDEFIELD_TYPE alpha_Xs;
  PDEFIELD_TYPE beta_Xs;
  PDEFIELD_TYPE tau_Xs;
    
  //// IKr
  PDEFIELD_TYPE L0 = 0.025;   // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE Q = 2.3;     // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE g_Kr = 29.8667;   // S_per_F (in i_Kr)
  PDEFIELD_TYPE i_Kr;
    
  PDEFIELD_TYPE V_half;
    
  PDEFIELD_TYPE Xr1_inf;
  PDEFIELD_TYPE alpha_Xr1;
  PDEFIELD_TYPE beta_Xr1;
  PDEFIELD_TYPE tau_Xr1;
    
  PDEFIELD_TYPE Xr2_infinity;
  PDEFIELD_TYPE alpha_Xr2;
  PDEFIELD_TYPE beta_Xr2;
  PDEFIELD_TYPE tau_Xr2;
    
  //// IK1
  PDEFIELD_TYPE alpha_K1;
  PDEFIELD_TYPE beta_K1;
  PDEFIELD_TYPE XK1_inf;
  PDEFIELD_TYPE g_K1 = 28.1492;   // S_per_F (in i_K1)
  PDEFIELD_TYPE i_K1;
    
  //// INaCa
  PDEFIELD_TYPE KmCa = 1.38;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE KmNai = 87.5;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE Ksat = 0.1;    // dimensionless (in i_NaCa)
  PDEFIELD_TYPE gamma = 0.35;   // dimensionless (in i_NaCa)
  PDEFIELD_TYPE alpha = 2.16659;
  PDEFIELD_TYPE kNaCa = 3917.0463; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaCaRedMed) * 6514.47574;   // A_per_F (in i_NaCa)
  PDEFIELD_TYPE i_NaCa;

  //// INaK
  PDEFIELD_TYPE Km_K = 1.0;    // millimolar (in i_NaK)
  PDEFIELD_TYPE Km_Na = 40.0;   // millimolar (in i_NaK)
  PDEFIELD_TYPE PNaK = 2.74240;// A_per_F (in i_NaK)
  PDEFIELD_TYPE i_NaK;
    
  //// IpCa
  PDEFIELD_TYPE KPCa = 0.0005;   // millimolar (in i_PCa)
  PDEFIELD_TYPE g_PCa = 0.4125;   // A_per_F (in i_PCa)
  PDEFIELD_TYPE i_PCa;
    
  //// Background currents
  PDEFIELD_TYPE g_b_Na = 1.14;         // S_per_F (in i_b_Na)
  PDEFIELD_TYPE i_b_Na;
    
  PDEFIELD_TYPE g_b_Ca = 0.8727264;    // S_per_F (in i_b_Ca)
  PDEFIELD_TYPE i_b_Ca;

  PDEFIELD_TYPE i_up;
  PDEFIELD_TYPE i_leak;
    
  //// Sarcoplasmic reticulum
  PDEFIELD_TYPE VmaxUp = 0.82205;
  PDEFIELD_TYPE Kup	= 4.40435e-4;
    
  PDEFIELD_TYPE V_leak = 4.48209e-4;
    
  // RyR
  PDEFIELD_TYPE g_irel_max = 55.808061;
  PDEFIELD_TYPE RyRa1 = 0.05169;
  PDEFIELD_TYPE RyRa2 = 0.050001;
  PDEFIELD_TYPE RyRahalf = 0.02632;
  PDEFIELD_TYPE RyRohalf = 0.00944;
  PDEFIELD_TYPE RyRchalf = 0.00167;
    
  PDEFIELD_TYPE RyRSRCass;
  PDEFIELD_TYPE i_rel;
    
  PDEFIELD_TYPE RyRainfss;
  PDEFIELD_TYPE RyRtauadapt = 1; //s
    
  PDEFIELD_TYPE RyRoinfss;
  PDEFIELD_TYPE RyRtauact;
    
  PDEFIELD_TYPE RyRcinfss;
  PDEFIELD_TYPE RyRtauinact;

  //// Ca2+ buffering
  PDEFIELD_TYPE Buf_C = 0.25;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Buf_SR = 10.0;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_C = 0.001;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_SR = 0.3;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Cai_bufc;
  PDEFIELD_TYPE Ca_SR_bufSR;
    
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
  

  //-----FIRST STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
   
  //// Nernst potential
  E_Na = R*T/F*log(Nao/y[17]);
  E_Ca = 0.5*R*T/F*log(Cao/y[2]);
  E_K  = R*T/F*log(Ko/Ki);
  E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*y[17]));
    
  //// INa adapted from DOI:10.3389/fphys.2018.00080
  i_Na        =  g_Na*pow((float)y[13],3.0f)*y[11]*y[12]*(y[0] - E_Na);
    
  m_inf       = 1 / (1 + exp((y[0]*1000 + 39)/-11.2));
  tau_m       = (0.00001 + 0.00013*exp(-pow((float)((y[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((y[0]*1000 + 42)/-5)));
  dydt[13]   = (m_inf-y[13])/tau_m;
    
  h_inf       = 1 / (1 + exp((y[0]*1000 + 66.5)/6.8));
  tau_h       = (0.00007 + 0.034 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(y[0]*1000 + 79)/14)));
  dydt[11]   = (h_inf-y[11])/tau_h;
    
  j_inf       = h_inf;
  tau_j       = 10*(0.0007 + 0.15 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(y[0]*1000 + 79)/14)));
  dydt[12]   = (j_inf-y[12])/tau_j;
    
    
  //// INaL
  tauINaL     = 200; //ms
  GNaLmax     = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
  Vh_hLate    = 87.61;
  i_NaL       = GNaLmax* pow((float)y[18],3.0f)*y[19]*(y[0]-E_Na);
    
  m_inf_L     = 1/(1+exp(-(y[0]*1000+42.85)/(5.264)));
  alpha_m_L   = 1/(1+exp((-60-y[0]*1000)/5));
  beta_m_L    = 0.1/(1+exp((y[0]*1000+35)/5))+0.1/(1+exp((y[0]*1000-50)/200));
  tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
  dydt[18]   = (m_inf_L-y[18])/tau_m_L*1000;
    
  h_inf_L     = 1/(1+exp((y[0]*1000+Vh_hLate)/(7.488)));
  tau_h_L     = 1 * tauINaL;
  dydt[19]   = (h_inf_L-y[19])/tau_h_L*1000;
    
  //// If adapted from DOI:10.3389/fphys.2018.00080
  i_fK        = fK*g_f*y[14]*(y[0] - E_K);
  i_fNa       = fNa*g_f*y[14]*(y[0] - E_Na);
  i_f         = i_fK + i_fNa;
    
  Xf_infinity = 1.0/(1.0 + exp((y[0]*1000 + 69)/8));
  tau_Xf      = 5600 / (1 + exp((y[0]*1000 + 65)/7) + exp(-(y[0]*1000 + 65)/19));
  dydt[14]   = 1000*(Xf_infinity-y[14])/tau_Xf;
      
    
  //// ICaL
  //Prevent division by 0
  if(y[0]< precision && y[0] > -precision) //hopital
    i_CaL =  g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T) *y[4]*y[5]*y[6]*y[7] / (2.0*F/(R*T)) * (y[2] - 0.341*Cao);
  else
    i_CaL = g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*y[0]*pow(F,2.0)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7];
    
  d_infinity  = 1.0/(1.0+exp(-(y[0]*1000.0+9.1)/7.0));
  alpha_d     = 0.25+1.4/(1.0+exp((-y[0]*1000.0-35.0)/13.0));
  beta_d      = 1.4/(1.0+exp((y[0]*1000.0+5.0)/5.0));
  gamma_d     = 1.0/(1.0+exp((-y[0]*1000.0+50.0)/20.0));
  tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
  dydt[4]    = (d_infinity-y[4])/tau_d;
    
  f1_inf      = 1.0/(1.0+exp((y[0]*1000.0+26.0)/3.0));
  if (f1_inf-y[5] > 0.0)
      constf1 = 1.0+1433.0*(y[2]-50.0*1.0e-6);
  else
      constf1 = 1.0;
  
  tau_f1      = (20.0+1102.5*exp(-pow((float)((y[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-y[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf1/1000.0;
  dydt[5]    = (f1_inf-y[5])/tau_f1;
   
  f2_inf      = 0.33+0.67/(1.0+exp((y[0]*1000.0+32.0)/4.0));
  tau_f2      = (600.0*exp(-pow((float)(y[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-y[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf2/1000.0;
  dydt[6]    = (f2_inf-y[6])/tau_f2;
   
  alpha_fCa   = 1.0/(1.0+pow((float)(y[2]/0.0006),8.0f));
  beta_fCa    = 0.1/(1.0+exp((y[2]-0.0009)/0.0001));
  gamma_fCa   = 0.3/(1.0+exp((y[2]-0.00075)/0.0008));
  fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
  if ((y[0] > -0.06) && (fCa_inf > y[7]))
    constfCa = 0.0;
  else
    constfCa = 1.0;
    
  dydt[7]    = constfCa*(fCa_inf-y[7])/tau_fCa;
    
  //// Ito
  i_to        = g_to*(y[0]-E_K)*y[15]*y[16];
    
  q_inf       = 1.0/(1.0+exp((y[0]*1000.0+53.0)/13.0));
  tau_q       = (6.06+39.102/(0.57*exp(-0.08*(y[0]*1000.0+44.0))+0.065*exp(0.1*(y[0]*1000.0+45.93))))/1000.0;
  dydt[15]   = (q_inf-y[15])/tau_q;
    
    
  r_inf       = 1.0/(1.0+exp(-(y[0]*1000.0-22.3)/18.75));
  tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(y[0]*1000.0+30.61))+0.369*exp(-0.12*(y[0]*1000.0+23.84))))/1000.0;
  dydt[16]   = (r_inf-y[16])/tau_r;
    
  //// IKs
  i_Ks        = g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
  Xs_infinity = 1.0/(1.0+exp((-y[0]*1000.0-20.0)/16.0));
  alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-y[0]*1000.0)/6.0));
  beta_Xs     = 1.0/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
  tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
  dydt[10]   = (Xs_infinity-y[10])/tau_Xs;
    
  //// IKr
  i_Kr         = g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4);
    
  V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
  Xr1_inf      = 1.0/(1.0+exp((V_half-y[0]*1000.0)/4.9));
  alpha_Xr1    = 450.0/(1.0+exp((-45.0-y[0]*1000.0)/10.0));
  beta_Xr1     = 6.0/(1.0+exp((30.0+y[0]*1000.0)/11.5));
  tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
  dydt[8]     = (Xr1_inf-y[8])/tau_Xr1;
    
  Xr2_infinity = 1.0/(1.0+exp((y[0]*1000.0+88.0)/50.0));
  alpha_Xr2    = 3.0/(1.0+exp((-60.0-y[0]*1000.0)/20.0));
  beta_Xr2     = 1.12/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
  tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
  dydt[9]    = (Xr2_infinity-y[9])/tau_Xr2;
    
  //// IK1
  alpha_K1    = 3.91/(1.0+exp(0.5942*(y[0]*1000.0-E_K*1000.0-200.0)));
  beta_K1     = (-1.509*exp(0.0002*(y[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(y[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(y[0]*1000.0-E_K*1000.0)));
  XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
  i_K1        = g_K1*XK1_inf*(y[0]-E_K)*sqrt(Ko/5.4);
    
  //// INaCa
  i_NaCa      = kNaCa*(exp(gamma*y[0]*F/(R*T))*pow((float)y[17],3.0f)*Cao-exp((gamma-1.0)*y[0]*F/(R*T))*pow((float)Nao,3.0f)*y[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*y[0]*F/(R*T))));
    
  //// INaK
  i_NaK       = PNaK*Ko/(Ko+Km_K)*y[17]/(y[17]+Km_Na)/(1.0+0.1245*exp(-0.1*y[0]*F/(R*T))+0.0353*exp(-y[0]*F/(R*T)));
    
  //// IpCa
  i_PCa       = g_PCa*y[2]/(y[2]+KPCa);
    
  //// Background currents
  i_b_Na      = g_b_Na*(y[0]-E_Na);
    
  i_b_Ca      = g_b_Ca*(y[0]-E_Ca);
   
  //// Sarcoplasmic reticulum
  i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)y[2],2.0f));
    
  i_leak      = (y[1]-y[2])*V_leak;
    
  dydt[3]    = 0;
    
  // RyR
    
  RyRSRCass   = (1 - 1/(1 +  exp((y[1]-0.3)/0.1)));
  i_rel       = g_irel_max*RyRSRCass*y[21]*y[22]*(y[1]-y[2]);
    
  RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*y[2]-(RyRahalf))/0.0082));
  dydt[20]   = (RyRainfss- y[20])/RyRtauadapt;
    
  RyRoinfss   = (1 - 1/(1 +  exp((1000*y[2]-(y[20]+ RyRohalf))/0.003)));
  if (RyRoinfss>= y[21])
    RyRtauact = 18.75e-3;       //s
  else
    RyRtauact = 0.1*18.75e-3;   //s
    
  dydt[21]    = (RyRoinfss- y[21])/(RyRtauact);
    
  RyRcinfss   = (1/(1 + exp((1000*y[2]-(y[20]+RyRchalf))/0.001)));
  if (RyRcinfss>= y[22])
    RyRtauinact = 2*87.5e-3;    //s
  else
    RyRtauinact = 87.5e-3;      //s
    
  dydt[22]    = (RyRcinfss- y[22])/(RyRtauinact);
    
    
    
    
  //// Ca2+ buffering
  Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(y[2]+Kbuf_C), 2.0f));
  Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(y[1]+Kbuf_SR), 2.0f));
    
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
    
  //// Membrane potential
  dydt[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);
  //printf("dydt[0] = %.9f\n", dydt[0]);

  //Modification of activation by Martijn de Jong, 02-03-2022
  if (celltype2){
    if (fmod(current_time,pacing_interval) <= pacing_duration/10)
      dydt[0] += fmod(current_time,pacing_interval)/(pacing_duration/10)*pacing_strength/Cm;
    else if (fmod(current_time,pacing_interval) > pacing_duration/10 && fmod(current_time,pacing_interval) <= 9*pacing_duration/10)
      dydt[0] += pacing_strength/Cm;
    else if (fmod(current_time,pacing_interval) > 9*pacing_duration/10 && fmod(current_time,pacing_interval) <= pacing_duration)
        dydt[0] += (1-(fmod(current_time,pacing_duration/10)-9*pacing_duration/10)/(pacing_duration/10))*pacing_strength/Cm;
  }  
  
}


__device__ void derivsFitzHughNagumo(PDEFIELD_TYPE current_time, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id, PDEFIELD_TYPE interval_beats, PDEFIELD_TYPE pulse_duration, PDEFIELD_TYPE pulse_strength,  PDEFIELD_TYPE a, PDEFIELD_TYPE b, PDEFIELD_TYPE tau){
  PDEFIELD_TYPE RIext = 0;
  if (fmod(current_time, interval_beats) < pulse_duration && celltype2)
    RIext = pulse_strength;

  dydt[0] = y[0] - pow(y[0],3)/3 - y[1] + RIext;
  dydt[1] = y[0]/tau + a/tau - b*y[1] / tau; 
}


__device__ void RungeKuttaStep(PDEFIELD_TYPE* y, PDEFIELD_TYPE *dydt, int layers, PDEFIELD_TYPE thetime, PDEFIELD_TYPE stepsize, PDEFIELD_TYPE* yout, PDEFIELD_TYPE *yerr, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id){/*
  //Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x, use
  //the fifth-order Dormand-Prince Runge-Kutta method to advance the solution over an interval h
  //and return the incremented variables as yout[1..n]. Also return an estimate of the local
  //truncation error in yout using the embedded fourth-order method. The user supplies the routine
  //derivs(t,y,dydt,celltype2,pacing_interval,pacing_duration,,acing_strength), which returns derivatives dydt at t.
  int i;
  static PDEFIELD_TYPE 
  c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,a21=0.2,a31=3.0/40.0,
  a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
  a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
  a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
  a71=35.0/384.0,a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
  a76=11.0/84.0,e1=71.0/57600.0,e3=-71.0/16695.0,e4=71.0/1920.0,
  e5=-17253.0/339200.0,e6=22.0/525.0,e7=-1.0/40.0;


  PDEFIELD_TYPE k2[ARRAY_SIZE];
  PDEFIELD_TYPE k3[ARRAY_SIZE];
  PDEFIELD_TYPE k4[ARRAY_SIZE];
  PDEFIELD_TYPE k5[ARRAY_SIZE];
  PDEFIELD_TYPE k6[ARRAY_SIZE];
  PDEFIELD_TYPE ytemp[ARRAY_SIZE];
  PDEFIELD_TYPE dydtnew[ARRAY_SIZE];
  for (i=0;i<layers;i++) //First step.
    ytemp[i]=y[i]+a21*stepsize*dydt[i];

  derivsFitzHughNagumo(thetime+c2*stepsize,ytemp,k2,celltype2,pacing_interval,pacing_duration,pacing_strength,id);// Second step.
  for (i=0;i<layers;i++)
    ytemp[i]=y[i]+stepsize*(a31*dydt[i]+a32*k2[i]);
  derivsFitzHughNagumo(thetime+c3*stepsize,ytemp,k3,celltype2,pacing_interval,pacing_duration,pacing_strength,id); //Third step.
  for (i=0;i<layers;i++)
    ytemp[i]=y[i]+stepsize*(a41*dydt[i]+a42*k2[i]+a43*k3[i]);
  derivsFitzHughNagumo(thetime+c4*stepsize,ytemp,k4,celltype2,pacing_interval,pacing_duration,pacing_strength,id); //Fourth step.
  for (i=0;i<layers;i++)
    ytemp[i]=y[i]+stepsize*(a51*dydt[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
  derivsFitzHughNagumo(thetime+c5*stepsize,ytemp,k5,celltype2,pacing_interval,pacing_duration,pacing_strength,id); //Fifth step.
  for (i=0;i<layers;i++)
    ytemp[i]=y[i]+stepsize*(a61*dydt[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
  PDEFIELD_TYPE timeplusdt = thetime+stepsize;
  derivsFitzHughNagumo(timeplusdt,ytemp,k6,celltype2,pacing_interval,pacing_duration,pacing_strength,id); //Sixth step.
  for (i=0;i<layers;i++) //Accumulate increments with proper weights.
    yout[i]=y[i]+stepsize*(a71*dydt[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
  derivsFitzHughNagumo(timeplusdt,yout,dydtnew,celltype2,pacing_interval,pacing_duration,pacing_strength,id);
  for (i=0;i<layers;i++) //Estimate error as difference between fourth- and fifth-order methods.
    yerr[i]=stepsize*(e1*dydt[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydtnew[i]);
*/}

__device__ void StepsizeControl(PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, int layers, PDEFIELD_TYPE *thetime, PDEFIELD_TYPE stepsize_try, PDEFIELD_TYPE eps, PDEFIELD_TYPE* yscal, PDEFIELD_TYPE* stepsize_did, PDEFIELD_TYPE* stepsize_next, PDEFIELD_TYPE dt, PDEFIELD_TYPE stepsize_min, bool overshot, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id){
  /* Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
  adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
  at the starting value of the independent variable x. Also input are the stepsize to be attempted
  htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
  scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
  actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied
  routine that computes the right-hand side derivatives. */
  int i;
  PDEFIELD_TYPE err,stepsize; //stepsize_temp;
  PDEFIELD_TYPE yerr[ARRAY_SIZE];
  PDEFIELD_TYPE ytemp[ARRAY_SIZE];
  PDEFIELD_TYPE scale,scaling_factor,maxy;;

  const PDEFIELD_TYPE alpha = 0.2;
  const PDEFIELD_TYPE Safety = 0.9;
  const PDEFIELD_TYPE minscale = 0.2;
  const PDEFIELD_TYPE maxscale = 10;
  const PDEFIELD_TYPE rtol = 1e-3;
  const PDEFIELD_TYPE atol = 1e-6;

  stepsize=stepsize_try; // Set stepsize to the initial trial value.
  //stepsize=0.000001;
  for(;;){
    RungeKuttaStep(y,dydt,layers,*thetime,stepsize,ytemp,yerr,celltype2,pacing_interval,pacing_duration,pacing_strength,id); // Take a step.  
    if (id == 2100){
      printf("At index %i, we attempt a step of %.10f nanoseconds and the time is %.9f\n", id, stepsize*1e9, *thetime);
    }

    err=0; //Evaluate accuracy.
    for (i=0;i<layers;i++){
      //compute the total euclidean scaled error
      maxy = fabs(y[i]);
      if (fabs(y[i]) < fabs(ytemp[i]))
        maxy = fabs(ytemp[i]);
      scaling_factor = atol+rtol*maxy;
      err += pow(yerr[i]/scaling_factor,2);
    }
    err = sqrt(err/layers); 
    err /= eps; // Scale relative to required tolerance.
    err = 0;
    if (err <= 1.0){  
      break; //Step succeeded. Compute size of next step.
    }

    scale=fmax(Safety*pow(err,-alpha),minscale);
    stepsize *= scale;
    

  }
  scale=Safety*pow(err,-alpha);
  if (scale<minscale) scale = minscale;
  if (scale>maxscale) scale = maxscale;
  *stepsize_next = stepsize*scale;
  *thetime += (*stepsize_did=stepsize);
  for (i=0;i<layers;i++) {
    y[i]=ytemp[i];
  }  
}

__global__ void ODEstepRKA(PDEFIELD_TYPE dt, PDEFIELD_TYPE thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype, PDEFIELD_TYPE* next_stepsize, PDEFIELD_TYPE stepsize_min, PDEFIELD_TYPE eps, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength){
  /* Ordinary Differential Equation step Runge Kutta Adaptive
  Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
  adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
  at the starting value of the independent variable x. Also input are the stepsize to be attempted
  htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
  scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
  actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied
  routine that computes the right-hand side derivatives. */
  

  PDEFIELD_TYPE begin_time,stepsize_next,stepsize_did,stepsize, end_time;
  PDEFIELD_TYPE yscal[ARRAY_SIZE];
  PDEFIELD_TYPE y[ARRAY_SIZE];
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE current_time;
  PDEFIELD_TYPE MaxTimeError = 1e-9;
  PDEFIELD_TYPE stepsize_overshot;
  bool overshot = false;
  bool celltype2 = false;
  int i;

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    if (celltype[id] < 1){
      for (i = 0; i < layers; i++) //fill with current PDE values
        alt_PDEvars[i*sizex*sizey + id]= PDEvars[i*sizex*sizey + id];
    }
    
    else{
      celltype2 = false; 
      if (celltype[id] == 2)
        celltype2 = true;
      begin_time = thetime;
      current_time = thetime;
      end_time = thetime + dt;
      stepsize=next_stepsize[id];
      for (i=0;i<layers;i++) 
        y[i]=PDEvars[i*sizex*sizey + id];
      while(fabs(current_time - begin_time - dt)>MaxTimeError){

        overshot = false;
        //derivsFitzHughNagumo(current_time,y,dydt,celltype2, pacing_interval,pacing_duration,pacing_strength,id);

        if (stepsize+current_time > end_time){
          stepsize_overshot = stepsize; 
          stepsize=end_time - current_time;// If stepsize can overshoot, decrease.
          overshot = true;
        }
        StepsizeControl(y,dydt,layers,&current_time,stepsize,eps,yscal,&stepsize_did,&stepsize_next, dt, stepsize_min, overshot, celltype2, pacing_interval,pacing_duration,pacing_strength,id);
        if (fabs(current_time - begin_time - dt)<MaxTimeError) { //Are we done?
          for (i=0;i<layers;i++) {
            alt_PDEvars[i*sizex*sizey + id]=y[i];
          }
          if(overshot && (stepsize_overshot < dt)) 
            next_stepsize[id] = stepsize_overshot;
          else
            next_stepsize[id] = stepsize_next;
        }
        else 
          stepsize = stepsize_next;
      }
    }
  }
}

__global__ void ODEstepFE(PDEFIELD_TYPE dt, PDEFIELD_TYPE thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype, PDEFIELD_TYPE* next_stepsize, PDEFIELD_TYPE stepsize_min, PDEFIELD_TYPE eps, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, PDEFIELD_TYPE interval_beats, PDEFIELD_TYPE pulse_duration, PDEFIELD_TYPE pulse_strength,  PDEFIELD_TYPE a, PDEFIELD_TYPE b, PDEFIELD_TYPE tau){
  /* Ordinary Differential Equation step Runge Kutta Adaptive
  Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
  adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
  at the starting value of the independent variable x. Also input are the stepsize to be attempted
  htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
  scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
  actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied
  routine that computes the right-hand side derivatives. */
  



  PDEFIELD_TYPE dtt = dt;
  PDEFIELD_TYPE begin_time,stepsize_next,stepsize_did,stepsize, end_time;
  PDEFIELD_TYPE yscal[ARRAY_SIZE];
  PDEFIELD_TYPE y[ARRAY_SIZE];
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE current_time;
  PDEFIELD_TYPE MaxTimeError = 1e-5;
  PDEFIELD_TYPE stepsize_overshot;
  bool overshot = false;
  bool celltype2 = false;
  int i;

  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    if (celltype[id] < 1){
      for (i = 0; i < layers; i++) //fill with current PDE values
        alt_PDEvars[i*sizex*sizey + id]= PDEvars[i*sizex*sizey + id];
    }
    
    else{
      celltype2 = false; 
      if (celltype[id] == 2)
        celltype2 = true;
      begin_time = thetime;
      current_time = thetime;
      end_time = thetime + dt;
      stepsize=next_stepsize[id];
      for (i=0;i<layers;i++) 
        y[i]=PDEvars[i*sizex*sizey + id];
      while(current_time - begin_time - dt<0){


        overshot = false;
        derivsFitzHughNagumo(current_time,y,dydt,celltype2, pacing_interval,pacing_duration,pacing_strength, id, interval_beats, pulse_duration, pulse_strength,  a, b, tau);
        if (id == 23885)
          printf("dydt[0] = %.10f, dydt[1] = %.10f, y[0] = %.10f, y[1] = %.10f \n", dydt[0], dydt[1], y[0], y[1]);

        current_time += dtt;
        if (fabs(current_time - begin_time - dt) < MaxTimeError) { //Are we done?
          for (i=0;i<layers;i++) {
            alt_PDEvars[i*sizex*sizey + id]=y[i]+dydt[i]*dtt; 
          }
        }
        else{   
          for (i=0;i<layers;i++) {
            y[i]=y[i]+dydt[i]*dtt;  
          }
        }
      }
    }
  }
}


__global__ void ForwardEulerStep(PDEFIELD_TYPE dt, PDEFIELD_TYPE thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype){
  
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int i;

  


  
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE y[ARRAY_SIZE];

  //Declare variables needed for Paci2020 model and assign the constants

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
    
      //// Nernst potential
  PDEFIELD_TYPE E_Na;
  PDEFIELD_TYPE E_Ca;
  PDEFIELD_TYPE E_K;
  PDEFIELD_TYPE PkNa = 0.03;   // dimensionless (in electric_potentials)
  PDEFIELD_TYPE E_Ks;
    
  //// INa adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_Na = 3671.2302; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaFRedMed)*6447.1896;
  PDEFIELD_TYPE i_Na;
    
  PDEFIELD_TYPE m_inf;
  PDEFIELD_TYPE tau_m;
    
  PDEFIELD_TYPE h_inf;
  PDEFIELD_TYPE tau_h;
    
  PDEFIELD_TYPE j_inf;
  PDEFIELD_TYPE tau_j;
    
    
  //// INaL
  PDEFIELD_TYPE myCoefTauM  = 1;
  PDEFIELD_TYPE tauINaL = 200; //ms
  PDEFIELD_TYPE GNaLmax = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
  PDEFIELD_TYPE Vh_hLate = 87.61;
  PDEFIELD_TYPE i_NaL;
    
  PDEFIELD_TYPE m_inf_L;
  PDEFIELD_TYPE alpha_m_L;
  PDEFIELD_TYPE beta_m_L;
  PDEFIELD_TYPE tau_m_L;
    
  PDEFIELD_TYPE h_inf_L;
  PDEFIELD_TYPE tau_h_L = 1 * tauINaL;
    
  //// If adapted from DOI:10.3389/fphys.2018.00080
  PDEFIELD_TYPE g_f = 1; //((time<tDrugApplication)*1+(time >= tDrugApplication)*IfRedMed)*22.2763088;
  PDEFIELD_TYPE fNa = 0.37;
  PDEFIELD_TYPE fK = 1 - fNa;
  PDEFIELD_TYPE i_fK;
  PDEFIELD_TYPE i_fNa;
  PDEFIELD_TYPE i_f;
    
  PDEFIELD_TYPE Xf_infinity;
  PDEFIELD_TYPE tau_Xf; 
    
      //// ICaL
  PDEFIELD_TYPE g_CaL = 8.635702e-5;   // metre_cube_per_F_per_s (in i_CaL)
  PDEFIELD_TYPE i_CaL;  
  PDEFIELD_TYPE precision = 0.0001;     
    
  PDEFIELD_TYPE d_infinity;
  PDEFIELD_TYPE alpha_d;
  PDEFIELD_TYPE beta_d;
  PDEFIELD_TYPE gamma_d;
  PDEFIELD_TYPE tau_d;
    
  PDEFIELD_TYPE f1_inf;
  PDEFIELD_TYPE constf1;
    
  PDEFIELD_TYPE tau_f1;
    
  PDEFIELD_TYPE f2_inf;
  PDEFIELD_TYPE constf2 = 1.0;
  PDEFIELD_TYPE tau_f2;
    
  PDEFIELD_TYPE alpha_fCa;
  PDEFIELD_TYPE beta_fCa;
  PDEFIELD_TYPE gamma_fCa;
  PDEFIELD_TYPE fCa_inf;
    
  PDEFIELD_TYPE constfCa;
    
  PDEFIELD_TYPE tau_fCa     = 0.002;   // second (in i_CaL_fCa_gate)
    
  //// Ito
  PDEFIELD_TYPE g_to = 29.9038;   // S_per_F (in i_to)  
  PDEFIELD_TYPE i_to;
    
  PDEFIELD_TYPE q_inf;
  PDEFIELD_TYPE tau_q;
    
    
  PDEFIELD_TYPE r_inf;
  PDEFIELD_TYPE tau_r;
    
  //// IKs
  PDEFIELD_TYPE g_Ks = 2.041;   // S_per_F (in i_Ks)
  PDEFIELD_TYPE i_Ks; // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
  PDEFIELD_TYPE Xs_infinity;
  PDEFIELD_TYPE alpha_Xs;
  PDEFIELD_TYPE beta_Xs;
  PDEFIELD_TYPE tau_Xs;
    
  //// IKr
  PDEFIELD_TYPE L0 = 0.025;   // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE Q = 2.3;     // dimensionless (in i_Kr_Xr1_gate)
  PDEFIELD_TYPE g_Kr = 29.8667;   // S_per_F (in i_Kr)
  PDEFIELD_TYPE i_Kr;
    
  PDEFIELD_TYPE V_half;
    
  PDEFIELD_TYPE Xr1_inf;
  PDEFIELD_TYPE alpha_Xr1;
  PDEFIELD_TYPE beta_Xr1;
  PDEFIELD_TYPE tau_Xr1;
    
  PDEFIELD_TYPE Xr2_infinity;
  PDEFIELD_TYPE alpha_Xr2;
  PDEFIELD_TYPE beta_Xr2;
  PDEFIELD_TYPE tau_Xr2;
    
  //// IK1
  PDEFIELD_TYPE alpha_K1;
  PDEFIELD_TYPE beta_K1;
  PDEFIELD_TYPE XK1_inf;
  PDEFIELD_TYPE g_K1 = 28.1492;   // S_per_F (in i_K1)
  PDEFIELD_TYPE i_K1;
    
  //// INaCa
  PDEFIELD_TYPE KmCa = 1.38;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE KmNai = 87.5;   // millimolar (in i_NaCa)
  PDEFIELD_TYPE Ksat = 0.1;    // dimensionless (in i_NaCa)
  PDEFIELD_TYPE gamma = 0.35;   // dimensionless (in i_NaCa)
  PDEFIELD_TYPE alpha = 2.16659;
  PDEFIELD_TYPE kNaCa = 3917.0463; //((time<tDrugApplication)*1+(time >= tDrugApplication)*INaCaRedMed) * 6514.47574;   // A_per_F (in i_NaCa)
  PDEFIELD_TYPE i_NaCa;

  //// INaK
  PDEFIELD_TYPE Km_K = 1.0;    // millimolar (in i_NaK)
  PDEFIELD_TYPE Km_Na = 40.0;   // millimolar (in i_NaK)
  PDEFIELD_TYPE PNaK = 2.74240;// A_per_F (in i_NaK)
  PDEFIELD_TYPE i_NaK;
    
  //// IpCa
  PDEFIELD_TYPE KPCa = 0.0005;   // millimolar (in i_PCa)
  PDEFIELD_TYPE g_PCa = 0.4125;   // A_per_F (in i_PCa)
  PDEFIELD_TYPE i_PCa;
    
  //// Background currents
  PDEFIELD_TYPE g_b_Na = 1.14;         // S_per_F (in i_b_Na)
  PDEFIELD_TYPE i_b_Na;
    
  PDEFIELD_TYPE g_b_Ca = 0.8727264;    // S_per_F (in i_b_Ca)
  PDEFIELD_TYPE i_b_Ca;

  PDEFIELD_TYPE i_up;
  PDEFIELD_TYPE i_leak;
    
  //// Sarcoplasmic reticulum
  PDEFIELD_TYPE VmaxUp = 0.82205;
  PDEFIELD_TYPE Kup	= 4.40435e-4;
    
  PDEFIELD_TYPE V_leak = 4.48209e-4;
    
  // RyR
  PDEFIELD_TYPE g_irel_max = 55.808061;
  PDEFIELD_TYPE RyRa1 = 0.05169;
  PDEFIELD_TYPE RyRa2 = 0.050001;
  PDEFIELD_TYPE RyRahalf = 0.02632;
  PDEFIELD_TYPE RyRohalf = 0.00944;
  PDEFIELD_TYPE RyRchalf = 0.00167;
    
  PDEFIELD_TYPE RyRSRCass;
  PDEFIELD_TYPE i_rel;
    
  PDEFIELD_TYPE RyRainfss;
  PDEFIELD_TYPE RyRtauadapt = 1; //s
    
  PDEFIELD_TYPE RyRoinfss;
  PDEFIELD_TYPE RyRtauact;
    
  PDEFIELD_TYPE RyRcinfss;
  PDEFIELD_TYPE RyRtauinact;

  //// Ca2+ buffering
  PDEFIELD_TYPE Buf_C = 0.25;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Buf_SR = 10.0;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_C = 0.001;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Kbuf_SR = 0.3;   // millimolar (in calcium_dynamics)
  PDEFIELD_TYPE Cai_bufc;
  PDEFIELD_TYPE Ca_SR_bufSR;
    
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
  
  for (int id = index; id < sizex*sizey; id += stride){
    if (celltype[id] < 1){
      for (int l = 0; l < layers; l++) //fill with current PDE values
        alt_PDEvars[l*sizex*sizey + id]= PDEvars[l*sizex*sizey + id];
    }
    else{   
      for (int l = 0; l < layers; l++) //fill with current PDE values
        y[l] = PDEvars[l*sizex*sizey + id];
      //-----FIRST STEP ------------------------------------------------------------------------------------------------------------------------------------------------------------------  
      
        //// Nernst potential
      E_Na = R*T/F*log(Nao/y[17]);
      E_Ca = 0.5*R*T/F*log(Cao/y[2]);
      E_K  = R*T/F*log(Ko/Ki);
      E_Ks = R*T/F*log((Ko+PkNa*Nao)/(Ki+PkNa*y[17]));
    
      //// INa adapted from DOI:10.3389/fphys.2018.00080
      i_Na        =  g_Na*pow((float)y[13],3.0f)*y[11]*y[12]*(y[0] - E_Na);
    
      m_inf       = 1 / (1 + exp((y[0]*1000 + 39)/-11.2));
      tau_m       = (0.00001 + 0.00013*exp(-pow((float)((y[0]*1000 + 48)/15),2.0f)) + 0.000045 / (1 + exp((y[0]*1000 + 42)/-5)));
      dydt[13]   = (m_inf-y[13])/tau_m;
    
      h_inf       = 1 / (1 + exp((y[0]*1000 + 66.5)/6.8));
      tau_h       = (0.00007 + 0.034 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(y[0]*1000 + 79)/14)));
      dydt[11]   = (h_inf-y[11])/tau_h;
    
      j_inf       = h_inf;
      tau_j       = 10*(0.0007 + 0.15 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(y[0]*1000 + 79)/14)));
      dydt[12]   = (j_inf-y[12])/tau_j;
    
    
      //// INaL
      tauINaL     = 200; //ms
      GNaLmax     = 17.25;//((time<tDrugApplication)*1+(time >= tDrugApplication)*INaLRedMed)* 2.3*7.5; //(S/F)
      Vh_hLate    = 87.61;
      i_NaL       = GNaLmax* pow((float)y[18],3.0f)*y[19]*(y[0]-E_Na);
    
      m_inf_L     = 1/(1+exp(-(y[0]*1000+42.85)/(5.264)));
      alpha_m_L   = 1/(1+exp((-60-y[0]*1000)/5));
      beta_m_L    = 0.1/(1+exp((y[0]*1000+35)/5))+0.1/(1+exp((y[0]*1000-50)/200));
      tau_m_L     = 1 * myCoefTauM*alpha_m_L*beta_m_L;
      dydt[18]   = (m_inf_L-y[18])/tau_m_L*1000;
    
      h_inf_L     = 1/(1+exp((y[0]*1000+Vh_hLate)/(7.488)));
      tau_h_L     = 1 * tauINaL;
      dydt[19]   = (h_inf_L-y[19])/tau_h_L*1000;
    
      //// If adapted from DOI:10.3389/fphys.2018.00080
      i_fK        = fK*g_f*y[14]*(y[0] - E_K);
      i_fNa       = fNa*g_f*y[14]*(y[0] - E_Na);
      i_f         = i_fK + i_fNa;
    
      Xf_infinity = 1.0/(1.0 + exp((y[0]*1000 + 69)/8));
      tau_Xf      = 5600 / (1 + exp((y[0]*1000 + 65)/7) + exp(-(y[0]*1000 + 65)/19));
      dydt[14]   = 1000*(Xf_infinity-y[14])/tau_Xf;
      
    
      //// ICaL
      //Prevent division by 0
      if(y[0]< precision && y[0] > -precision) //hopital
        i_CaL =  g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T) *y[4]*y[5]*y[6]*y[7] / (2.0*F/(R*T)) * (y[2] - 0.341*Cao);
      else
        i_CaL = g_CaL*4.0*y[0]*pow((float)F,2.0f)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7]; //((time<tDrugApplication)*1+(time >= tDrugApplication)*ICaLRedMed)*g_CaL*4.0*y[0]*pow(F,2.0)/(R*T)*(y[2]*exp(2.0*y[0]*F/(R*T))-0.341*Cao)/(exp(2.0*y[0]*F/(R*T))-1.0)*y[4]*y[5]*y[6]*y[7];
    
      d_infinity  = 1.0/(1.0+exp(-(y[0]*1000.0+9.1)/7.0));
      alpha_d     = 0.25+1.4/(1.0+exp((-y[0]*1000.0-35.0)/13.0));
      beta_d      = 1.4/(1.0+exp((y[0]*1000.0+5.0)/5.0));
      gamma_d     = 1.0/(1.0+exp((-y[0]*1000.0+50.0)/20.0));
      tau_d       = (alpha_d*beta_d+gamma_d)*1.0/1000.0;
      dydt[4]    = (d_infinity-y[4])/tau_d;
    
      f1_inf      = 1.0/(1.0+exp((y[0]*1000.0+26.0)/3.0));
      if (f1_inf-y[5] > 0.0)
          constf1 = 1.0+1433.0*(y[2]-50.0*1.0e-6);
      else
          constf1 = 1.0;
    
      tau_f1      = (20.0+1102.5*exp(-pow((float)((y[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-y[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf1/1000.0;
      dydt[5]    = (f1_inf-y[5])/tau_f1;
    
      f2_inf      = 0.33+0.67/(1.0+exp((y[0]*1000.0+32.0)/4.0));
      tau_f2      = (600.0*exp(-pow((float)(y[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-y[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf2/1000.0;
      dydt[6]    = (f2_inf-y[6])/tau_f2;
    
      alpha_fCa   = 1.0/(1.0+pow((float)(y[2]/0.0006),8.0f));
      beta_fCa    = 0.1/(1.0+exp((y[2]-0.0009)/0.0001));
      gamma_fCa   = 0.3/(1.0+exp((y[2]-0.00075)/0.0008));
      fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
      if ((y[0] > -0.06) && (fCa_inf > y[7]))
          constfCa = 0.0;
      else
          constfCa = 1.0;
    
      dydt[7]    = constfCa*(fCa_inf-y[7])/tau_fCa;
    
      //// Ito
      i_to        = g_to*(y[0]-E_K)*y[15]*y[16];
    
      q_inf       = 1.0/(1.0+exp((y[0]*1000.0+53.0)/13.0));
      tau_q       = (6.06+39.102/(0.57*exp(-0.08*(y[0]*1000.0+44.0))+0.065*exp(0.1*(y[0]*1000.0+45.93))))/1000.0;
      dydt[15]   = (q_inf-y[15])/tau_q;
    
    
      r_inf       = 1.0/(1.0+exp(-(y[0]*1000.0-22.3)/18.75));
      tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(y[0]*1000.0+30.61))+0.369*exp(-0.12*(y[0]*1000.0+23.84))))/1000.0;
      dydt[16]   = (r_inf-y[16])/tau_r;
    
      //// IKs
      i_Ks        = g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
      Xs_infinity = 1.0/(1.0+exp((-y[0]*1000.0-20.0)/16.0));
      alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-y[0]*1000.0)/6.0));
      beta_Xs     = 1.0/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
      tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
      dydt[10]   = (Xs_infinity-y[10])/tau_Xs;
    
      //// IKr
      i_Kr         = g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4);
    
      V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
      Xr1_inf      = 1.0/(1.0+exp((V_half-y[0]*1000.0)/4.9));
      alpha_Xr1    = 450.0/(1.0+exp((-45.0-y[0]*1000.0)/10.0));
      beta_Xr1     = 6.0/(1.0+exp((30.0+y[0]*1000.0)/11.5));
      tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
      dydt[8]     = (Xr1_inf-y[8])/tau_Xr1;
    
      Xr2_infinity = 1.0/(1.0+exp((y[0]*1000.0+88.0)/50.0));
      alpha_Xr2    = 3.0/(1.0+exp((-60.0-y[0]*1000.0)/20.0));
      beta_Xr2     = 1.12/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
      tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
      dydt[9]    = (Xr2_infinity-y[9])/tau_Xr2;
    
      //// IK1
      alpha_K1    = 3.91/(1.0+exp(0.5942*(y[0]*1000.0-E_K*1000.0-200.0)));
      beta_K1     = (-1.509*exp(0.0002*(y[0]*1000.0-E_K*1000.0+100.0))+exp(0.5886*(y[0]*1000.0-E_K*1000.0-10.0)))/(1.0+exp(0.4547*(y[0]*1000.0-E_K*1000.0)));
      XK1_inf     = alpha_K1/(alpha_K1+beta_K1);
      i_K1        = g_K1*XK1_inf*(y[0]-E_K)*sqrt(Ko/5.4);
    
      //// INaCa
      i_NaCa      = kNaCa*(exp(gamma*y[0]*F/(R*T))*pow((float)y[17],3.0f)*Cao-exp((gamma-1.0)*y[0]*F/(R*T))*pow((float)Nao,3.0f)*y[2]*alpha)/((pow((float)KmNai,3.0f)+pow((float)Nao,3.0f))*(KmCa+Cao)*(1.0+Ksat*exp((gamma-1.0)*y[0]*F/(R*T))));
    
      //// INaK
      i_NaK       = PNaK*Ko/(Ko+Km_K)*y[17]/(y[17]+Km_Na)/(1.0+0.1245*exp(-0.1*y[0]*F/(R*T))+0.0353*exp(-y[0]*F/(R*T)));
    
      //// IpCa
      i_PCa       = g_PCa*y[2]/(y[2]+KPCa);
    
      //// Background currents
      i_b_Na      = g_b_Na*(y[0]-E_Na);
    
      i_b_Ca      = g_b_Ca*(y[0]-E_Ca);
    
      //// Sarcoplasmic reticulum
      i_up        = VmaxUp/(1.0+pow((float)Kup,2.0f)/pow((float)y[2],2.0f));
    
      i_leak      = (y[1]-y[2])*V_leak;
    
      dydt[3]    = 0;
    
      // RyR
    
      RyRSRCass   = (1 - 1/(1 +  exp((y[1]-0.3)/0.1)));
      i_rel       = g_irel_max*RyRSRCass*y[21]*y[22]*(y[1]-y[2]);
    
      RyRainfss   = RyRa1-RyRa2/(1 + exp((1000*y[2]-(RyRahalf))/0.0082));
      dydt[20]   = (RyRainfss- y[20])/RyRtauadapt;
    
      RyRoinfss   = (1 - 1/(1 +  exp((1000*y[2]-(y[20]+ RyRohalf))/0.003)));
      if (RyRoinfss>= y[21])
        RyRtauact = 18.75e-3;       //s
      else
        RyRtauact = 0.1*18.75e-3;   //s
    
      dydt[21]    = (RyRoinfss- y[21])/(RyRtauact);
    
      RyRcinfss   = (1/(1 + exp((1000*y[2]-(y[20]+RyRchalf))/0.001)));
      if (RyRcinfss>= y[22])
        RyRtauinact = 2*87.5e-3;    //s
      else
        RyRtauinact = 87.5e-3;      //s
    
      dydt[22]    = (RyRcinfss- y[22])/(RyRtauinact);
    
    
    
    
      //// Ca2+ buffering
      Cai_bufc    = 1.0/(1.0+Buf_C*Kbuf_C/pow((float)(y[2]+Kbuf_C), 2.0f));
      Ca_SR_bufSR = 1.0/(1.0+Buf_SR*Kbuf_SR/pow((float)(y[1]+Kbuf_SR), 2.0f));
    
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
    
      //// Membrane potential
      dydt[0] = -(i_K1+i_to+i_Kr+i_Ks+i_CaL+i_NaK+i_Na+i_NaL+i_NaCa+i_PCa+i_f+i_b_Na+i_b_Ca-i_stim);

      //-----WRITE NEW VALUES ------------------------------------------------------------------------------------------------------------------------------------------------------------------      
      #pragma unroll
      for (i=0;i<layers;i++) //Accumulate increments with proper weights.
        alt_PDEvars[i*sizex*sizey + id]=y[i]+dydt[i]*dt;
    } 

  }
}




__global__ void CopyAltToOriginalPDEvars(int sizex, int sizey, int layers, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < layers*sizex*sizey; id += stride){
    PDEvars[id] = alt_PDEvars[id]; 
  }
}

__global__ void CopyOriginalToAltPDEvars(int sizex, int sizey, int layers, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < layers*sizex*sizey; id += stride){
    alt_PDEvars[id] = PDEvars[id]; 
  }
}

__global__ void FlipSigns(int sizex, int sizey, int layers, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < layers*sizex*sizey; id += stride){
    alt_PDEvars[id] = -PDEvars[id]; 
  }
}


void PDE::cuPDEsteps(CellularPotts * cpm, int repeat){
  //copy current couplingcoefficient matrix and celltype matrix from host to device
  couplingcoefficient = cpm->getCouplingCoefficient();
  //int** cellnumber = cpm -> getSigma(); 
  cudaError_t errSync;
  cudaError_t errAsync;
  celltype = cpm->getTau();
  cudaMemcpy(d_couplingcoefficient, couplingcoefficient[0], sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  cudaMemcpy(d_celltype, celltype[0], sizex*sizey*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_PDEvars, PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);

  //setup matrices for upperdiagonal, diagonal and lower diagonal for both the horizontal and vertical direction, since these remain the same during once MCS
  InitializeDiagonals<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, lowerH, upperH, diagH, lowerV, upperV, diagV, d_couplingcoefficient);
  cudaDeviceSynchronize();
  errSync  = cudaGetLastError();
  errAsync = cudaDeviceSynchronize();
  if (errSync != cudaSuccess) 
    printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
  if (errAsync != cudaSuccess)
    printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));

  for (int iteration = 0; iteration < repeat; iteration++){
    //Do an ODE step of size dt/2
    ODEstepFE<<<par.number_of_cores, par.threads_per_core>>>(dt/2, thetime, layers, sizex, sizey, d_PDEvars, d_alt_PDEvars, d_celltype, next_stepsize, min_stepsize, par.eps, pacing_interval, par.pacing_duration, par.pacing_strength, par.FHN_interval_beats, par.FHN_pulse_duration, par.FHN_pulse_strength, par.FHN_a, par.FHN_b, par.FHN_tau);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));

    cudaMemcpy(alt_PDEvars, d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    for (int i=0; i < sizex*sizey; i++){
      if (!(alt_PDEvars[i]>-100000 && alt_PDEvars[i] < 100000)){
        cout << "Error at i = " << i << ". Abort.\n";
        exit(1);
      }
    }
    cout << "After first FE step, alt_PDEvars[23885] = " << alt_PDEvars[23885] << endl;


    //Do a vertical ADI sweep of size dt/2
    InitializeVerticalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BV, d_couplingcoefficient, d_alt_PDEvars);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
    statusV = cusparseSgtsvInterleavedBatch(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, pbufferV);
    if (statusV != CUSPARSE_STATUS_SUCCESS)
    {
      cout << statusV << endl;
    }
    cudaDeviceSynchronize();
    NewPDEfieldV0<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, BV, d_PDEvars); //////
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
    NewPDEfieldOthers<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, layers, BV, d_PDEvars, d_alt_PDEvars); //////
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync)); 
    

    cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    cout << "After first ADI step, PDEvars[23885] = " << PDEvars[23885] << endl;

    //increase time by dt/2
    thetime = thetime + dt/2;  
    //Do an ODE step of size dt/2
    ODEstepFE<<<par.number_of_cores, par.threads_per_core>>>(dt/2, thetime, layers, sizex, sizey, d_PDEvars, d_alt_PDEvars, d_celltype, next_stepsize, min_stepsize, par.eps, pacing_interval, par.pacing_duration, par.pacing_strength, par.FHN_interval_beats, par.FHN_pulse_duration, par.FHN_pulse_strength, par.FHN_a, par.FHN_b, par.FHN_tau);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));   
      
      
    cudaMemcpy(alt_PDEvars, d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    cout << "After second FE step, alt_PDEvars[23885] = " << alt_PDEvars[23885] << endl;

    //Do a horizontal ADI sweep of size dt/2
    InitializeHorizontalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BH, d_couplingcoefficient, d_alt_PDEvars);
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
    statusH = cusparseSgtsvInterleavedBatch(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, pbufferH);
    if (statusH != CUSPARSE_STATUS_SUCCESS)
    {
      cout << statusH << endl;
    }
    NewPDEfieldH0<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, BH, d_PDEvars);    
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
    NewPDEfieldOthers<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, layers, BV, d_PDEvars, d_alt_PDEvars); //////
    errSync  = cudaGetLastError();
      errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));  
    
      //increase time by dt/2
    thetime = thetime + dt/2; 


    cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    cout << "After first ADI step, PDEvars[23885] = " << PDEvars[23885] << endl;

    
    cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    ofstream myfile;
    myfile.open ("data.txt", std::ios_base::app);
    myfile << thetime << ",";
    for (int i = 0; i < layers; i++)
      myfile << PDEvars[23885 + i*sizex*sizey] << ",";
    myfile << endl;
    cout << "PDEvars[" << int(sizex/2*sizey + 0.5 * sizey) << "] = " << PDEvars[int(sizex/2*sizey + 0.5 * sizey)] << endl;
    if (!(PDEvars[int(sizex/2*sizey + 0.5 * sizey) ]>-100000 && PDEvars[int(sizex/2*sizey + 0.5 * sizey)] < 100000)){
      cout << "We encountered a NaN error. Abort the program. \n";
      exit(1);
    }
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
	  sum+=PDEvars[l*sizex*sizey + (x+1)*sizey+y];
	  sum+=PDEvars[l*sizex*sizey + (x-1)*sizey+y];
	  sum+=PDEvars[l*sizex*sizey + x*sizey+y+1];
	  sum+=PDEvars[l*sizex*sizey + x*sizey+y-1];
	  sum-=4*PDEvars[l*sizex*sizey + x*sizey+y];
	  alt_PDEvars[l*sizex*sizey + x*sizey+y]=PDEvars[l*sizex*sizey + x*sizey+y]+sum*dt*par.diff_coeff[l]/dx2;
      }
    }
    PDEFIELD_TYPE *tmp;
    tmp=PDEvars;
    PDEvars=alt_PDEvars;
    alt_PDEvars=tmp;
  
    thetime+=dt;
  }
}

/*double PDE::GetChemAmount(const int layer) {
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
}*/

// private
void PDE::NoFluxBoundaries(void) {
  // all gradients at the edges become zero, 
  // so nothing flows out
  // Note that four corners points are not defined (0.)
  // but they aren't used in the calculations
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l*sizex*sizey + x*sizey+0]=PDEvars[l*sizex*sizey + x*sizey+1];
      PDEvars[l*sizex*sizey + x*sizey+sizey-1]=PDEvars[l*sizex*sizey + x*sizey+sizey-2];
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l*sizex*sizey + 0*sizey+y]=PDEvars[l*sizex*sizey + 1*sizey+y];
      PDEvars[l*sizex*sizey + (sizex-1)*sizey+y]=PDEvars[l*sizex*sizey + (sizex-2)*sizey+y];
    }
  }
}


// private
void PDE::AbsorbingBoundaries(void) {
  // all boundaries are sinks, 
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l*sizex*sizey + x*sizey+0]=0.;
      PDEvars[l*sizex*sizey + x*sizey+sizey-1]=0.;
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l*sizex*sizey + 0*sizey+y]=0.;
      PDEvars[l*sizex*sizey + (sizex-1)*sizey+y]=0.;
    }
  }
}

// private
void PDE::PeriodicBoundaries(void) {
  // periodic...
  for (int l=0;l<layers;l++) {
    for (int x=0;x<sizex;x++) {
      PDEvars[l*sizex*sizey + x*sizey]=PDEvars[l*sizex*sizey + x*sizey+sizey-2];
      PDEvars[l*sizex*sizey + x*sizey+sizey-1]=PDEvars[l*sizex*sizey + x+sizey+1];
    }
    for (int y=0;y<sizey;y++) {
      PDEvars[l*sizex*sizey + y]=PDEvars[l*sizex*sizey + (sizex-2)*sizey+y];
      PDEvars[l*sizex*sizey + (sizex-1)+y]=PDEvars[l*sizex*sizey + 1*sizey+y];
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
      PDEvars[first_grad_layer*sizex*sizey + x*sizey+y]=(PDEvars[layer*sizex*sizey+(x+1)*sizey+y]-PDEvars[layer*sizex*sizey + (x-1)*sizey+y])/2.;
    } 
  }
  // GradY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      PDEvars[(first_grad_layer+1)*sizex*sizey + x*sizey+y]=(PDEvars[layer*sizex*sizey+x*sizey+y+1]-PDEvars[layer*sizex*sizey+x*sizey+y-1])/2.;
    } 
  }
  // GradXX
  for (int y=0;y<sizey;y++) {
    for (int x=1;x<sizex-1;x++) {
      PDEvars[(first_grad_layer+2)*sizex*sizey + x*sizey+y]=PDEvars[layer*sizex*sizey+ (x+1)*sizey+y]-PDEvars[layer*sizex*sizey + (x-1)*sizey+y]-2*PDEvars[layer*sizex*sizey+x*sizey+y];
    } 
  }

  // GradYY
  for (int x=0;x<sizex;x++) {
    for (int y=1;y<sizey-1;y++) {
      PDEvars[(first_grad_layer+3)*sizex*sizey + x*sizey+y]=PDEvars[layer*sizex*sizey + x*sizey-1]-PDEvars[layer*sizex*sizey+x*sizey+y+1]-2*PDEvars[layer*sizex*sizey + x*sizey+y];
    } 
  }
}

void PDE::PlotVectorField(Graphics &g, int stride, int linelength, int first_grad_layer) {
  // Plot vector field assuming it's in layer 1 and 2
  for (int x=1;x<sizex-1;x+=stride) {
    for (int y=1;y<sizey-1;y+=stride) {
      
      // calculate line
      int x1,y1,x2,y2;
      
      x1=(int)(x-linelength*PDEvars[(first_grad_layer)*sizex*sizey + x*sizey + y]);
      y1=(int)(y-linelength*PDEvars[(first_grad_layer+1)*sizex*sizey + x*sizey + y]);
      x2=(int)(x+linelength*PDEvars[(first_grad_layer)*sizex*sizey + x*sizey + y]);
      y2=(int)(y+linelength*PDEvars[(first_grad_layer+1)*sizex*sizey + x*sizey + y]);
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
  double val = PDEvars[layer*sizex*sizey+x*sizey+y];
  if (val > -100){
    graphics->Rectangle(MapColour(val), x, y);
    return false;
  }
  return true;
}


void PDE::SetSpeciesName(int l, const char *name) {
    species_names[l]=string(name);
}


/*void PDE::InitLinearYGradient(int spec, double conc_top, double conc_bottom) {
    for (int y=0;y<sizey;y++) {
      double val=(double)conc_top+y*((double)(conc_bottom-conc_top)/(double)sizey);
    for (int x=0;x<sizex;x++) {
      PDEvars[spec][x][y]=val;
    }
    cerr << y << " " << val << endl;
  }
}*/
