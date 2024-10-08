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
#include "random.hpp"
#include "ca.hpp"
#include "pde.hpp"
#include "conrec.hpp"
#include "graph.hpp"
#include <cusparse.h>
#include <cuda_profiler_api.h>
#define ARRAY_SIZE 33

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
  if (par.activation_times)
    Activation_times_array = new PDEFIELD_TYPE[sizex*sizey];
    Activation_times_array_written = new bool[sizex*sizey];
  if (par.SF_one_pixel)
    copy_for_SF = new PDEFIELD_TYPE[ARRAY_SIZE]; //for SF computation
  if (par.SF_all){
    SF_start_array = new bool[sizex*sizey];
    SF_end_array = new bool[sizex*sizey];
    SF_Q_tot_array = new PDEFIELD_TYPE[sizex*sizey];
  }
  SF_locator_one = par.sizey * par.SF_x + par.SF_y;
  min_stepsize = par.min_stepsize;
  btype = 1;
  dx2 = par.dx*par.dx;
  dt = par.dt;
  ddt = par.ddt;
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
    free(PDEvars);
  }
  if (alt_PDEvars) {
    free(alt_PDEvars);
  }
  if (FHN_a) {
    free(FHN_a);
  }
  if (FHN_b) {
    free(FHN_b);
  }
  if (FHN_tau) {
    free(FHN_tau);
  }
  if (couplingcoefficient) {
    free(couplingcoefficient);
  }
  if (sigmafield) {
    free(sigmafield);
  }
  if (celltype) {
    free(celltype);
  }
  cudaFree(d_PDEvars);
  cudaFree(d_alt_PDEvars);
  cudaFree(d_couplingcoefficient);
  cudaFree(d_celltype);
  cudaFree(d_sigmafield);
  cudaFree(d_FHN_a);
  cudaFree(d_FHN_b);
  cudaFree(d_FHN_tau);
  cudaFree(upperH);
  cudaFree(diagH);
  cudaFree(lowerH);
  cudaFree(BH);
  cudaFree(upperV);
  cudaFree(diagV);
  cudaFree(lowerV);
  cudaFree(BV);
  cudaFree(next_stepsize);

}



void PDE::InitializePDEvars(CellularPotts *cpm, int* celltypes){
  Successful_activation = false;
  PDEFIELD_TYPE PDEinit1[ARRAY_SIZE];
  PDEFIELD_TYPE PDEinit2[ARRAY_SIZE];
  int* mask = cpm->getMask()[0];
   //For Paci2018
  /*PDEinit[0] = -0.070;
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

  //Paci2018 at rest
  /*
  PDEinit[0] = -0.0793716;
  PDEinit[1] = 0.399439; 
  PDEinit[2] = 2.23913e-05;
  PDEinit[3] = 0;
  PDEinit[4] = 4.36816e-05; 
  PDEinit[5] = 0.954067;
  PDEinit[6] = 0.999986; 
  PDEinit[7] = 0.99867;
  PDEinit[8] = 0.00524384;
  PDEinit[9] = 0.456945;
  PDEinit[10] = 0.0238964;
  PDEinit[11] = 0.868928;
  PDEinit[12] = 0.843866; 
  PDEinit[13] = 0.0264785; 
  PDEinit[14] = 0.193978; 
  PDEinit[15] = 0.883421;
  PDEinit[16] = 0.00439947;
  PDEinit[17] = 9.1711; 
  PDEinit[18] = 0.000969319;
  PDEinit[19] = 0.210587; 
  PDEinit[20] = 0.136443; 
  PDEinit[21] = 2.52234e-44; 
  PDEinit[22] = 0.997029;*/

  
  //For FHN
  /*
  PDEinit[0] = FHN_0;
  PDEinit[1] = FHN_1;
  */

  //Maleckar2009
  
  PDEinit1[0] = -74.031982;
  PDEinit1[1] = 130.022096;
  PDEinit1[2] = 8.516766;
  PDEinit1[3] = 0.003289;
  PDEinit1[4] =  0.877202; 
  PDEinit1[5] = 0.873881;
  PDEinit1[6] = 7.1e-5;
  PDEinit1[7] = 0.000014;
  PDEinit1[8] = 0.998597;
  PDEinit1[9] = 0.998586;
  PDEinit1[10] = 5.560224;
  PDEinit1[11] = 129.485991;
  PDEinit1[12] = 0.001089; 
  PDEinit1[13] = 0.948597;
  PDEinit1[14] = 0.000367;
  PDEinit1[15] = 0.96729;
  PDEinit1[16] = 0.004374;
  PDEinit1[17] = 0.000053;
  PDEinit1[18] = 1.815768;
  PDEinit1[19] = 6.5e-5;
  PDEinit1[20] = 0.026766;
  PDEinit1[21] = 0.012922;
  PDEinit1[22] = 0.190369;
  PDEinit1[23] = 0.714463;
  PDEinit1[24] = 1.38222;
  PDEinit1[25] = 0.632613;
  PDEinit1[26] = 0.649195;
  PDEinit1[27] = 0.431547;
  PDEinit1[28] = 0.470055;
  PDEinit1[29] = 0.002814;
  for (int j = 30; j < ARRAY_SIZE; j++)
    PDEinit1[j] = 0;
  
  /*
  //Grandi2011
  PDEinit1[0] = -7.34336366728778671e+01;
  PDEinit1[1] =  3.94923428392655786e-03;
  PDEinit1[2] =  1.35538532457244482e-01;
  PDEinit1[3] =  1.03674364292988680e-01;
  PDEinit1[4] =  1.90759804527589089e-01;
  PDEinit1[5] =  1.35640688636079511e-02;
  PDEinit1[6] =  2.14063418881809235e-02;
  PDEinit1[7] =  4.45327242854324807e-03;
  PDEinit1[8] =  1.27856586024588575e-01;
  PDEinit1[9] =  5.69999505293381902e-03;
  PDEinit1[10] =  1.83143535034222225e-02;
  PDEinit1[11] =  2.10808768153058460e-04;
  PDEinit1[12] =  3.25814677291117296e-04;
  PDEinit1[13] =  2.33018340557575125e-04;
  PDEinit1[14] =  3.61396062660070427e+00;
  PDEinit1[15] =  7.88607791910409195e-01;
  PDEinit1[16] =  9.15153381546177336e+00;
  PDEinit1[17] =  9.15182798281732346e+00;
  PDEinit1[18] =  5.02305826642838293e-01;
  PDEinit1[19] =  1.13337536953687845e+00;
  PDEinit1[20] =  7.02128101897185673e-04;
  PDEinit1[21] =  2.16850216379767157e-05;
  PDEinit1[22] =  9.98384427312367095e-01;
  PDEinit1[23] =  4.49572164109603364e-02;
  PDEinit1[24] =  3.28512098597005947e-02;
  PDEinit1[25] = 120.0;
  PDEinit1[26] =  1.31290096227093382e-03;
  PDEinit1[27] =  7.49436760722081534e-03;
  PDEinit1[28] =  9.15199678386256998e+00;
  PDEinit1[29] =  3.93548562883350357e-04;
  PDEinit1[30] =  9.58234428284286399e-01;
  PDEinit1[31] =  3.15482710277587786e-01;
  PDEinit1[32] =  2.48034071360795916e-01;
  PDEinit1[33] =  1.89326933812916480e-02;
  PDEinit1[34] =  3.79829335413739144e-02;
  PDEinit1[35] =  1.01974216400706526e-02;
  PDEinit1[36] =  1.37939236359928058e-03;
  PDEinit1[37] =  9.45874848392074696e-01;
  PDEinit1[38] =  5.01323282772066123e-07;
  PDEinit1[39] =  2.01567245823636694e-06;
  PDEinit1[40] =  8.00819151705148946e-01;
  for (int j = 41; j < ARRAY_SIZE; j++)
    PDEinit1[j] = 0;
  */
  //Fabbri-Severi 2017
  PDEinit2[0] = -47.787168;
  PDEinit2[1] = 6.226104e-5; 
  PDEinit2[2] = 5;
  PDEinit2[3] = 0.009508;
  PDEinit2[4] = 0.447724;
  PDEinit2[5] = 0.003058;
  PDEinit2[6] = 0.001921;
  PDEinit2[7] = 0.846702;
  PDEinit2[8] = 0.844449;
  PDEinit2[9] = 0.268909;
  PDEinit2[10] = 0.020484;
  PDEinit2[11] = 0.9308;
  PDEinit2[12] = 6.181512e-9; 
  PDEinit2[13] = 4.595622e-10; 
  PDEinit2[14] = 0.069199; 
  PDEinit2[15] = 0.409551;
  PDEinit2[16] = 0.435148;
  PDEinit2[17] = 9.15641e-6;
  PDEinit2[18] = 0.653777;
  PDEinit2[19] = 0.217311; 
  PDEinit2[20] = 0.158521;
  PDEinit2[21] = 0.017929; 
  PDEinit2[22] = 0.259947;
  PDEinit2[23] = 0.138975;
  PDEinit2[24] = 0.011845;
  PDEinit2[25] = 0.845304;
  PDEinit2[26] = 0.430836;
  PDEinit2[27] = 0.014523;
  PDEinit2[28] = 0.283185;
  PDEinit2[29] = 0.011068;
  PDEinit2[30] = 0.709051;
  PDEinit2[31] = 0.1162;
  PDEinit2[32] = 0.00277;
  for (int j = 33; j < ARRAY_SIZE; j++)
    PDEinit2[j] = 0;

  

  for (int layer = 0; layer<ARRAY_SIZE; layer++){
    for (int i = layer*sizex*sizey; i<(layer+1)*sizex*sizey; i++){
      if (celltypes[i-layer*sizey*sizex] == 1){
        PDEvars[i] = PDEinit1[layer];
      }
      else if (celltypes[i-layer*sizey*sizex] == 2){
        PDEvars[i] = PDEinit2[layer];
      }
      else{ 
        PDEvars[i] = 0;
      }
      //if (layer == 0 && (i%(sizey) == 200))
      //  PDEvars[i] = 100000000;
    }
  }
  /*cout << "Celltype of interest = " << celltypes[int(406.5 * 410)] << endl;
  for (int k = 0; k < ARRAY_SIZE; k++)
    cout << "pos[" << k  << "] = " << PDEvars[int(406.5 * 410)+k*sizex*sizey] << endl;*/
  
  


}

void PDE::InitializeFHNvarsCells(int nr_cells, PDEFIELD_TYPE* FHN_a, PDEFIELD_TYPE* FHN_b, PDEFIELD_TYPE* FHN_tau, PDEFIELD_TYPE FHN_a_var, PDEFIELD_TYPE FHN_b_var, PDEFIELD_TYPE FHN_tau_var, PDEFIELD_TYPE FHN_a_base, PDEFIELD_TYPE FHN_b_base, PDEFIELD_TYPE FHN_tau_base){
  

  PDEFIELD_TYPE multiplication_factor;

  //randomly choose the base FHN_variables with a factor in the range (1-FHN_x_var, 1+FHN_x_var)
  for (int i = 0; i < nr_cells; i++){
    multiplication_factor = RANDOM()*2*FHN_a_var+1-FHN_a_var;
    FHN_a[i] = FHN_a_base * multiplication_factor;
    multiplication_factor = RANDOM()*2*FHN_b_var+1-FHN_b_var;
    FHN_b[i] = FHN_b_base * multiplication_factor;
    multiplication_factor = RANDOM()*2*FHN_tau_var+1-FHN_tau_var;
    FHN_tau[i] = FHN_tau_base * multiplication_factor;
  }
  
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

  // number of contouring levelss
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
  PDEFIELD_TYPE decay_rate = par.decay_rate[0];
  PDEFIELD_TYPE secr_rate = par.secr_rate[0];
  
  
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

void PDE::InitializeCuda(CellularPotts *cpm, int n_init_cells){
  cout << "Start cuda init" << endl;
  //AllocateTridiagonalvars(sizex, sizey);

  cudaMalloc((void**) &d_couplingcoefficient, sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMalloc((void**) &d_celltype, sizex*sizey*sizeof(int));
  cudaMalloc((void**) &d_sigmafield, sizex*sizey*sizeof(int));

  cudaMalloc((void**) &d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMemcpy(d_PDEvars, PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  cudaMalloc((void**) &d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMemcpy(d_alt_PDEvars, alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);


  //Needed for ADI steps
  gpuErrchk(cudaMallocManaged(&upperH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&diagH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&lowerH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&BH, sizex*sizey*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&upperV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&diagV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&lowerV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&BV, sizey*sizex*sizeof(PDEFIELD_TYPE)));
  gpuErrchk(cudaMallocManaged(&next_stepsize, sizey*sizex*sizeof(PDEFIELD_TYPE)));

  handleH = 0;
  pbuffersizeH = 0;
  pbufferH = NULL;
  statusH=cusparseCreate(&handleH);
  #ifdef PDEFIELD_DOUBLE
    cusparseDgtsvInterleavedBatch_bufferSizeExt(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, &pbuffersizeH); //Compute required buffersize for horizontal sweep
  #else
    cusparseSgtsvInterleavedBatch_bufferSizeExt(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, &pbuffersizeH); //Compute required buffersize for horizontal sweep
  #endif
  gpuErrchk(cudaMalloc( &pbufferH, sizeof(char)* pbuffersizeH));
  

  handleV = 0;
  pbuffersizeV = 0;
  pbufferV = NULL;
  statusV=cusparseCreate(&handleV);
  #ifdef PDEFIELD_DOUBLE
    cusparseDgtsvInterleavedBatch_bufferSizeExt(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, &pbuffersizeV); //Compute required buffersize for vertical sweep
  #else
    cusparseSgtsvInterleavedBatch_bufferSizeExt(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, &pbuffersizeV); //Compute required buffersize for vertical sweep
  #endif
  gpuErrchk(cudaMalloc( &pbufferV, sizeof(char)* pbuffersizeV));

  
  //Needed for FHN variation
  //cout << "par.n_init_cells = " << n_init_cells << endl;
  //cout << "FHN_a[0] = " << FHN_a[0] << endl;
  //cudaMalloc((void**) &d_FHN_a, n_init_cells*sizeof(PDEFIELD_TYPE));
  //cudaMalloc((void**) &d_FHN_b, n_init_cells*sizeof(PDEFIELD_TYPE));
  //cudaMalloc((void**) &d_FHN_tau, n_init_cells*sizeof(PDEFIELD_TYPE));
  //cudaMemcpy(d_FHN_a, FHN_a, n_init_cells*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  //cudaMemcpy(d_FHN_b, FHN_b, n_init_cells*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  //cudaMemcpy(d_FHN_tau, FHN_tau, n_init_cells*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
  cout << "End cuda init" << endl;

}

__global__ void InitializeLastStepsize(PDEFIELD_TYPE min_stepsize, PDEFIELD_TYPE* next_stepsize, int sizex, int sizey){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    next_stepsize[id] = min_stepsize;
  }
}


void PDE::InitializeActivationTimes(void){
  for (int i = 0; i < par.sizex*par.sizey; i++){
    Activation_times_array[i] = 0;
    Activation_times_array_written[i] = false;
  }
  cudaMalloc((void**) &d_Activation_times_array, sizex*sizey*sizeof(PDEFIELD_TYPE));
  cudaMemcpy(d_Activation_times_array, Activation_times_array, sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice);
}

void PDE::InitializeSFComputation(CellularPotts *cpm){
  cout << "start init" << endl;
  int** celltypes = cpm->getTau();
  cout << "Point 2" << endl;
  for (int i = 0; i < par.sizex*par.sizey; i++){
    SF_Q_tot_array [i] = 0;
    //cout << "i = " << i << " and therefore, x = " << i/par.sizey << " and y = " << i % par.sizey << endl;
    if (celltypes[0][i] == 1){
      //cout << "i = " << i << " and therefore, x = " << i/par.sizey << " and y = " << i % par.sizey << endl;
      SF_start_array[i] = false;
      SF_end_array[i] = false;
    }
    else{
      //cout << "i = " << i << " and therefore, x = " << i/par.sizey << " and y = " << i % par.sizey << endl;
      SF_start_array[i] = true;
      SF_end_array[i] = true;
    }
  }
  cout << "end init" << endl;

}

void PDE::InitializePDEs(CellularPotts *cpm){
  int** celltypes = cpm->getTau();
  //InitializePDEvars(cpm, par.FHN_start_0, par.FHN_start_1);
  InitializePDEvars(cpm, celltypes[0]);
  InitializeCuda(cpm, par.n_init_cells+1);
  if (par.activation_times){
    InitializeActivationTimes();
  }
  if (par.SF_all){
    InitializeSFComputation(cpm);
  }
  //InitializeLastStepsize<<<par.number_of_cores, par.threads_per_core>>>(min_stepsize, next_stepsize, sizex, sizey);
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



__device__ void derivsPaci(PDEFIELD_TYPE current_time, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id){

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

__device__ void derivsPaci_RL(PDEFIELD_TYPE current_time, PDEFIELD_TYPE dt, PDEFIELD_TYPE* y, PDEFIELD_TYPE* y_new, PDEFIELD_TYPE* dydt, bool celltype2, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id){

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
  y_new[13] = m_inf + (y[13] - m_inf) * exp(-dt/tau_m);
    
  h_inf       = 1 / (1 + exp((y[0]*1000 + 66.5)/6.8));
  tau_h       = (0.00007 + 0.034 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.0002 / (1 + exp(-(y[0]*1000 + 79)/14)));
  dydt[11]   = (h_inf-y[11])/tau_h;
  y_new[11] = h_inf + (y[11] - h_inf) * exp(-dt/tau_h);
    
  j_inf       = h_inf;
  tau_j       = 10*(0.0007 + 0.15 / (1 + exp((y[0]*1000 + 41)/5.5) + exp(-(y[0]*1000 + 41)/14)) + 0.002 / (1 + exp(-(y[0]*1000 + 79)/14)));
  dydt[12]   = (j_inf-y[12])/tau_j;
  y_new[12] = j_inf + (y[12] - j_inf) * exp(-dt/tau_j);
    
    
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
  y_new[18] = m_inf_L + (y[18] - m_inf_L) * exp(-dt/(tau_m_L/1000));
    
  h_inf_L     = 1/(1+exp((y[0]*1000+Vh_hLate)/(7.488)));
  tau_h_L     = 1 * tauINaL;
  dydt[19]   = (h_inf_L-y[19])/tau_h_L*1000;
  y_new[19] = h_inf_L + (y[19] - h_inf_L) * exp(-dt/(tau_h_L/1000));
    
  //// If adapted from DOI:10.3389/fphys.2018.00080
  i_fK        = fK*g_f*y[14]*(y[0] - E_K);
  i_fNa       = fNa*g_f*y[14]*(y[0] - E_Na);
  i_f         = i_fK + i_fNa;
    
  Xf_infinity = 1.0/(1.0 + exp((y[0]*1000 + 69)/8));
  tau_Xf      = 5600 / (1 + exp((y[0]*1000 + 65)/7) + exp(-(y[0]*1000 + 65)/19));
  dydt[14]   = 1000*(Xf_infinity-y[14])/tau_Xf;
  y_new[14] = Xf_infinity + (y[14] - Xf_infinity) * exp(-dt/(tau_Xf/1000));
      
    
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
  y_new[4] = d_infinity + (y[4] - d_infinity) * exp(-dt/tau_d);
    
  f1_inf      = 1.0/(1.0+exp((y[0]*1000.0+26.0)/3.0));
  if (f1_inf-y[5] > 0.0)
      constf1 = 1.0+1433.0*(y[2]-50.0*1.0e-6);
  else
      constf1 = 1.0;
  
  tau_f1      = (20.0+1102.5*exp(-pow((float)((y[0]*1000.0+27.0)/15.0),2.0f))+200.0/(1.0+exp((13.0-y[0]*1000.0)/10.0))+180.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf1/1000.0;
  dydt[5]    = (f1_inf-y[5])/tau_f1;
  y_new[5] = f1_inf + (y[5] - f1_inf) * exp(-dt/tau_f1);
   
  f2_inf      = 0.33+0.67/(1.0+exp((y[0]*1000.0+32.0)/4.0));
  tau_f2      = (600.0*exp(-pow((float)(y[0]*1000.0+25.0),2.0f)/170.0)+31.0/(1.0+exp((25.0-y[0]*1000.0)/10.0))+16.0/(1.0+exp((30.0+y[0]*1000.0)/10.0)))*constf2/1000.0;
  dydt[6]    = (f2_inf-y[6])/tau_f2;
  y_new[6] = f2_inf + (y[6] - f2_inf) * exp(-dt/tau_f2);
   
  alpha_fCa   = 1.0/(1.0+pow((float)(y[2]/0.0006),8.0f));
  beta_fCa    = 0.1/(1.0+exp((y[2]-0.0009)/0.0001));
  gamma_fCa   = 0.3/(1.0+exp((y[2]-0.00075)/0.0008));
  fCa_inf     = (alpha_fCa+beta_fCa+gamma_fCa)/1.3156;
    
  if ((y[0] > -0.06) && (fCa_inf > y[7]))
    constfCa = 0.0;
  else
    constfCa = 1.0;
    
  dydt[7]    = constfCa*(fCa_inf-y[7])/tau_fCa;
  y_new[7] = fCa_inf + (y[7] - fCa_inf) * exp(-dt/(tau_fCa/constfCa));
    
  //// Ito
  i_to        = g_to*(y[0]-E_K)*y[15]*y[16];
    
  q_inf       = 1.0/(1.0+exp((y[0]*1000.0+53.0)/13.0));
  tau_q       = (6.06+39.102/(0.57*exp(-0.08*(y[0]*1000.0+44.0))+0.065*exp(0.1*(y[0]*1000.0+45.93))))/1000.0;
  dydt[15]   = (q_inf-y[15])/tau_q;
  y_new[15] = q_inf + (y[15] - q_inf) * exp(-dt/tau_q);
    
  r_inf       = 1.0/(1.0+exp(-(y[0]*1000.0-22.3)/18.75));
  tau_r       = (2.75352+14.40516/(1.037*exp(0.09*(y[0]*1000.0+30.61))+0.369*exp(-0.12*(y[0]*1000.0+23.84))))/1000.0;
  dydt[16]   = (r_inf-y[16])/tau_r;
  y_new[16] = r_inf + (y[16] - r_inf) * exp(-dt/tau_r);
    
  //// IKs
  i_Ks        = g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0f)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f))); // ((time<tDrugApplication)*1+(time >= tDrugApplication)*IKsRedMed)*g_Ks*(y[0]-E_Ks)*pow((float)y[10],2.0)*(1.0+0.6/(1.0+pow((float)(3.8*0.00001/y[2]),1.4f)));
    
  Xs_infinity = 1.0/(1.0+exp((-y[0]*1000.0-20.0)/16.0));
  alpha_Xs    = 1100.0/sqrt(1.0+exp((-10.0-y[0]*1000.0)/6.0));
  beta_Xs     = 1.0/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
  tau_Xs      = 1.0*alpha_Xs*beta_Xs/1000.0;
  dydt[10]   = (Xs_infinity-y[10])/tau_Xs;
  y_new[10] = Xs_infinity + (y[10] - Xs_infinity) * exp(-dt/tau_Xs);
    
  //// IKr
  i_Kr         = g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4); //((time<tDrugApplication)*1+(time >= tDrugApplication)*IKrRedMed)*g_Kr*(y[0]-E_K)*y[8]*y[9]*sqrt(Ko/5.4);
    
  V_half       = 1000.0*(-R*T/(F*Q)*log(pow((float)(1.0+Cao/2.6),4.0f)/(pow((float)(1.0+Cao/0.58),4.0f)*L0))-0.019);
    
  Xr1_inf      = 1.0/(1.0+exp((V_half-y[0]*1000.0)/4.9));
  alpha_Xr1    = 450.0/(1.0+exp((-45.0-y[0]*1000.0)/10.0));
  beta_Xr1     = 6.0/(1.0+exp((30.0+y[0]*1000.0)/11.5));
  tau_Xr1      = 1.0*alpha_Xr1*beta_Xr1/1000.0;
  dydt[8]     = (Xr1_inf-y[8])/tau_Xr1;
  y_new[8] = Xr1_inf + (y[8] - Xr1_inf) * exp(-dt/tau_Xr1);
    
  Xr2_infinity = 1.0/(1.0+exp((y[0]*1000.0+88.0)/50.0));
  alpha_Xr2    = 3.0/(1.0+exp((-60.0-y[0]*1000.0)/20.0));
  beta_Xr2     = 1.12/(1.0+exp((-60.0+y[0]*1000.0)/20.0));
  tau_Xr2      = 1.0*alpha_Xr2*beta_Xr2/1000.0;
  dydt[9]    = (Xr2_infinity-y[9])/tau_Xr2;
  y_new[9] = Xr2_infinity + (y[9] - Xr2_infinity) * exp(-dt/tau_Xr2);
    
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

  //Modification of activation by Martijn de Jong, 02-03-2022
  if (celltype2){
    if (fmod(current_time,pacing_interval) <= pacing_duration/10)
      dydt[0] += fmod(current_time,pacing_interval)/(pacing_duration/10)*pacing_strength/Cm;
    else if (fmod(current_time,pacing_interval) > pacing_duration/10 && fmod(current_time,pacing_interval) <= 9*pacing_duration/10)
      dydt[0] += pacing_strength/Cm;
    else if (fmod(current_time,pacing_interval) > 9*pacing_duration/10 && fmod(current_time,pacing_interval) <= pacing_duration)
        dydt[0] += (1-(fmod(current_time,pacing_duration/10)-9*pacing_duration/10)/(pacing_duration/10))*pacing_strength/Cm;
  }
  
  for (int i = 0; i < 23; i++){
    if (!(i == 4  || i == 5 || i == 6 || i == 7 || i == 8 || i == 9 || i == 10 || i == 11 || i == 12 || i == 13 || i == 14 || i == 15 || i == 16 || i == 18 || i == 19)){
      y_new[i] = y[i]+dt*dydt[i]; 
    }
  }
}



__device__ void derivsFitzHughNagumo(PDEFIELD_TYPE current_time, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, bool celltype2, int* sigmafield, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength, int id, PDEFIELD_TYPE interval_beats, PDEFIELD_TYPE pulse_duration, PDEFIELD_TYPE pulse_strength,  PDEFIELD_TYPE a, PDEFIELD_TYPE b, PDEFIELD_TYPE tau, PDEFIELD_TYPE* FHN_a, PDEFIELD_TYPE* FHN_b, PDEFIELD_TYPE* FHN_tau ){
  
    a = 0.1;
    PDEFIELD_TYPE epsilon = 10;
    PDEFIELD_TYPE beta = -1.0;
    PDEFIELD_TYPE RIext = 0;
    PDEFIELD_TYPE c = 0.191;
    PDEFIELD_TYPE timescale = 40;
  
  int sigma = sigmafield[id];
  if (fmod(current_time, interval_beats) < pulse_duration && celltype2)
    RIext = pulse_strength;


  dydt[0] = timescale*(epsilon*(y[0]*(1-y[0])*(y[0]-beta)) - y[1] + RIext);
  dydt[1] = timescale*((y[0] -a*y[1]+c));      

}

__device__ void derivsFabbriSeveri(PDEFIELD_TYPE VOI, PDEFIELD_TYPE* STATES, PDEFIELD_TYPE* RATES, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE I_f_factor, PDEFIELD_TYPE I_Kr_factor, int id){
  /*
   There are a total of 101 entries in the algebraic variable array.
   There are a total of 33 entries in each of the rate and state variable arrays.
   There are a total of 116 entries in the constant variable array.
 */
  /*
  * VOI is time in component environment (second).
  * CONSTANTS_FS[0] is R in component Membrane (joule_per_kilomole_kelvin).
  * CONSTANTS_FS[1] is T in component Membrane (kelvin).
  * CONSTANTS_FS[2] is F in component Membrane (coulomb_per_mole).
  * CONSTANTS_FS[3] is C in component Membrane (microF).
  * CONSTANTS_FS[91] is RTONF in component Membrane (millivolt).
  * ALGEBRAIC[59] is i_f in component i_f (nanoA).
  * ALGEBRAIC[61] is i_NaK in component i_NaK (nanoA).
  * ALGEBRAIC[75] is i_NaCa in component i_NaCa (nanoA).
  * ALGEBRAIC[79] is i_Na in component i_Na (nanoA).
  * ALGEBRAIC[89] is i_Kr in component i_Kr (nanoA).
  * ALGEBRAIC[95] is i_Ks in component i_Ks (nanoA).
  * ALGEBRAIC[87] is i_to in component i_to (nanoA).
  * ALGEBRAIC[83] is i_CaL in component i_CaL (nanoA).
  * ALGEBRAIC[84] is i_CaT in component i_CaT (nanoA).
  * ALGEBRAIC[98] is i_KACh in component i_KACh (nanoA).
  * ALGEBRAIC[85] is i_Kur in component i_Kur (nanoA).
  * ALGEBRAIC[9] is V in component Membrane (millivolt).
  * CONSTANTS_FS[4] is clamp_mode in component Membrane (dimensionless).
  * ALGEBRAIC[5] is V_clamp in component Voltage_clamp (millivolt).
  * STATES[0] is V_ode in component Membrane (millivolt).
  * ALGEBRAIC[100] is i_tot in component Membrane (nanoA).
  * CONSTANTS_FS[5] is t_holding in component Voltage_clamp (second).
  * CONSTANTS_FS[6] is t_test in component Voltage_clamp (second).
  * CONSTANTS_FS[7] is V_test in component Voltage_clamp (millivolt).
  * CONSTANTS_FS[8] is V_holding in component Voltage_clamp (millivolt).
  * CONSTANTS_FS[9] is ACh in component Rate_modulation_experiments (millimolar).
  * CONSTANTS_FS[10] is Iso_1_uM in component Rate_modulation_experiments (dimensionless).
  * ALGEBRAIC[18] is Nai in component Nai_concentration (millimolar).
  * CONSTANTS_FS[11] is Nao in component Ionic_values (millimolar).
  * CONSTANTS_FS[12] is Ki in component Ionic_values (millimolar).
  * CONSTANTS_FS[13] is Ko in component Ionic_values (millimolar).
  * STATES[1] is Ca_sub in component Ca_dynamics (millimolar).
  * CONSTANTS_FS[14] is Cao in component Ionic_values (millimolar).
  * ALGEBRAIC[36] is E_Na in component Ionic_values (millivolt).
  * CONSTANTS_FS[96] is E_K in component Ionic_values (millivolt).
  * ALGEBRAIC[0] is E_Ca in component Ionic_values (millivolt).
  * CONSTANTS_FS[110] is V_sub in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[112] is V_i in component Cell_parameters (millimetre3).
  * ALGEBRAIC[49] is i_fNa in component i_f (nanoA).
  * ALGEBRAIC[82] is i_siNa in component i_CaL (nanoA).
  * STATES[2] is Nai_ in component Nai_concentration (millimolar).
  * CONSTANTS_FS[15] is Nai_clamp in component Nai_concentration (dimensionless).
  * ALGEBRAIC[55] is i_fK in component i_f (nanoA).
  * CONSTANTS_FS[16] is g_f in component i_f (microS).
  * CONSTANTS_FS[92] is G_f in component i_f (microS).
  * CONSTANTS_FS[103] is g_f_Na in component i_f (microS).
  * CONSTANTS_FS[100] is G_f_Na in component i_f (microS).
  * CONSTANTS_FS[101] is g_f_K in component i_f (microS).
  * CONSTANTS_FS[97] is G_f_K in component i_f (microS).
  * CONSTANTS_FS[17] is Km_f in component i_f (millimolar).
  * CONSTANTS_FS[18] is alpha in component i_f (dimensionless).
  * STATES[3] is y in component i_f_y_gate (dimensionless).
  * CONSTANTS_FS[19] is blockade in component i_f (dimensionless).
  * ALGEBRAIC[10] is tau_y in component i_f_y_gate (second).
  * ALGEBRAIC[29] is y_infinity in component i_f_y_gate (dimensionless).
  * CONSTANTS_FS[95] is ACh_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[99] is Iso_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[20] is y_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[21] is Km_Kp in component i_NaK (millimolar).
  * CONSTANTS_FS[22] is Km_Nap in component i_NaK (millimolar).
  * CONSTANTS_FS[23] is i_NaK_max in component i_NaK (nanoA).
  * CONSTANTS_FS[102] is Iso_increase in component i_NaK (dimensionless).
  * CONSTANTS_FS[24] is K_NaCa in component i_NaCa (nanoA).
  * ALGEBRAIC[72] is x1 in component i_NaCa (dimensionless).
  * ALGEBRAIC[68] is x2 in component i_NaCa (dimensionless).
  * ALGEBRAIC[73] is x3 in component i_NaCa (dimensionless).
  * ALGEBRAIC[74] is x4 in component i_NaCa (dimensionless).
  * ALGEBRAIC[63] is k41 in component i_NaCa (dimensionless).
  * CONSTANTS_FS[104] is k34 in component i_NaCa (dimensionless).
  * ALGEBRAIC[71] is k23 in component i_NaCa (dimensionless).
  * ALGEBRAIC[70] is k21 in component i_NaCa (dimensionless).
  * ALGEBRAIC[67] is k32 in component i_NaCa (dimensionless).
  * ALGEBRAIC[62] is k43 in component i_NaCa (dimensionless).
  * ALGEBRAIC[65] is k12 in component i_NaCa (dimensionless).
  * ALGEBRAIC[66] is k14 in component i_NaCa (dimensionless).
  * CONSTANTS_FS[25] is Qci in component i_NaCa (dimensionless).
  * CONSTANTS_FS[26] is Qn in component i_NaCa (dimensionless).
  * CONSTANTS_FS[27] is Qco in component i_NaCa (dimensionless).
  * CONSTANTS_FS[28] is K3ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[29] is Kci in component i_NaCa (millimolar).
  * CONSTANTS_FS[30] is K1ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[31] is K2ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[32] is Kcni in component i_NaCa (millimolar).
  * CONSTANTS_FS[33] is K3no in component i_NaCa (millimolar).
  * CONSTANTS_FS[34] is K1no in component i_NaCa (millimolar).
  * CONSTANTS_FS[35] is K2no in component i_NaCa (millimolar).
  * CONSTANTS_FS[36] is Kco in component i_NaCa (millimolar).
  * ALGEBRAIC[69] is do in component i_NaCa (dimensionless).
  * ALGEBRAIC[64] is di in component i_NaCa (dimensionless).
  * CONSTANTS_FS[37] is blockade_NaCa in component i_NaCa (dimensionless).
  * ALGEBRAIC[77] is i_Na_ in component i_Na (nanoA).
  * ALGEBRAIC[78] is i_Na_L in component i_Na (nanoA).
  * CONSTANTS_FS[38] is g_Na in component i_Na (microS).
  * CONSTANTS_FS[39] is g_Na_L in component i_Na (microS).
  * ALGEBRAIC[76] is E_mh in component i_Na (millivolt).
  * STATES[4] is m in component i_Na_m_gate (dimensionless).
  * STATES[5] is h in component i_Na_h_gate (dimensionless).
  * ALGEBRAIC[46] is alpha_m in component i_Na_m_gate (per_second).
  * ALGEBRAIC[52] is beta_m in component i_Na_m_gate (per_second).
  * ALGEBRAIC[11] is m_infinity in component i_Na_m_gate (dimensionless).
  * ALGEBRAIC[57] is tau_m in component i_Na_m_gate (second).
  * CONSTANTS_FS[40] is delta_m in component i_Na_m_gate (millivolt).
  * ALGEBRAIC[30] is E0_m in component i_Na_m_gate (millivolt).
  * ALGEBRAIC[31] is alpha_h in component i_Na_h_gate (per_second).
  * ALGEBRAIC[47] is beta_h in component i_Na_h_gate (per_second).
  * ALGEBRAIC[12] is h_infinity in component i_Na_h_gate (dimensionless).
  * ALGEBRAIC[53] is tau_h in component i_Na_h_gate (second).
  * ALGEBRAIC[80] is i_siCa in component i_CaL (nanoA).
  * ALGEBRAIC[81] is i_siK in component i_CaL (nanoA).
  * CONSTANTS_FS[106] is ACh_block in component i_CaL (dimensionless).
  * CONSTANTS_FS[41] is P_CaL in component i_CaL (nanoA_per_millimolar).
  * STATES[6] is dL in component i_CaL_dL_gate (dimensionless).
  * STATES[7] is fL in component i_CaL_fL_gate (dimensionless).
  * STATES[8] is fCa in component i_CaL_fCa_gate (dimensionless).
  * CONSTANTS_FS[105] is Iso_increase in component i_CaL (dimensionless).
  * ALGEBRAIC[13] is dL_infinity in component i_CaL_dL_gate (dimensionless).
  * ALGEBRAIC[60] is tau_dL in component i_CaL_dL_gate (second).
  * ALGEBRAIC[48] is alpha_dL in component i_CaL_dL_gate (per_second).
  * ALGEBRAIC[58] is beta_dL in component i_CaL_dL_gate (per_second).
  * ALGEBRAIC[32] is adVm in component i_CaL_dL_gate (millivolt).
  * ALGEBRAIC[54] is bdVm in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[42] is k_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[43] is V_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[107] is Iso_shift_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[108] is Iso_slope_dL in component i_CaL_dL_gate (dimensionless).
  * ALGEBRAIC[14] is fL_infinity in component i_CaL_fL_gate (dimensionless).
  * ALGEBRAIC[33] is tau_fL in component i_CaL_fL_gate (second).
  * CONSTANTS_FS[44] is shift_fL in component i_CaL_fL_gate (millivolt).
  * CONSTANTS_FS[45] is k_fL in component i_CaL_fL_gate (millivolt).
  * CONSTANTS_FS[46] is alpha_fCa in component i_CaL_fCa_gate (per_second).
  * ALGEBRAIC[1] is fCa_infinity in component i_CaL_fCa_gate (dimensionless).
  * ALGEBRAIC[7] is tau_fCa in component i_CaL_fCa_gate (second).
  * CONSTANTS_FS[47] is Km_fCa in component i_CaL_fCa_gate (millimolar).
  * CONSTANTS_FS[48] is P_CaT in component i_CaT (nanoA_per_millimolar).
  * STATES[9] is dT in component i_CaT_dT_gate (dimensionless).
  * STATES[10] is fT in component i_CaT_fT_gate (dimensionless).
  * ALGEBRAIC[15] is dT_infinity in component i_CaT_dT_gate (dimensionless).
  * ALGEBRAIC[34] is tau_dT in component i_CaT_dT_gate (second).
  * ALGEBRAIC[16] is fT_infinity in component i_CaT_fT_gate (dimensionless).
  * ALGEBRAIC[35] is tau_fT in component i_CaT_fT_gate (second).
  * CONSTANTS_FS[49] is offset_fT in component i_CaT_fT_gate (second).
  * ALGEBRAIC[86] is j_SRCarel in component Ca_SR_release (millimolar_per_second).
  * STATES[11] is R in component Ca_SR_release (dimensionless).
  * STATES[12] is O in component Ca_SR_release (dimensionless).
  * STATES[13] is I in component Ca_SR_release (dimensionless).
  * STATES[14] is RI in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[50] is ks in component Ca_SR_release (per_second).
  * CONSTANTS_FS[51] is MaxSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[52] is MinSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[53] is EC50_SR in component Ca_SR_release (millimolar).
  * CONSTANTS_FS[54] is HSR in component Ca_SR_release (dimensionless).
  * ALGEBRAIC[8] is koSRCa in component Ca_SR_release (per_millimolar2_second).
  * ALGEBRAIC[17] is kiSRCa in component Ca_SR_release (per_millimolar_second).
  * CONSTANTS_FS[55] is koCa in component Ca_SR_release (per_millimolar2_second).
  * CONSTANTS_FS[56] is kiCa in component Ca_SR_release (per_millimolar_second).
  * ALGEBRAIC[2] is kCaSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[57] is kim in component Ca_SR_release (per_second).
  * CONSTANTS_FS[58] is kom in component Ca_SR_release (per_second).
  * STATES[15] is Ca_jsr in component Ca_dynamics (millimolar).
  * ALGEBRAIC[3] is diff in component Ca_SR_release (millimolar).
  * ALGEBRAIC[4] is P_tot in component Ca_SR_release (dimensionless).
  * ALGEBRAIC[88] is j_Ca_dif in component Ca_intracellular_fluxes (millimolar_per_second).
  * ALGEBRAIC[91] is j_up in component Ca_intracellular_fluxes (millimolar_per_second).
  * ALGEBRAIC[94] is j_tr in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[59] is tau_dif_Ca in component Ca_intracellular_fluxes (second).
  * CONSTANTS_FS[60] is tau_tr in component Ca_intracellular_fluxes (second).
  * CONSTANTS_FS[98] is P_up in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[61] is P_up_basal in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[93] is b_up in component Ca_intracellular_fluxes (dimensionless).
  * CONSTANTS_FS[62] is K_up in component Ca_intracellular_fluxes (millimolar).
  * STATES[16] is Ca_nsr in component Ca_dynamics (millimolar).
  * STATES[17] is Cai in component Ca_dynamics (millimolar).
  * CONSTANTS_FS[63] is slope_up in component Ca_intracellular_fluxes (millimolar).
  * CONSTANTS_FS[64] is TC_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[65] is TMC_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[66] is CM_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[67] is CQ_tot in component Ca_buffering (millimolar).
  * ALGEBRAIC[93] is delta_fTC in component Ca_buffering (per_second).
  * ALGEBRAIC[96] is delta_fTMC in component Ca_buffering (per_second).
  * ALGEBRAIC[90] is delta_fCMs in component Ca_buffering (per_second).
  * ALGEBRAIC[99] is delta_fCMi in component Ca_buffering (per_second).
  * ALGEBRAIC[97] is delta_fCQ in component Ca_buffering (per_second).
  * ALGEBRAIC[6] is delta_fTMM in component Ca_buffering (per_second).
  * STATES[18] is fTMM in component Ca_buffering (dimensionless).
  * STATES[19] is fCMi in component Ca_buffering (dimensionless).
  * STATES[20] is fCMs in component Ca_buffering (dimensionless).
  * STATES[21] is fTC in component Ca_buffering (dimensionless).
  * STATES[22] is fTMC in component Ca_buffering (dimensionless).
  * STATES[23] is fCQ in component Ca_buffering (dimensionless).
  * CONSTANTS_FS[68] is kf_TC in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[69] is kf_TMM in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[70] is kf_TMC in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[71] is kf_CM in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[72] is kf_CQ in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[73] is kb_TC in component Ca_buffering (per_second).
  * CONSTANTS_FS[74] is kb_TMC in component Ca_buffering (per_second).
  * CONSTANTS_FS[75] is kb_TMM in component Ca_buffering (per_second).
  * CONSTANTS_FS[76] is kb_CM in component Ca_buffering (per_second).
  * CONSTANTS_FS[77] is kb_CQ in component Ca_buffering (per_second).
  * CONSTANTS_FS[78] is Mgi in component Ca_buffering (millimolar).
  * CONSTANTS_FS[111] is V_jsr in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[113] is V_nsr in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[109] is V_cell in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[79] is V_jsr_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[80] is V_i_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[81] is V_nsr_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[82] is R_cell in component Cell_parameters (micrometre).
  * CONSTANTS_FS[83] is L_cell in component Cell_parameters (micrometre).
  * CONSTANTS_FS[84] is L_sub in component Cell_parameters (micrometre).
  * CONSTANTS_FS[85] is g_Kur in component i_Kur (microS).
  * STATES[24] is r_Kur in component i_Kur_rKur_gate (dimensionless).
  * STATES[25] is s_Kur in component i_Kur_sKur_gate (dimensionless).
  * ALGEBRAIC[37] is tau_r_Kur in component i_Kur_rKur_gate (second).
  * ALGEBRAIC[19] is r_Kur_infinity in component i_Kur_rKur_gate (dimensionless).
  * ALGEBRAIC[38] is tau_s_Kur in component i_Kur_sKur_gate (second).
  * ALGEBRAIC[20] is s_Kur_infinity in component i_Kur_sKur_gate (dimensionless).
  * CONSTANTS_FS[86] is g_to in component i_to (microS).
  * STATES[26] is q in component i_to_q_gate (dimensionless).
  * STATES[27] is r in component i_to_r_gate (dimensionless).
  * ALGEBRAIC[21] is q_infinity in component i_to_q_gate (dimensionless).
  * ALGEBRAIC[39] is tau_q in component i_to_q_gate (second).
  * ALGEBRAIC[22] is r_infinity in component i_to_r_gate (dimensionless).
  * ALGEBRAIC[40] is tau_r in component i_to_r_gate (second).
  * CONSTANTS_FS[87] is g_Kr in component i_Kr (microS).
  * STATES[28] is paS in component i_Kr_pa_gate (dimensionless).
  * STATES[29] is paF in component i_Kr_pa_gate (dimensionless).
  * STATES[30] is piy in component i_Kr_pi_gate (dimensionless).
  * ALGEBRAIC[23] is pa_infinity in component i_Kr_pa_gate (dimensionless).
  * ALGEBRAIC[24] is alfapaF in component i_Kr_pa_gate (per_second).
  * ALGEBRAIC[25] is betapaF in component i_Kr_pa_gate (per_second).
  * ALGEBRAIC[41] is tau_paS in component i_Kr_pa_gate (second).
  * ALGEBRAIC[42] is tau_paF in component i_Kr_pa_gate (second).
  * ALGEBRAIC[43] is pi_infinity in component i_Kr_pi_gate (dimensionless).
  * ALGEBRAIC[26] is tau_pi in component i_Kr_pi_gate (second).
  * CONSTANTS_FS[94] is g_Ks in component i_Ks (microS).
  * CONSTANTS_FS[88] is g_Ks_ in component i_Ks (microS).
  * ALGEBRAIC[92] is E_Ks in component i_Ks (millivolt).
  * STATES[31] is n in component i_Ks_n_gate (dimensionless).
  * ALGEBRAIC[27] is n_infinity in component i_Ks_n_gate (dimensionless).
  * ALGEBRAIC[56] is tau_n in component i_Ks_n_gate (second).
  * CONSTANTS_FS[114] is Iso_shift in component i_Ks_n_gate (millivolt).
  * ALGEBRAIC[44] is alpha_n in component i_Ks_n_gate (per_second).
  * ALGEBRAIC[50] is beta_n in component i_Ks_n_gate (per_second).
  * CONSTANTS_FS[89] is g_KACh in component i_KACh (microS).
  * STATES[32] is a in component i_KACh_a_gate (dimensionless).
  * CONSTANTS_FS[90] is ACh_on in component i_KACh (dimensionless).
  * CONSTANTS_FS[115] is alpha_a in component i_KACh_a_gate (per_second).
  * ALGEBRAIC[28] is beta_a in component i_KACh_a_gate (per_second).
  * ALGEBRAIC[45] is a_infinity in component i_KACh_a_gate (dimensionless).
  * ALGEBRAIC[51] is tau_a in component i_KACh_a_gate (second).
  * RATES[0] is d/dt V_ode in component Membrane (millivolt).
  * RATES[2] is d/dt Nai_ in component Nai_concentration (millimolar).
  * RATES[3] is d/dt y in component i_f_y_gate (dimensionless).
  * RATES[4] is d/dt m in component i_Na_m_gate (dimensionless).
  * RATES[5] is d/dt h in component i_Na_h_gate (dimensionless).
  * RATES[6] is d/dt dL in component i_CaL_dL_gate (dimensionless).
  * RATES[7] is d/dt fL in component i_CaL_fL_gate (dimensionless).
  * RATES[8] is d/dt fCa in component i_CaL_fCa_gate (dimensionless).
  * RATES[9] is d/dt dT in component i_CaT_dT_gate (dimensionless).
  * RATES[10] is d/dt fT in component i_CaT_fT_gate (dimensionless).
  * RATES[11] is d/dt R in component Ca_SR_release (dimensionless).
  * RATES[12] is d/dt O in component Ca_SR_release (dimensionless).
  * RATES[13] is d/dt I in component Ca_SR_release (dimensionless).
  * RATES[14] is d/dt RI in component Ca_SR_release (dimensionless).
  * RATES[21] is d/dt fTC in component Ca_buffering (dimensionless).
  * RATES[22] is d/dt fTMC in component Ca_buffering (dimensionless).
  * RATES[18] is d/dt fTMM in component Ca_buffering (dimensionless).
  * RATES[19] is d/dt fCMi in component Ca_buffering (dimensionless).
  * RATES[20] is d/dt fCMs in component Ca_buffering (dimensionless).
  * RATES[23] is d/dt fCQ in component Ca_buffering (dimensionless).
  * RATES[17] is d/dt Cai in component Ca_dynamics (millimolar).
  * RATES[1] is d/dt Ca_sub in component Ca_dynamics (millimolar).
  * RATES[16] is d/dt Ca_nsr in component Ca_dynamics (millimolar).
  * RATES[15] is d/dt Ca_jsr in component Ca_dynamics (millimolar).
  * RATES[24] is d/dt r_Kur in component i_Kur_rKur_gate (dimensionless).
  * RATES[25] is d/dt s_Kur in component i_Kur_sKur_gate (dimensionless).
  * RATES[26] is d/dt q in component i_to_q_gate (dimensionless).
  * RATES[27] is d/dt r in component i_to_r_gate (dimensionless).
  * RATES[28] is d/dt paS in component i_Kr_pa_gate (dimensionless).
  * RATES[29] is d/dt paF in component i_Kr_pa_gate (dimensionless).
  * RATES[30] is d/dt piy in component i_Kr_pi_gate (dimensionless).
  * RATES[31] is d/dt n in component i_Ks_n_gate (dimensionless).
  * RATES[32] is d/dt a in component i_KACh_a_gate (dimensionless).
  */
  PDEFIELD_TYPE CONSTANTS_FS[116];
  PDEFIELD_TYPE ALGEBRAIC[101];

  CONSTANTS_FS[0] = 8314.472;
  CONSTANTS_FS[1] = 310;
  CONSTANTS_FS[2] = 96485.3415;
  CONSTANTS_FS[3] = 5.7e-5;
  CONSTANTS_FS[4] = 0;
  CONSTANTS_FS[5] = 0.5;
  CONSTANTS_FS[6] = 0.5;
  CONSTANTS_FS[7] = -35;
  CONSTANTS_FS[8] = -45;
  CONSTANTS_FS[9] = 0;
  CONSTANTS_FS[10] = 0;
  CONSTANTS_FS[11] = 140;
  CONSTANTS_FS[12] = 140;
  CONSTANTS_FS[13] = 5.4;
  CONSTANTS_FS[14] = 1.8;
  CONSTANTS_FS[15] = 1;
  CONSTANTS_FS[16] = 0.00427;
  CONSTANTS_FS[17] = 45;
  CONSTANTS_FS[18] = 0.5927;
  CONSTANTS_FS[19] = 0;
  CONSTANTS_FS[20] = 0;
  CONSTANTS_FS[21] = 1.4;
  CONSTANTS_FS[22] = 14;
  CONSTANTS_FS[23] = 0.08105;
  CONSTANTS_FS[24] = 3.343;
  CONSTANTS_FS[25] = 0.1369;
  CONSTANTS_FS[26] = 0.4315;
  CONSTANTS_FS[27] = 0;
  CONSTANTS_FS[28] = 26.44;
  CONSTANTS_FS[29] = 0.0207;
  CONSTANTS_FS[30] = 395.3;
  CONSTANTS_FS[31] = 2.289;
  CONSTANTS_FS[32] = 26.44;
  CONSTANTS_FS[33] = 4.663;
  CONSTANTS_FS[34] = 1628;
  CONSTANTS_FS[35] = 561.4;
  CONSTANTS_FS[36] = 3.663;
  CONSTANTS_FS[37] = 0;
  CONSTANTS_FS[38] = 0.0223;
  CONSTANTS_FS[39] = 0;
  CONSTANTS_FS[40] = 1e-5;
  CONSTANTS_FS[41] = 0.4578;
  CONSTANTS_FS[42] = 4.3371;
  CONSTANTS_FS[43] = -16.4508;
  CONSTANTS_FS[44] = 0;
  CONSTANTS_FS[45] = 0;
  CONSTANTS_FS[46] = 0.0075;
  CONSTANTS_FS[47] = 0.000338;
  CONSTANTS_FS[48] = 0.04132;
  CONSTANTS_FS[49] = 0;
  CONSTANTS_FS[50] = 148041085.1;
  CONSTANTS_FS[51] = 15;
  CONSTANTS_FS[52] = 1;
  CONSTANTS_FS[53] = 0.45;
  CONSTANTS_FS[54] = 2.5;
  CONSTANTS_FS[55] = 10000;
  CONSTANTS_FS[56] = 500;
  CONSTANTS_FS[57] = 5;
  CONSTANTS_FS[58] = 660;
  CONSTANTS_FS[59] = 5.469e-5;
  CONSTANTS_FS[60] = 0.04;
  CONSTANTS_FS[61] = 5;
  CONSTANTS_FS[62] = 0.000286113;
  CONSTANTS_FS[63] = 5e-5;
  CONSTANTS_FS[64] = 0.031;
  CONSTANTS_FS[65] = 0.062;
  CONSTANTS_FS[66] = 0.045;
  CONSTANTS_FS[67] = 10;
  CONSTANTS_FS[68] = 88800;
  CONSTANTS_FS[69] = 2277;
  CONSTANTS_FS[70] = 227700;
  CONSTANTS_FS[71] = 1.642e6;
  CONSTANTS_FS[72] = 175.4;
  CONSTANTS_FS[73] = 446;
  CONSTANTS_FS[74] = 7.51;
  CONSTANTS_FS[75] = 751;
  CONSTANTS_FS[76] = 542;
  CONSTANTS_FS[77] = 445;
  CONSTANTS_FS[78] = 2.5;
  CONSTANTS_FS[79] = 0.0012;
  CONSTANTS_FS[80] = 0.46;
  CONSTANTS_FS[81] = 0.0116;
  CONSTANTS_FS[82] = 3.9;
  CONSTANTS_FS[83] = 67;
  CONSTANTS_FS[84] = 0.02;
  CONSTANTS_FS[85] = 0.1539e-3;
  CONSTANTS_FS[86] = 3.5e-3;
  CONSTANTS_FS[87] = 0.00424;
  CONSTANTS_FS[88] = 0.00065;
  CONSTANTS_FS[89] = 0.00345;
  CONSTANTS_FS[90] = 1;
  CONSTANTS_FS[91] = ( CONSTANTS_FS[0]*CONSTANTS_FS[1])/CONSTANTS_FS[2];
  CONSTANTS_FS[92] = CONSTANTS_FS[16]/(CONSTANTS_FS[13]/(CONSTANTS_FS[13]+CONSTANTS_FS[17]));
  CONSTANTS_FS[93] = (CONSTANTS_FS[10]>0.00000 ? - 0.250000 : CONSTANTS_FS[9]>0.00000 ? ( 0.700000*CONSTANTS_FS[9])/(9.00000e-05+CONSTANTS_FS[9]) : 0.00000);
  CONSTANTS_FS[94] = (CONSTANTS_FS[10]>0.00000 ?  1.20000*CONSTANTS_FS[88] : CONSTANTS_FS[88]);
  CONSTANTS_FS[95] = (CONSTANTS_FS[9]>0.00000 ? - 1.00000 - ( 9.89800*pow( 1.00000*CONSTANTS_FS[9], 0.618000))/(pow( 1.00000*CONSTANTS_FS[9], 0.618000)+0.00122423) : 0.00000);
  CONSTANTS_FS[96] =  CONSTANTS_FS[91]*log(CONSTANTS_FS[13]/CONSTANTS_FS[12]);
  CONSTANTS_FS[97] = CONSTANTS_FS[92]/(CONSTANTS_FS[18]+1.00000);
  CONSTANTS_FS[98] =  CONSTANTS_FS[61]*(1.00000 - CONSTANTS_FS[93]);
  CONSTANTS_FS[99] = (CONSTANTS_FS[10]>0.00000 ? 7.50000 : 0.00000);
  CONSTANTS_FS[100] =  CONSTANTS_FS[18]*CONSTANTS_FS[97];
  CONSTANTS_FS[101] = ( CONSTANTS_FS[97]*CONSTANTS_FS[13])/(CONSTANTS_FS[13]+CONSTANTS_FS[17]);
  CONSTANTS_FS[102] = (CONSTANTS_FS[10]>0.00000 ? 1.20000 : 1.00000);
  CONSTANTS_FS[103] = ( CONSTANTS_FS[100]*CONSTANTS_FS[13])/(CONSTANTS_FS[13]+CONSTANTS_FS[17]);
  CONSTANTS_FS[104] = CONSTANTS_FS[11]/(CONSTANTS_FS[33]+CONSTANTS_FS[11]);
  CONSTANTS_FS[105] = (CONSTANTS_FS[10]>0.00000 ? 1.23000 : 1.00000);
  CONSTANTS_FS[106] = ( 0.310000*CONSTANTS_FS[9])/(CONSTANTS_FS[9]+9.00000e-05);
  CONSTANTS_FS[107] = (CONSTANTS_FS[10]>0.00000 ? - 8.00000 : 0.00000);
  CONSTANTS_FS[108] = (CONSTANTS_FS[10]>0.00000 ? - 27.0000 : 0.00000);
  CONSTANTS_FS[109] =  1.00000e-09* 3.14159265358979*pow(CONSTANTS_FS[82], 2.00000)*CONSTANTS_FS[83];
  CONSTANTS_FS[110] =  1.00000e-09*2.00000* 3.14159265358979*CONSTANTS_FS[84]*(CONSTANTS_FS[82] - CONSTANTS_FS[84]/2.00000)*CONSTANTS_FS[83];
  CONSTANTS_FS[111] =  CONSTANTS_FS[79]*CONSTANTS_FS[109];
  CONSTANTS_FS[112] =  CONSTANTS_FS[80]*CONSTANTS_FS[109] - CONSTANTS_FS[110];
  CONSTANTS_FS[113] =  CONSTANTS_FS[81]*CONSTANTS_FS[109];
  CONSTANTS_FS[114] = (CONSTANTS_FS[10]>0.00000 ? - 14.0000 : 0.00000);
  CONSTANTS_FS[115] = (3.59880 - 0.0256410)/(1.00000+1.21550e-06/pow( 1.00000*CONSTANTS_FS[9], 1.69510))+0.0256410;


  ALGEBRAIC[6] =  CONSTANTS_FS[69]*CONSTANTS_FS[78]*(1.00000 - (STATES[22]+STATES[18])) -  CONSTANTS_FS[75]*STATES[18];
  RATES[18] = ALGEBRAIC[6];
  ALGEBRAIC[1] = CONSTANTS_FS[47]/(CONSTANTS_FS[47]+STATES[1]);
  ALGEBRAIC[7] = ( 0.00100000*ALGEBRAIC[1])/CONSTANTS_FS[46];
  RATES[8] = (ALGEBRAIC[1] - STATES[8])/ALGEBRAIC[7];
  ALGEBRAIC[2] = CONSTANTS_FS[51] - (CONSTANTS_FS[51] - CONSTANTS_FS[52])/(1.00000+pow(CONSTANTS_FS[53]/STATES[15], CONSTANTS_FS[54]));
  ALGEBRAIC[8] = CONSTANTS_FS[55]/ALGEBRAIC[2];
  ALGEBRAIC[17] =  CONSTANTS_FS[56]*ALGEBRAIC[2];
  RATES[11] = ( CONSTANTS_FS[57]*STATES[14] -  ALGEBRAIC[17]*STATES[1]*STATES[11]) - ( ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[11] -  CONSTANTS_FS[58]*STATES[12]);
  RATES[12] = ( ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[11] -  CONSTANTS_FS[58]*STATES[12]) - ( ALGEBRAIC[17]*STATES[1]*STATES[12] -  CONSTANTS_FS[57]*STATES[13]);
  RATES[13] = ( ALGEBRAIC[17]*STATES[1]*STATES[12] -  CONSTANTS_FS[57]*STATES[13]) - ( CONSTANTS_FS[58]*STATES[13] -  ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[14]);
  RATES[14] = ( CONSTANTS_FS[58]*STATES[13] -  ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[14]) - ( CONSTANTS_FS[57]*STATES[14] -  ALGEBRAIC[17]*STATES[1]*STATES[11]);
  ALGEBRAIC[5] = (VOI>CONSTANTS_FS[5]&&VOI<CONSTANTS_FS[5]+CONSTANTS_FS[6] ? CONSTANTS_FS[7] : CONSTANTS_FS[8]);
  ALGEBRAIC[9] = (CONSTANTS_FS[4]>=1.00000 ? ALGEBRAIC[5] : STATES[0]);
  ALGEBRAIC[10] = 1.00000/(( 0.360000*(((ALGEBRAIC[9]+148.800) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]))/(exp( 0.0660000*(((ALGEBRAIC[9]+148.800) - CONSTANTS_FS[95]) - CONSTANTS_FS[99])) - 1.00000)+( 0.100000*(((ALGEBRAIC[9]+87.3000) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]))/(1.00000 - exp( - 0.200000*(((ALGEBRAIC[9]+87.3000) - CONSTANTS_FS[95]) - CONSTANTS_FS[99])))) - 0.0540000;
  ALGEBRAIC[29] = (ALGEBRAIC[9]<- (((80.0000 - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20]) ? 0.0132900+0.999210/(1.00000+exp(((((ALGEBRAIC[9]+97.1340) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20])/8.17520)) :  0.000250100*exp(- (((ALGEBRAIC[9] - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20])/12.8610));
  RATES[3] = (ALGEBRAIC[29] - STATES[3])/ALGEBRAIC[10];
  ALGEBRAIC[14] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+37.4000+CONSTANTS_FS[44])/(5.30000+CONSTANTS_FS[45])));
  ALGEBRAIC[33] =  0.00100000*(44.3000+ 230.000*exp(- pow((ALGEBRAIC[9]+36.0000)/10.0000, 2.00000)));
  RATES[7] = (ALGEBRAIC[14] - STATES[7])/ALGEBRAIC[33];
  ALGEBRAIC[15] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+38.3000)/5.50000));
  ALGEBRAIC[34] = 0.00100000/( 1.06800*exp((ALGEBRAIC[9]+38.3000)/30.0000)+ 1.06800*exp(- (ALGEBRAIC[9]+38.3000)/30.0000));
  RATES[9] = (ALGEBRAIC[15] - STATES[9])/ALGEBRAIC[34];
  ALGEBRAIC[16] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+58.7000)/3.80000));
  ALGEBRAIC[35] = 1.00000/( 16.6700*exp(- (ALGEBRAIC[9]+75.0000)/83.3000)+ 16.6700*exp((ALGEBRAIC[9]+75.0000)/15.3800))+CONSTANTS_FS[49];
  RATES[10] = (ALGEBRAIC[16] - STATES[10])/ALGEBRAIC[35];
  ALGEBRAIC[37] = 0.00900000/(1.00000+exp((ALGEBRAIC[9]+5.00000)/12.0000))+0.000500000;
  ALGEBRAIC[19] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+6.00000)/- 8.60000));
  RATES[24] = (ALGEBRAIC[19] - STATES[24])/ALGEBRAIC[37];
  ALGEBRAIC[38] = 0.590000/(1.00000+exp((ALGEBRAIC[9]+60.0000)/10.0000))+3.05000;
  ALGEBRAIC[20] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+7.50000)/10.0000));
  RATES[25] = (ALGEBRAIC[20] - STATES[25])/ALGEBRAIC[38];
  ALGEBRAIC[21] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+49.0000)/13.0000));
  ALGEBRAIC[39] =  0.00100000*0.600000*(65.1700/( 0.570000*exp( - 0.0800000*(ALGEBRAIC[9]+44.0000))+ 0.0650000*exp( 0.100000*(ALGEBRAIC[9]+45.9300)))+10.1000);
  RATES[26] = (ALGEBRAIC[21] - STATES[26])/ALGEBRAIC[39];
  ALGEBRAIC[22] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9] - 19.3000)/15.0000));
  ALGEBRAIC[40] =  0.00100000*0.660000*1.40000*(15.5900/( 1.03700*exp( 0.0900000*(ALGEBRAIC[9]+30.6100))+ 0.369000*exp( - 0.120000*(ALGEBRAIC[9]+23.8400)))+2.98000);
  RATES[27] = (ALGEBRAIC[22] - STATES[27])/ALGEBRAIC[40];
  ALGEBRAIC[23] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+10.0144)/7.66070));
  ALGEBRAIC[41] = 0.846554/( 4.20000*exp(ALGEBRAIC[9]/17.0000)+ 0.150000*exp(- ALGEBRAIC[9]/21.6000));
  RATES[28] = (ALGEBRAIC[23] - STATES[28])/ALGEBRAIC[41];
  ALGEBRAIC[42] = 1.00000/( 30.0000*exp(ALGEBRAIC[9]/10.0000)+exp(- ALGEBRAIC[9]/12.0000));
  RATES[29] = (ALGEBRAIC[23] - STATES[29])/ALGEBRAIC[42];
  ALGEBRAIC[43] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+28.6000)/17.1000));
  ALGEBRAIC[26] = 1.00000/( 100.000*exp(- ALGEBRAIC[9]/54.6450)+ 656.000*exp(ALGEBRAIC[9]/106.157));
  RATES[30] = (ALGEBRAIC[43] - STATES[30])/ALGEBRAIC[26];
  ALGEBRAIC[28] =  10.0000*exp( 0.0133000*(ALGEBRAIC[9]+40.0000));
  ALGEBRAIC[45] = CONSTANTS_FS[115]/(CONSTANTS_FS[115]+ALGEBRAIC[28]);
  ALGEBRAIC[51] = 1.00000/(CONSTANTS_FS[115]+ALGEBRAIC[28]);
  RATES[32] = (ALGEBRAIC[45] - STATES[32])/ALGEBRAIC[51];
  ALGEBRAIC[12] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+69.8040)/4.45650));
  ALGEBRAIC[31] =  20.0000*exp( - 0.125000*(ALGEBRAIC[9]+75.0000));
  ALGEBRAIC[47] = 2000.00/( 320.000*exp( - 0.100000*(ALGEBRAIC[9]+75.0000))+1.00000);
  ALGEBRAIC[53] = 1.00000/(ALGEBRAIC[31]+ALGEBRAIC[47]);
  RATES[5] = (ALGEBRAIC[12] - STATES[5])/ALGEBRAIC[53];
  ALGEBRAIC[27] =  pow((1.00000/(1.00000+exp(- ((ALGEBRAIC[9]+0.638300) - CONSTANTS_FS[114])/10.7071))), 1.0 / 2);
  ALGEBRAIC[44] = 28.0000/(1.00000+exp(- ((ALGEBRAIC[9] - 40.0000) - CONSTANTS_FS[114])/3.00000));
  ALGEBRAIC[50] =  1.00000*exp(- ((ALGEBRAIC[9] - CONSTANTS_FS[114]) - 5.00000)/25.0000);
  ALGEBRAIC[56] = 1.00000/(ALGEBRAIC[44]+ALGEBRAIC[50]);
  RATES[31] = (ALGEBRAIC[27] - STATES[31])/ALGEBRAIC[56];
  ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+42.0504)/8.31060));
  ALGEBRAIC[30] = ALGEBRAIC[9]+41.0000;
  ALGEBRAIC[46] = (fabs(ALGEBRAIC[30])<CONSTANTS_FS[40] ? 2000.00 : ( 200.000*ALGEBRAIC[30])/(1.00000 - exp( - 0.100000*ALGEBRAIC[30])));
  ALGEBRAIC[52] =  8000.00*exp( - 0.0560000*(ALGEBRAIC[9]+66.0000));
  ALGEBRAIC[57] = 1.00000/(ALGEBRAIC[46]+ALGEBRAIC[52]);
  RATES[4] = (ALGEBRAIC[11] - STATES[4])/ALGEBRAIC[57];
  ALGEBRAIC[13] = 1.00000/(1.00000+exp(- ((ALGEBRAIC[9] - CONSTANTS_FS[43]) - CONSTANTS_FS[107])/( CONSTANTS_FS[42]*(1.00000+CONSTANTS_FS[108]/100.000))));
  ALGEBRAIC[32] = (ALGEBRAIC[9]==- 41.8000 ? - 41.8000 : ALGEBRAIC[9]==0.00000 ? 0.00000 : ALGEBRAIC[9]==- 6.80000 ? - 6.80001 : ALGEBRAIC[9]);
  ALGEBRAIC[48] = ( - 0.0283900*(ALGEBRAIC[32]+41.8000))/(exp(- (ALGEBRAIC[32]+41.8000)/2.50000) - 1.00000) - ( 0.0849000*(ALGEBRAIC[32]+6.80000))/(exp(- (ALGEBRAIC[32]+6.80000)/4.80000) - 1.00000);
  ALGEBRAIC[54] = (ALGEBRAIC[9]==- 1.80000 ? - 1.80001 : ALGEBRAIC[9]);
  ALGEBRAIC[58] = ( 0.0114300*(ALGEBRAIC[54]+1.80000))/(exp((ALGEBRAIC[54]+1.80000)/2.50000) - 1.00000);
  ALGEBRAIC[60] = 0.00100000/(ALGEBRAIC[48]+ALGEBRAIC[58]);
  RATES[6] = (ALGEBRAIC[13] - STATES[6])/ALGEBRAIC[60];
  ALGEBRAIC[18] = STATES[2];
  ALGEBRAIC[36] =  CONSTANTS_FS[91]*log(CONSTANTS_FS[11]/ALGEBRAIC[18]);
  ALGEBRAIC[61] =  CONSTANTS_FS[102]*CONSTANTS_FS[23]*pow(1.00000+pow(CONSTANTS_FS[21]/CONSTANTS_FS[13], 1.20000), - 1.00000)*pow(1.00000+pow(CONSTANTS_FS[22]/ALGEBRAIC[18], 1.30000), - 1.00000)*pow(1.00000+exp(- ((ALGEBRAIC[9] - ALGEBRAIC[36])+110.000)/20.0000), - 1.00000);
  ALGEBRAIC[63] = exp(( - CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91]));
  ALGEBRAIC[69] = 1.00000+ (CONSTANTS_FS[14]/CONSTANTS_FS[36])*(1.00000+exp(( CONSTANTS_FS[27]*ALGEBRAIC[9])/CONSTANTS_FS[91]))+ (CONSTANTS_FS[11]/CONSTANTS_FS[34])*(1.00000+ (CONSTANTS_FS[11]/CONSTANTS_FS[35])*(1.00000+CONSTANTS_FS[11]/CONSTANTS_FS[33]));
  ALGEBRAIC[71] = ( (( (CONSTANTS_FS[11]/CONSTANTS_FS[34])*CONSTANTS_FS[11])/CONSTANTS_FS[35])*(1.00000+CONSTANTS_FS[11]/CONSTANTS_FS[33])*exp(( - CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91])))/ALGEBRAIC[69];
  ALGEBRAIC[70] = ( (CONSTANTS_FS[14]/CONSTANTS_FS[36])*exp(( CONSTANTS_FS[27]*ALGEBRAIC[9])/CONSTANTS_FS[91]))/ALGEBRAIC[69];
  ALGEBRAIC[67] = exp(( CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91]));
  ALGEBRAIC[62] = ALGEBRAIC[18]/(CONSTANTS_FS[28]+ALGEBRAIC[18]);
  ALGEBRAIC[72] =  ALGEBRAIC[63]*CONSTANTS_FS[104]*(ALGEBRAIC[71]+ALGEBRAIC[70])+ ALGEBRAIC[70]*ALGEBRAIC[67]*(ALGEBRAIC[62]+ALGEBRAIC[63]);
  ALGEBRAIC[64] = 1.00000+ (STATES[1]/CONSTANTS_FS[29])*(1.00000+exp(( - CONSTANTS_FS[25]*ALGEBRAIC[9])/CONSTANTS_FS[91])+ALGEBRAIC[18]/CONSTANTS_FS[32])+ (ALGEBRAIC[18]/CONSTANTS_FS[30])*(1.00000+ (ALGEBRAIC[18]/CONSTANTS_FS[31])*(1.00000+ALGEBRAIC[18]/CONSTANTS_FS[28]));
  ALGEBRAIC[65] = ( (STATES[1]/CONSTANTS_FS[29])*exp(( - CONSTANTS_FS[25]*ALGEBRAIC[9])/CONSTANTS_FS[91]))/ALGEBRAIC[64];
  ALGEBRAIC[66] = ( (( (ALGEBRAIC[18]/CONSTANTS_FS[30])*ALGEBRAIC[18])/CONSTANTS_FS[31])*(1.00000+ALGEBRAIC[18]/CONSTANTS_FS[28])*exp(( CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91])))/ALGEBRAIC[64];
  ALGEBRAIC[68] =  ALGEBRAIC[67]*ALGEBRAIC[62]*(ALGEBRAIC[66]+ALGEBRAIC[65])+ ALGEBRAIC[63]*ALGEBRAIC[65]*(CONSTANTS_FS[104]+ALGEBRAIC[67]);
  ALGEBRAIC[73] =  ALGEBRAIC[66]*ALGEBRAIC[62]*(ALGEBRAIC[71]+ALGEBRAIC[70])+ ALGEBRAIC[65]*ALGEBRAIC[71]*(ALGEBRAIC[62]+ALGEBRAIC[63]);
  ALGEBRAIC[74] =  ALGEBRAIC[71]*CONSTANTS_FS[104]*(ALGEBRAIC[66]+ALGEBRAIC[65])+ ALGEBRAIC[66]*ALGEBRAIC[70]*(CONSTANTS_FS[104]+ALGEBRAIC[67]);
  ALGEBRAIC[75] = ( (1.00000 - CONSTANTS_FS[37])*CONSTANTS_FS[24]*( ALGEBRAIC[68]*ALGEBRAIC[70] -  ALGEBRAIC[72]*ALGEBRAIC[65]))/(ALGEBRAIC[72]+ALGEBRAIC[68]+ALGEBRAIC[73]+ALGEBRAIC[74]);
  ALGEBRAIC[76] =  CONSTANTS_FS[91]*log((CONSTANTS_FS[11]+ 0.120000*CONSTANTS_FS[13])/(ALGEBRAIC[18]+ 0.120000*CONSTANTS_FS[12]));
  ALGEBRAIC[77] =  CONSTANTS_FS[38]*pow(STATES[4], 3.00000)*STATES[5]*(ALGEBRAIC[9] - ALGEBRAIC[76]);
  ALGEBRAIC[78] =  CONSTANTS_FS[39]*pow(STATES[4], 3.00000)*(ALGEBRAIC[9] - ALGEBRAIC[76]);
  ALGEBRAIC[79] = ALGEBRAIC[77]+ALGEBRAIC[78];
  ALGEBRAIC[49] =  STATES[3]*CONSTANTS_FS[103]*(ALGEBRAIC[9] - ALGEBRAIC[36])*(1.00000 - CONSTANTS_FS[19]);
  ALGEBRAIC[82] =  (( 1.85000e-05*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))))*(ALGEBRAIC[18] -  CONSTANTS_FS[11]*exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  RATES[2] = ( (1.00000 - CONSTANTS_FS[15])*- 1.00000*(ALGEBRAIC[79]+ALGEBRAIC[49]+ALGEBRAIC[82]+ 3.00000*ALGEBRAIC[61]+ 3.00000*ALGEBRAIC[75]))/( 1.00000*(CONSTANTS_FS[112]+CONSTANTS_FS[110])*CONSTANTS_FS[2]);
  ALGEBRAIC[90] =  CONSTANTS_FS[71]*STATES[1]*(1.00000 - STATES[20]) -  CONSTANTS_FS[76]*STATES[20];
  RATES[20] = ALGEBRAIC[90];
  ALGEBRAIC[84] =  (( 2.00000*CONSTANTS_FS[48]*ALGEBRAIC[9])/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*ALGEBRAIC[9]*2.00000)/CONSTANTS_FS[91]))))*(STATES[1] -  CONSTANTS_FS[14]*exp(( - 2.00000*ALGEBRAIC[9])/CONSTANTS_FS[91]))*STATES[9]*STATES[10];
  ALGEBRAIC[80] =  (( 2.00000*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000)*2.00000)/CONSTANTS_FS[91]))))*(STATES[1] -  CONSTANTS_FS[14]*exp(( - 2.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  ALGEBRAIC[86] =  CONSTANTS_FS[50]*STATES[12]*(STATES[15] - STATES[1]);
  ALGEBRAIC[88] = (STATES[1] - STATES[17])/CONSTANTS_FS[59];
  RATES[1] = ( ALGEBRAIC[86]*CONSTANTS_FS[111])/CONSTANTS_FS[110] - (((ALGEBRAIC[80]+ALGEBRAIC[84]) -  2.00000*ALGEBRAIC[75])/( 2.00000*CONSTANTS_FS[2]*CONSTANTS_FS[110])+ALGEBRAIC[88]+ CONSTANTS_FS[66]*ALGEBRAIC[90]);
  ALGEBRAIC[93] =  CONSTANTS_FS[68]*STATES[17]*(1.00000 - STATES[21]) -  CONSTANTS_FS[73]*STATES[21];
  RATES[21] = ALGEBRAIC[93];
  ALGEBRAIC[91] = CONSTANTS_FS[98]/(1.00000+exp((- STATES[17]+CONSTANTS_FS[62])/CONSTANTS_FS[63]));
  ALGEBRAIC[94] = (STATES[16] - STATES[15])/CONSTANTS_FS[60];
  RATES[16] = ALGEBRAIC[91] - ( ALGEBRAIC[94]*CONSTANTS_FS[111])/CONSTANTS_FS[113];
  ALGEBRAIC[96] =  CONSTANTS_FS[70]*STATES[17]*(1.00000 - (STATES[22]+STATES[18])) -  CONSTANTS_FS[74]*STATES[22];
  RATES[22] = ALGEBRAIC[96];
  ALGEBRAIC[97] =  CONSTANTS_FS[72]*STATES[15]*(1.00000 - STATES[23]) -  CONSTANTS_FS[77]*STATES[23];
  RATES[23] = ALGEBRAIC[97];
  RATES[15] = ALGEBRAIC[94] - (ALGEBRAIC[86]+ CONSTANTS_FS[67]*ALGEBRAIC[97]);
  ALGEBRAIC[99] =  CONSTANTS_FS[71]*STATES[17]*(1.00000 - STATES[19]) -  CONSTANTS_FS[76]*STATES[19];
  RATES[19] = ALGEBRAIC[99];
  RATES[17] = ( 1.00000*( ALGEBRAIC[88]*CONSTANTS_FS[110] -  ALGEBRAIC[91]*CONSTANTS_FS[113]))/CONSTANTS_FS[112] - ( CONSTANTS_FS[66]*ALGEBRAIC[99]+ CONSTANTS_FS[64]*ALGEBRAIC[93]+ CONSTANTS_FS[65]*ALGEBRAIC[96]);
  ALGEBRAIC[55] =  STATES[3]*CONSTANTS_FS[101]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*(1.00000 - CONSTANTS_FS[19]);
  ALGEBRAIC[59] = ALGEBRAIC[49]+ALGEBRAIC[55];
  ALGEBRAIC[89] =  CONSTANTS_FS[87]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*( 0.900000*STATES[29]+ 0.100000*STATES[28])*STATES[30];
  ALGEBRAIC[92] =  CONSTANTS_FS[91]*log((CONSTANTS_FS[13]+ 0.120000*CONSTANTS_FS[11])/(CONSTANTS_FS[12]+ 0.120000*ALGEBRAIC[18]));
  ALGEBRAIC[95] =  CONSTANTS_FS[94]*(ALGEBRAIC[9] - ALGEBRAIC[92])*pow(STATES[31], 2.00000);
  ALGEBRAIC[87] =  CONSTANTS_FS[86]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*STATES[26]*STATES[27];
  ALGEBRAIC[81] =  (( 0.000365000*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))))*(CONSTANTS_FS[12] -  CONSTANTS_FS[13]*exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  ALGEBRAIC[83] =  (ALGEBRAIC[80]+ALGEBRAIC[81]+ALGEBRAIC[82])*(1.00000 - CONSTANTS_FS[106])*1.00000*CONSTANTS_FS[105];
  ALGEBRAIC[98] = (CONSTANTS_FS[9]>0.00000 ?  CONSTANTS_FS[90]*CONSTANTS_FS[89]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*(1.00000+exp((ALGEBRAIC[9]+20.0000)/20.0000))*STATES[32] : 0.00000);
  ALGEBRAIC[85] =  CONSTANTS_FS[85]*STATES[24]*STATES[25]*(ALGEBRAIC[9] - CONSTANTS_FS[96]);
  ALGEBRAIC[100] = I_f_factor*ALGEBRAIC[59]+I_Kr_factor*ALGEBRAIC[89]+ALGEBRAIC[95]+ALGEBRAIC[87]+ALGEBRAIC[61]+ALGEBRAIC[75]+ALGEBRAIC[79]+ALGEBRAIC[83]+ALGEBRAIC[84]+ALGEBRAIC[98]+ALGEBRAIC[85];
  RATES[0] = - ALGEBRAIC[100]/CONSTANTS_FS[3];

}

__device__ void RushLarsenFabbriSeveri(PDEFIELD_TYPE VOI, PDEFIELD_TYPE* STATES, PDEFIELD_TYPE* STATES_NEW, PDEFIELD_TYPE* RATES, PDEFIELD_TYPE dt, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE I_f_factor, PDEFIELD_TYPE I_Kr_factor, int id){
  /*
   There are a total of 101 entries in the algebraic variable array.
   There are a total of 33 entries in each of the rate and state variable arrays.
   There are a total of 116 entries in the constant variable array.
 */
  /*
  * VOI is time in component environment (second).
  * CONSTANTS_FS[0] is R in component Membrane (joule_per_kilomole_kelvin).
  * CONSTANTS_FS[1] is T in component Membrane (kelvin).
  * CONSTANTS_FS[2] is F in component Membrane (coulomb_per_mole).
  * CONSTANTS_FS[3] is C in component Membrane (microF).
  * CONSTANTS_FS[91] is RTONF in component Membrane (millivolt).
  * ALGEBRAIC[59] is i_f in component i_f (nanoA).
  * ALGEBRAIC[61] is i_NaK in component i_NaK (nanoA).
  * ALGEBRAIC[75] is i_NaCa in component i_NaCa (nanoA).
  * ALGEBRAIC[79] is i_Na in component i_Na (nanoA).
  * ALGEBRAIC[89] is i_Kr in component i_Kr (nanoA).
  * ALGEBRAIC[95] is i_Ks in component i_Ks (nanoA).
  * ALGEBRAIC[87] is i_to in component i_to (nanoA).
  * ALGEBRAIC[83] is i_CaL in component i_CaL (nanoA).
  * ALGEBRAIC[84] is i_CaT in component i_CaT (nanoA).
  * ALGEBRAIC[98] is i_KACh in component i_KACh (nanoA).
  * ALGEBRAIC[85] is i_Kur in component i_Kur (nanoA).
  * ALGEBRAIC[9] is V in component Membrane (millivolt).
  * CONSTANTS_FS[4] is clamp_mode in component Membrane (dimensionless).
  * ALGEBRAIC[5] is V_clamp in component Voltage_clamp (millivolt).
  * STATES[0] is V_ode in component Membrane (millivolt).
  * ALGEBRAIC[100] is i_tot in component Membrane (nanoA).
  * CONSTANTS_FS[5] is t_holding in component Voltage_clamp (second).
  * CONSTANTS_FS[6] is t_test in component Voltage_clamp (second).
  * CONSTANTS_FS[7] is V_test in component Voltage_clamp (millivolt).
  * CONSTANTS_FS[8] is V_holding in component Voltage_clamp (millivolt).
  * CONSTANTS_FS[9] is ACh in component Rate_modulation_experiments (millimolar).
  * CONSTANTS_FS[10] is Iso_1_uM in component Rate_modulation_experiments (dimensionless).
  * ALGEBRAIC[18] is Nai in component Nai_concentration (millimolar).
  * CONSTANTS_FS[11] is Nao in component Ionic_values (millimolar).
  * CONSTANTS_FS[12] is Ki in component Ionic_values (millimolar).
  * CONSTANTS_FS[13] is Ko in component Ionic_values (millimolar).
  * STATES[1] is Ca_sub in component Ca_dynamics (millimolar).
  * CONSTANTS_FS[14] is Cao in component Ionic_values (millimolar).
  * ALGEBRAIC[36] is E_Na in component Ionic_values (millivolt).
  * CONSTANTS_FS[96] is E_K in component Ionic_values (millivolt).
  * ALGEBRAIC[0] is E_Ca in component Ionic_values (millivolt).
  * CONSTANTS_FS[110] is V_sub in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[112] is V_i in component Cell_parameters (millimetre3).
  * ALGEBRAIC[49] is i_fNa in component i_f (nanoA).
  * ALGEBRAIC[82] is i_siNa in component i_CaL (nanoA).
  * STATES[2] is Nai_ in component Nai_concentration (millimolar).
  * CONSTANTS_FS[15] is Nai_clamp in component Nai_concentration (dimensionless).
  * ALGEBRAIC[55] is i_fK in component i_f (nanoA).
  * CONSTANTS_FS[16] is g_f in component i_f (microS).
  * CONSTANTS_FS[92] is G_f in component i_f (microS).
  * CONSTANTS_FS[103] is g_f_Na in component i_f (microS).
  * CONSTANTS_FS[100] is G_f_Na in component i_f (microS).
  * CONSTANTS_FS[101] is g_f_K in component i_f (microS).
  * CONSTANTS_FS[97] is G_f_K in component i_f (microS).
  * CONSTANTS_FS[17] is Km_f in component i_f (millimolar).
  * CONSTANTS_FS[18] is alpha in component i_f (dimensionless).
  * STATES[3] is y in component i_f_y_gate (dimensionless).
  * CONSTANTS_FS[19] is blockade in component i_f (dimensionless).
  * ALGEBRAIC[10] is tau_y in component i_f_y_gate (second).
  * ALGEBRAIC[29] is y_infinity in component i_f_y_gate (dimensionless).
  * CONSTANTS_FS[95] is ACh_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[99] is Iso_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[20] is y_shift in component i_f_y_gate (millivolt).
  * CONSTANTS_FS[21] is Km_Kp in component i_NaK (millimolar).
  * CONSTANTS_FS[22] is Km_Nap in component i_NaK (millimolar).
  * CONSTANTS_FS[23] is i_NaK_max in component i_NaK (nanoA).
  * CONSTANTS_FS[102] is Iso_increase in component i_NaK (dimensionless).
  * CONSTANTS_FS[24] is K_NaCa in component i_NaCa (nanoA).
  * ALGEBRAIC[72] is x1 in component i_NaCa (dimensionless).
  * ALGEBRAIC[68] is x2 in component i_NaCa (dimensionless).
  * ALGEBRAIC[73] is x3 in component i_NaCa (dimensionless).
  * ALGEBRAIC[74] is x4 in component i_NaCa (dimensionless).
  * ALGEBRAIC[63] is k41 in component i_NaCa (dimensionless).
  * CONSTANTS_FS[104] is k34 in component i_NaCa (dimensionless).
  * ALGEBRAIC[71] is k23 in component i_NaCa (dimensionless).
  * ALGEBRAIC[70] is k21 in component i_NaCa (dimensionless).
  * ALGEBRAIC[67] is k32 in component i_NaCa (dimensionless).
  * ALGEBRAIC[62] is k43 in component i_NaCa (dimensionless).
  * ALGEBRAIC[65] is k12 in component i_NaCa (dimensionless).
  * ALGEBRAIC[66] is k14 in component i_NaCa (dimensionless).
  * CONSTANTS_FS[25] is Qci in component i_NaCa (dimensionless).
  * CONSTANTS_FS[26] is Qn in component i_NaCa (dimensionless).
  * CONSTANTS_FS[27] is Qco in component i_NaCa (dimensionless).
  * CONSTANTS_FS[28] is K3ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[29] is Kci in component i_NaCa (millimolar).
  * CONSTANTS_FS[30] is K1ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[31] is K2ni in component i_NaCa (millimolar).
  * CONSTANTS_FS[32] is Kcni in component i_NaCa (millimolar).
  * CONSTANTS_FS[33] is K3no in component i_NaCa (millimolar).
  * CONSTANTS_FS[34] is K1no in component i_NaCa (millimolar).
  * CONSTANTS_FS[35] is K2no in component i_NaCa (millimolar).
  * CONSTANTS_FS[36] is Kco in component i_NaCa (millimolar).
  * ALGEBRAIC[69] is do in component i_NaCa (dimensionless).
  * ALGEBRAIC[64] is di in component i_NaCa (dimensionless).
  * CONSTANTS_FS[37] is blockade_NaCa in component i_NaCa (dimensionless).
  * ALGEBRAIC[77] is i_Na_ in component i_Na (nanoA).
  * ALGEBRAIC[78] is i_Na_L in component i_Na (nanoA).
  * CONSTANTS_FS[38] is g_Na in component i_Na (microS).
  * CONSTANTS_FS[39] is g_Na_L in component i_Na (microS).
  * ALGEBRAIC[76] is E_mh in component i_Na (millivolt).
  * STATES[4] is m in component i_Na_m_gate (dimensionless).
  * STATES[5] is h in component i_Na_h_gate (dimensionless).
  * ALGEBRAIC[46] is alpha_m in component i_Na_m_gate (per_second).
  * ALGEBRAIC[52] is beta_m in component i_Na_m_gate (per_second).
  * ALGEBRAIC[11] is m_infinity in component i_Na_m_gate (dimensionless).
  * ALGEBRAIC[57] is tau_m in component i_Na_m_gate (second).
  * CONSTANTS_FS[40] is delta_m in component i_Na_m_gate (millivolt).
  * ALGEBRAIC[30] is E0_m in component i_Na_m_gate (millivolt).
  * ALGEBRAIC[31] is alpha_h in component i_Na_h_gate (per_second).
  * ALGEBRAIC[47] is beta_h in component i_Na_h_gate (per_second).
  * ALGEBRAIC[12] is h_infinity in component i_Na_h_gate (dimensionless).
  * ALGEBRAIC[53] is tau_h in component i_Na_h_gate (second).
  * ALGEBRAIC[80] is i_siCa in component i_CaL (nanoA).
  * ALGEBRAIC[81] is i_siK in component i_CaL (nanoA).
  * CONSTANTS_FS[106] is ACh_block in component i_CaL (dimensionless).
  * CONSTANTS_FS[41] is P_CaL in component i_CaL (nanoA_per_millimolar).
  * STATES[6] is dL in component i_CaL_dL_gate (dimensionless).
  * STATES[7] is fL in component i_CaL_fL_gate (dimensionless).
  * STATES[8] is fCa in component i_CaL_fCa_gate (dimensionless).
  * CONSTANTS_FS[105] is Iso_increase in component i_CaL (dimensionless).
  * ALGEBRAIC[13] is dL_infinity in component i_CaL_dL_gate (dimensionless).
  * ALGEBRAIC[60] is tau_dL in component i_CaL_dL_gate (second).
  * ALGEBRAIC[48] is alpha_dL in component i_CaL_dL_gate (per_second).
  * ALGEBRAIC[58] is beta_dL in component i_CaL_dL_gate (per_second).
  * ALGEBRAIC[32] is adVm in component i_CaL_dL_gate (millivolt).
  * ALGEBRAIC[54] is bdVm in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[42] is k_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[43] is V_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[107] is Iso_shift_dL in component i_CaL_dL_gate (millivolt).
  * CONSTANTS_FS[108] is Iso_slope_dL in component i_CaL_dL_gate (dimensionless).
  * ALGEBRAIC[14] is fL_infinity in component i_CaL_fL_gate (dimensionless).
  * ALGEBRAIC[33] is tau_fL in component i_CaL_fL_gate (second).
  * CONSTANTS_FS[44] is shift_fL in component i_CaL_fL_gate (millivolt).
  * CONSTANTS_FS[45] is k_fL in component i_CaL_fL_gate (millivolt).
  * CONSTANTS_FS[46] is alpha_fCa in component i_CaL_fCa_gate (per_second).
  * ALGEBRAIC[1] is fCa_infinity in component i_CaL_fCa_gate (dimensionless).
  * ALGEBRAIC[7] is tau_fCa in component i_CaL_fCa_gate (second).
  * CONSTANTS_FS[47] is Km_fCa in component i_CaL_fCa_gate (millimolar).
  * CONSTANTS_FS[48] is P_CaT in component i_CaT (nanoA_per_millimolar).
  * STATES[9] is dT in component i_CaT_dT_gate (dimensionless).
  * STATES[10] is fT in component i_CaT_fT_gate (dimensionless).
  * ALGEBRAIC[15] is dT_infinity in component i_CaT_dT_gate (dimensionless).
  * ALGEBRAIC[34] is tau_dT in component i_CaT_dT_gate (second).
  * ALGEBRAIC[16] is fT_infinity in component i_CaT_fT_gate (dimensionless).
  * ALGEBRAIC[35] is tau_fT in component i_CaT_fT_gate (second).
  * CONSTANTS_FS[49] is offset_fT in component i_CaT_fT_gate (second).
  * ALGEBRAIC[86] is j_SRCarel in component Ca_SR_release (millimolar_per_second).
  * STATES[11] is R in component Ca_SR_release (dimensionless).
  * STATES[12] is O in component Ca_SR_release (dimensionless).
  * STATES[13] is I in component Ca_SR_release (dimensionless).
  * STATES[14] is RI in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[50] is ks in component Ca_SR_release (per_second).
  * CONSTANTS_FS[51] is MaxSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[52] is MinSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[53] is EC50_SR in component Ca_SR_release (millimolar).
  * CONSTANTS_FS[54] is HSR in component Ca_SR_release (dimensionless).
  * ALGEBRAIC[8] is koSRCa in component Ca_SR_release (per_millimolar2_second).
  * ALGEBRAIC[17] is kiSRCa in component Ca_SR_release (per_millimolar_second).
  * CONSTANTS_FS[55] is koCa in component Ca_SR_release (per_millimolar2_second).
  * CONSTANTS_FS[56] is kiCa in component Ca_SR_release (per_millimolar_second).
  * ALGEBRAIC[2] is kCaSR in component Ca_SR_release (dimensionless).
  * CONSTANTS_FS[57] is kim in component Ca_SR_release (per_second).
  * CONSTANTS_FS[58] is kom in component Ca_SR_release (per_second).
  * STATES[15] is Ca_jsr in component Ca_dynamics (millimolar).
  * ALGEBRAIC[3] is diff in component Ca_SR_release (millimolar).
  * ALGEBRAIC[4] is P_tot in component Ca_SR_release (dimensionless).
  * ALGEBRAIC[88] is j_Ca_dif in component Ca_intracellular_fluxes (millimolar_per_second).
  * ALGEBRAIC[91] is j_up in component Ca_intracellular_fluxes (millimolar_per_second).
  * ALGEBRAIC[94] is j_tr in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[59] is tau_dif_Ca in component Ca_intracellular_fluxes (second).
  * CONSTANTS_FS[60] is tau_tr in component Ca_intracellular_fluxes (second).
  * CONSTANTS_FS[98] is P_up in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[61] is P_up_basal in component Ca_intracellular_fluxes (millimolar_per_second).
  * CONSTANTS_FS[93] is b_up in component Ca_intracellular_fluxes (dimensionless).
  * CONSTANTS_FS[62] is K_up in component Ca_intracellular_fluxes (millimolar).
  * STATES[16] is Ca_nsr in component Ca_dynamics (millimolar).
  * STATES[17] is Cai in component Ca_dynamics (millimolar).
  * CONSTANTS_FS[63] is slope_up in component Ca_intracellular_fluxes (millimolar).
  * CONSTANTS_FS[64] is TC_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[65] is TMC_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[66] is CM_tot in component Ca_buffering (millimolar).
  * CONSTANTS_FS[67] is CQ_tot in component Ca_buffering (millimolar).
  * ALGEBRAIC[93] is delta_fTC in component Ca_buffering (per_second).
  * ALGEBRAIC[96] is delta_fTMC in component Ca_buffering (per_second).
  * ALGEBRAIC[90] is delta_fCMs in component Ca_buffering (per_second).
  * ALGEBRAIC[99] is delta_fCMi in component Ca_buffering (per_second).
  * ALGEBRAIC[97] is delta_fCQ in component Ca_buffering (per_second).
  * ALGEBRAIC[6] is delta_fTMM in component Ca_buffering (per_second).
  * STATES[18] is fTMM in component Ca_buffering (dimensionless).
  * STATES[19] is fCMi in component Ca_buffering (dimensionless).
  * STATES[20] is fCMs in component Ca_buffering (dimensionless).
  * STATES[21] is fTC in component Ca_buffering (dimensionless).
  * STATES[22] is fTMC in component Ca_buffering (dimensionless).
  * STATES[23] is fCQ in component Ca_buffering (dimensionless).
  * CONSTANTS_FS[68] is kf_TC in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[69] is kf_TMM in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[70] is kf_TMC in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[71] is kf_CM in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[72] is kf_CQ in component Ca_buffering (per_millimolar_second).
  * CONSTANTS_FS[73] is kb_TC in component Ca_buffering (per_second).
  * CONSTANTS_FS[74] is kb_TMC in component Ca_buffering (per_second).
  * CONSTANTS_FS[75] is kb_TMM in component Ca_buffering (per_second).
  * CONSTANTS_FS[76] is kb_CM in component Ca_buffering (per_second).
  * CONSTANTS_FS[77] is kb_CQ in component Ca_buffering (per_second).
  * CONSTANTS_FS[78] is Mgi in component Ca_buffering (millimolar).
  * CONSTANTS_FS[111] is V_jsr in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[113] is V_nsr in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[109] is V_cell in component Cell_parameters (millimetre3).
  * CONSTANTS_FS[79] is V_jsr_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[80] is V_i_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[81] is V_nsr_part in component Cell_parameters (dimensionless).
  * CONSTANTS_FS[82] is R_cell in component Cell_parameters (micrometre).
  * CONSTANTS_FS[83] is L_cell in component Cell_parameters (micrometre).
  * CONSTANTS_FS[84] is L_sub in component Cell_parameters (micrometre).
  * CONSTANTS_FS[85] is g_Kur in component i_Kur (microS).
  * STATES[24] is r_Kur in component i_Kur_rKur_gate (dimensionless).
  * STATES[25] is s_Kur in component i_Kur_sKur_gate (dimensionless).
  * ALGEBRAIC[37] is tau_r_Kur in component i_Kur_rKur_gate (second).
  * ALGEBRAIC[19] is r_Kur_infinity in component i_Kur_rKur_gate (dimensionless).
  * ALGEBRAIC[38] is tau_s_Kur in component i_Kur_sKur_gate (second).
  * ALGEBRAIC[20] is s_Kur_infinity in component i_Kur_sKur_gate (dimensionless).
  * CONSTANTS_FS[86] is g_to in component i_to (microS).
  * STATES[26] is q in component i_to_q_gate (dimensionless).
  * STATES[27] is r in component i_to_r_gate (dimensionless).
  * ALGEBRAIC[21] is q_infinity in component i_to_q_gate (dimensionless).
  * ALGEBRAIC[39] is tau_q in component i_to_q_gate (second).
  * ALGEBRAIC[22] is r_infinity in component i_to_r_gate (dimensionless).
  * ALGEBRAIC[40] is tau_r in component i_to_r_gate (second).
  * CONSTANTS_FS[87] is g_Kr in component i_Kr (microS).
  * STATES[28] is paS in component i_Kr_pa_gate (dimensionless).
  * STATES[29] is paF in component i_Kr_pa_gate (dimensionless).
  * STATES[30] is piy in component i_Kr_pi_gate (dimensionless).
  * ALGEBRAIC[23] is pa_infinity in component i_Kr_pa_gate (dimensionless).
  * ALGEBRAIC[24] is alfapaF in component i_Kr_pa_gate (per_second).
  * ALGEBRAIC[25] is betapaF in component i_Kr_pa_gate (per_second).
  * ALGEBRAIC[41] is tau_paS in component i_Kr_pa_gate (second).
  * ALGEBRAIC[42] is tau_paF in component i_Kr_pa_gate (second).
  * ALGEBRAIC[43] is pi_infinity in component i_Kr_pi_gate (dimensionless).
  * ALGEBRAIC[26] is tau_pi in component i_Kr_pi_gate (second).
  * CONSTANTS_FS[94] is g_Ks in component i_Ks (microS).
  * CONSTANTS_FS[88] is g_Ks_ in component i_Ks (microS).
  * ALGEBRAIC[92] is E_Ks in component i_Ks (millivolt).
  * STATES[31] is n in component i_Ks_n_gate (dimensionless).
  * ALGEBRAIC[27] is n_infinity in component i_Ks_n_gate (dimensionless).
  * ALGEBRAIC[56] is tau_n in component i_Ks_n_gate (second).
  * CONSTANTS_FS[114] is Iso_shift in component i_Ks_n_gate (millivolt).
  * ALGEBRAIC[44] is alpha_n in component i_Ks_n_gate (per_second).
  * ALGEBRAIC[50] is beta_n in component i_Ks_n_gate (per_second).
  * CONSTANTS_FS[89] is g_KACh in component i_KACh (microS).
  * STATES[32] is a in component i_KACh_a_gate (dimensionless).
  * CONSTANTS_FS[90] is ACh_on in component i_KACh (dimensionless).
  * CONSTANTS_FS[115] is alpha_a in component i_KACh_a_gate (per_second).
  * ALGEBRAIC[28] is beta_a in component i_KACh_a_gate (per_second).
  * ALGEBRAIC[45] is a_infinity in component i_KACh_a_gate (dimensionless).
  * ALGEBRAIC[51] is tau_a in component i_KACh_a_gate (second).
  * RATES[0] is d/dt V_ode in component Membrane (millivolt).
  * RATES[2] is d/dt Nai_ in component Nai_concentration (millimolar).
  * RATES[3] is d/dt y in component i_f_y_gate (dimensionless).
  * RATES[4] is d/dt m in component i_Na_m_gate (dimensionless).
  * RATES[5] is d/dt h in component i_Na_h_gate (dimensionless).
  * RATES[6] is d/dt dL in component i_CaL_dL_gate (dimensionless).
  * RATES[7] is d/dt fL in component i_CaL_fL_gate (dimensionless).
  * RATES[8] is d/dt fCa in component i_CaL_fCa_gate (dimensionless).
  * RATES[9] is d/dt dT in component i_CaT_dT_gate (dimensionless).
  * RATES[10] is d/dt fT in component i_CaT_fT_gate (dimensionless).
  * RATES[11] is d/dt R in component Ca_SR_release (dimensionless).
  * RATES[12] is d/dt O in component Ca_SR_release (dimensionless).
  * RATES[13] is d/dt I in component Ca_SR_release (dimensionless).
  * RATES[14] is d/dt RI in component Ca_SR_release (dimensionless).
  * RATES[21] is d/dt fTC in component Ca_buffering (dimensionless).
  * RATES[22] is d/dt fTMC in component Ca_buffering (dimensionless).
  * RATES[18] is d/dt fTMM in component Ca_buffering (dimensionless).
  * RATES[19] is d/dt fCMi in component Ca_buffering (dimensionless).
  * RATES[20] is d/dt fCMs in component Ca_buffering (dimensionless).
  * RATES[23] is d/dt fCQ in component Ca_buffering (dimensionless).
  * RATES[17] is d/dt Cai in component Ca_dynamics (millimolar).
  * RATES[1] is d/dt Ca_sub in component Ca_dynamics (millimolar).
  * RATES[16] is d/dt Ca_nsr in component Ca_dynamics (millimolar).
  * RATES[15] is d/dt Ca_jsr in component Ca_dynamics (millimolar).
  * RATES[24] is d/dt r_Kur in component i_Kur_rKur_gate (dimensionless).
  * RATES[25] is d/dt s_Kur in component i_Kur_sKur_gate (dimensionless).
  * RATES[26] is d/dt q in component i_to_q_gate (dimensionless).
  * RATES[27] is d/dt r in component i_to_r_gate (dimensionless).
  * RATES[28] is d/dt paS in component i_Kr_pa_gate (dimensionless).
  * RATES[29] is d/dt paF in component i_Kr_pa_gate (dimensionless).
  * RATES[30] is d/dt piy in component i_Kr_pi_gate (dimensionless).
  * RATES[31] is d/dt n in component i_Ks_n_gate (dimensionless).
  * RATES[32] is d/dt a in component i_KACh_a_gate (dimensionless).
  */
  PDEFIELD_TYPE CONSTANTS_FS[116];
  PDEFIELD_TYPE ALGEBRAIC[101];

  CONSTANTS_FS[0] = 8314.472;
  CONSTANTS_FS[1] = 310;
  CONSTANTS_FS[2] = 96485.3415;
  CONSTANTS_FS[3] = 5.7e-5;
  CONSTANTS_FS[4] = 0;
  CONSTANTS_FS[5] = 0.5;
  CONSTANTS_FS[6] = 0.5;
  CONSTANTS_FS[7] = -35;
  CONSTANTS_FS[8] = -45;
  CONSTANTS_FS[9] = 0;
  CONSTANTS_FS[10] = 0;
  CONSTANTS_FS[11] = 140;
  CONSTANTS_FS[12] = 140;
  CONSTANTS_FS[13] = 5.4;
  CONSTANTS_FS[14] = 1.8;
  CONSTANTS_FS[15] = 1;
  CONSTANTS_FS[16] = 0.00427;
  CONSTANTS_FS[17] = 45;
  CONSTANTS_FS[18] = 0.5927;
  CONSTANTS_FS[19] = 0;
  CONSTANTS_FS[20] = 0;
  CONSTANTS_FS[21] = 1.4;
  CONSTANTS_FS[22] = 14;
  CONSTANTS_FS[23] = 0.08105;
  CONSTANTS_FS[24] = 3.343;
  CONSTANTS_FS[25] = 0.1369;
  CONSTANTS_FS[26] = 0.4315;
  CONSTANTS_FS[27] = 0;
  CONSTANTS_FS[28] = 26.44;
  CONSTANTS_FS[29] = 0.0207;
  CONSTANTS_FS[30] = 395.3;
  CONSTANTS_FS[31] = 2.289;
  CONSTANTS_FS[32] = 26.44;
  CONSTANTS_FS[33] = 4.663;
  CONSTANTS_FS[34] = 1628;
  CONSTANTS_FS[35] = 561.4;
  CONSTANTS_FS[36] = 3.663;
  CONSTANTS_FS[37] = 0;
  CONSTANTS_FS[38] = 0.0223;
  CONSTANTS_FS[39] = 0;
  CONSTANTS_FS[40] = 1e-5;
  CONSTANTS_FS[41] = 0.4578;
  CONSTANTS_FS[42] = 4.3371;
  CONSTANTS_FS[43] = -16.4508;
  CONSTANTS_FS[44] = 0;
  CONSTANTS_FS[45] = 0;
  CONSTANTS_FS[46] = 0.0075;
  CONSTANTS_FS[47] = 0.000338;
  CONSTANTS_FS[48] = 0.04132;
  CONSTANTS_FS[49] = 0;
  CONSTANTS_FS[50] = 148041085.1;
  CONSTANTS_FS[51] = 15;
  CONSTANTS_FS[52] = 1;
  CONSTANTS_FS[53] = 0.45;
  CONSTANTS_FS[54] = 2.5;
  CONSTANTS_FS[55] = 10000;
  CONSTANTS_FS[56] = 500;
  CONSTANTS_FS[57] = 5;
  CONSTANTS_FS[58] = 660;
  CONSTANTS_FS[59] = 5.469e-5;
  CONSTANTS_FS[60] = 0.04;
  CONSTANTS_FS[61] = 5;
  CONSTANTS_FS[62] = 0.000286113;
  CONSTANTS_FS[63] = 5e-5;
  CONSTANTS_FS[64] = 0.031;
  CONSTANTS_FS[65] = 0.062;
  CONSTANTS_FS[66] = 0.045;
  CONSTANTS_FS[67] = 10;
  CONSTANTS_FS[68] = 88800;
  CONSTANTS_FS[69] = 2277;
  CONSTANTS_FS[70] = 227700;
  CONSTANTS_FS[71] = 1.642e6;
  CONSTANTS_FS[72] = 175.4;
  CONSTANTS_FS[73] = 446;
  CONSTANTS_FS[74] = 7.51;
  CONSTANTS_FS[75] = 751;
  CONSTANTS_FS[76] = 542;
  CONSTANTS_FS[77] = 445;
  CONSTANTS_FS[78] = 2.5;
  CONSTANTS_FS[79] = 0.0012;
  CONSTANTS_FS[80] = 0.46;
  CONSTANTS_FS[81] = 0.0116;
  CONSTANTS_FS[82] = 3.9;
  CONSTANTS_FS[83] = 67;
  CONSTANTS_FS[84] = 0.02;
  CONSTANTS_FS[85] = 0.1539e-3;
  CONSTANTS_FS[86] = 3.5e-3;
  CONSTANTS_FS[87] = 0.00424;
  CONSTANTS_FS[88] = 0.00065;
  CONSTANTS_FS[89] = 0.00345;
  CONSTANTS_FS[90] = 1;
  CONSTANTS_FS[91] = ( CONSTANTS_FS[0]*CONSTANTS_FS[1])/CONSTANTS_FS[2];
  CONSTANTS_FS[92] = CONSTANTS_FS[16]/(CONSTANTS_FS[13]/(CONSTANTS_FS[13]+CONSTANTS_FS[17]));
  CONSTANTS_FS[93] = (CONSTANTS_FS[10]>0.00000 ? - 0.250000 : CONSTANTS_FS[9]>0.00000 ? ( 0.700000*CONSTANTS_FS[9])/(9.00000e-05+CONSTANTS_FS[9]) : 0.00000);
  CONSTANTS_FS[94] = (CONSTANTS_FS[10]>0.00000 ?  1.20000*CONSTANTS_FS[88] : CONSTANTS_FS[88]);
  CONSTANTS_FS[95] = (CONSTANTS_FS[9]>0.00000 ? - 1.00000 - ( 9.89800*pow( 1.00000*CONSTANTS_FS[9], 0.618000))/(pow( 1.00000*CONSTANTS_FS[9], 0.618000)+0.00122423) : 0.00000);
  CONSTANTS_FS[96] =  CONSTANTS_FS[91]*log(CONSTANTS_FS[13]/CONSTANTS_FS[12]);
  CONSTANTS_FS[97] = CONSTANTS_FS[92]/(CONSTANTS_FS[18]+1.00000);
  CONSTANTS_FS[98] =  CONSTANTS_FS[61]*(1.00000 - CONSTANTS_FS[93]);
  CONSTANTS_FS[99] = (CONSTANTS_FS[10]>0.00000 ? 7.50000 : 0.00000);
  CONSTANTS_FS[100] =  CONSTANTS_FS[18]*CONSTANTS_FS[97];
  CONSTANTS_FS[101] = ( CONSTANTS_FS[97]*CONSTANTS_FS[13])/(CONSTANTS_FS[13]+CONSTANTS_FS[17]);
  CONSTANTS_FS[102] = (CONSTANTS_FS[10]>0.00000 ? 1.20000 : 1.00000);
  CONSTANTS_FS[103] = ( CONSTANTS_FS[100]*CONSTANTS_FS[13])/(CONSTANTS_FS[13]+CONSTANTS_FS[17]);
  CONSTANTS_FS[104] = CONSTANTS_FS[11]/(CONSTANTS_FS[33]+CONSTANTS_FS[11]);
  CONSTANTS_FS[105] = (CONSTANTS_FS[10]>0.00000 ? 1.23000 : 1.00000);
  CONSTANTS_FS[106] = ( 0.310000*CONSTANTS_FS[9])/(CONSTANTS_FS[9]+9.00000e-05);
  CONSTANTS_FS[107] = (CONSTANTS_FS[10]>0.00000 ? - 8.00000 : 0.00000);
  CONSTANTS_FS[108] = (CONSTANTS_FS[10]>0.00000 ? - 27.0000 : 0.00000);
  CONSTANTS_FS[109] =  1.00000e-09* 3.14159265358979*pow(CONSTANTS_FS[82], 2.00000)*CONSTANTS_FS[83];
  CONSTANTS_FS[110] =  1.00000e-09*2.00000* 3.14159265358979*CONSTANTS_FS[84]*(CONSTANTS_FS[82] - CONSTANTS_FS[84]/2.00000)*CONSTANTS_FS[83];
  CONSTANTS_FS[111] =  CONSTANTS_FS[79]*CONSTANTS_FS[109];
  CONSTANTS_FS[112] =  CONSTANTS_FS[80]*CONSTANTS_FS[109] - CONSTANTS_FS[110];
  CONSTANTS_FS[113] =  CONSTANTS_FS[81]*CONSTANTS_FS[109];
  CONSTANTS_FS[114] = (CONSTANTS_FS[10]>0.00000 ? - 14.0000 : 0.00000);
  CONSTANTS_FS[115] = (3.59880 - 0.0256410)/(1.00000+1.21550e-06/pow( 1.00000*CONSTANTS_FS[9], 1.69510))+0.0256410;


  ALGEBRAIC[6] =  CONSTANTS_FS[69]*CONSTANTS_FS[78]*(1.00000 - (STATES[22]+STATES[18])) -  CONSTANTS_FS[75]*STATES[18];
  RATES[18] = ALGEBRAIC[6];
  ALGEBRAIC[1] = CONSTANTS_FS[47]/(CONSTANTS_FS[47]+STATES[1]);
  ALGEBRAIC[7] = ( 0.00100000*ALGEBRAIC[1])/CONSTANTS_FS[46];
  RATES[8] = (ALGEBRAIC[1] - STATES[8])/ALGEBRAIC[7];
  STATES_NEW[8] = ALGEBRAIC[1] + (STATES[8] - ALGEBRAIC[1]) * exp(-dt/ALGEBRAIC[7]);
  ALGEBRAIC[2] = CONSTANTS_FS[51] - (CONSTANTS_FS[51] - CONSTANTS_FS[52])/(1.00000+pow(CONSTANTS_FS[53]/STATES[15], CONSTANTS_FS[54]));
  ALGEBRAIC[8] = CONSTANTS_FS[55]/ALGEBRAIC[2];
  ALGEBRAIC[17] =  CONSTANTS_FS[56]*ALGEBRAIC[2];
  RATES[11] = ( CONSTANTS_FS[57]*STATES[14] -  ALGEBRAIC[17]*STATES[1]*STATES[11]) - ( ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[11] -  CONSTANTS_FS[58]*STATES[12]);
  RATES[12] = ( ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[11] -  CONSTANTS_FS[58]*STATES[12]) - ( ALGEBRAIC[17]*STATES[1]*STATES[12] -  CONSTANTS_FS[57]*STATES[13]);
  RATES[13] = ( ALGEBRAIC[17]*STATES[1]*STATES[12] -  CONSTANTS_FS[57]*STATES[13]) - ( CONSTANTS_FS[58]*STATES[13] -  ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[14]);
  RATES[14] = ( CONSTANTS_FS[58]*STATES[13] -  ALGEBRAIC[8]*pow(STATES[1], 2.00000)*STATES[14]) - ( CONSTANTS_FS[57]*STATES[14] -  ALGEBRAIC[17]*STATES[1]*STATES[11]);
  ALGEBRAIC[5] = (VOI>CONSTANTS_FS[5]&&VOI<CONSTANTS_FS[5]+CONSTANTS_FS[6] ? CONSTANTS_FS[7] : CONSTANTS_FS[8]);
  ALGEBRAIC[9] = (CONSTANTS_FS[4]>=1.00000 ? ALGEBRAIC[5] : STATES[0]);
  ALGEBRAIC[10] = 1.00000/(( 0.360000*(((ALGEBRAIC[9]+148.800) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]))/(exp( 0.0660000*(((ALGEBRAIC[9]+148.800) - CONSTANTS_FS[95]) - CONSTANTS_FS[99])) - 1.00000)+( 0.100000*(((ALGEBRAIC[9]+87.3000) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]))/(1.00000 - exp( - 0.200000*(((ALGEBRAIC[9]+87.3000) - CONSTANTS_FS[95]) - CONSTANTS_FS[99])))) - 0.0540000;
  ALGEBRAIC[29] = (ALGEBRAIC[9]<- (((80.0000 - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20]) ? 0.0132900+0.999210/(1.00000+exp(((((ALGEBRAIC[9]+97.1340) - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20])/8.17520)) :  0.000250100*exp(- (((ALGEBRAIC[9] - CONSTANTS_FS[95]) - CONSTANTS_FS[99]) - CONSTANTS_FS[20])/12.8610));
  RATES[3] = (ALGEBRAIC[29] - STATES[3])/ALGEBRAIC[10];
  STATES_NEW[3] = ALGEBRAIC[29] + (STATES[3] - ALGEBRAIC[29]) * exp(-dt/ALGEBRAIC[10]);
  ALGEBRAIC[14] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+37.4000+CONSTANTS_FS[44])/(5.30000+CONSTANTS_FS[45])));
  ALGEBRAIC[33] =  0.00100000*(44.3000+ 230.000*exp(- pow((ALGEBRAIC[9]+36.0000)/10.0000, 2.00000)));
  RATES[7] = (ALGEBRAIC[14] - STATES[7])/ALGEBRAIC[33];
  STATES_NEW[7] = ALGEBRAIC[14] + (STATES[7] - ALGEBRAIC[14]) * exp(-dt/ALGEBRAIC[33]);
  ALGEBRAIC[15] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+38.3000)/5.50000));
  ALGEBRAIC[34] = 0.00100000/( 1.06800*exp((ALGEBRAIC[9]+38.3000)/30.0000)+ 1.06800*exp(- (ALGEBRAIC[9]+38.3000)/30.0000));
  RATES[9] = (ALGEBRAIC[15] - STATES[9])/ALGEBRAIC[34];
  STATES_NEW[9] = ALGEBRAIC[15] + (STATES[9] - ALGEBRAIC[15]) * exp(-dt/ALGEBRAIC[34]);
  ALGEBRAIC[16] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+58.7000)/3.80000));
  ALGEBRAIC[35] = 1.00000/( 16.6700*exp(- (ALGEBRAIC[9]+75.0000)/83.3000)+ 16.6700*exp((ALGEBRAIC[9]+75.0000)/15.3800))+CONSTANTS_FS[49];
  RATES[10] = (ALGEBRAIC[16] - STATES[10])/ALGEBRAIC[35];
  STATES_NEW[10] = ALGEBRAIC[16] + (STATES[10] - ALGEBRAIC[16]) * exp(-dt/ALGEBRAIC[35]);
  ALGEBRAIC[37] = 0.00900000/(1.00000+exp((ALGEBRAIC[9]+5.00000)/12.0000))+0.000500000;
  ALGEBRAIC[19] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+6.00000)/- 8.60000));
  RATES[24] = (ALGEBRAIC[19] - STATES[24])/ALGEBRAIC[37];
  STATES_NEW[24] = ALGEBRAIC[19] + (STATES[24] - ALGEBRAIC[19]) * exp(-dt/ALGEBRAIC[37]);
  ALGEBRAIC[38] = 0.590000/(1.00000+exp((ALGEBRAIC[9]+60.0000)/10.0000))+3.05000;
  ALGEBRAIC[20] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+7.50000)/10.0000));
  RATES[25] = (ALGEBRAIC[20] - STATES[25])/ALGEBRAIC[38];
  STATES_NEW[25] = ALGEBRAIC[20] + (STATES[25] - ALGEBRAIC[20]) * exp(-dt/ALGEBRAIC[38]);
  ALGEBRAIC[21] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+49.0000)/13.0000));
  ALGEBRAIC[39] =  0.00100000*0.600000*(65.1700/( 0.570000*exp( - 0.0800000*(ALGEBRAIC[9]+44.0000))+ 0.0650000*exp( 0.100000*(ALGEBRAIC[9]+45.9300)))+10.1000);
  RATES[26] = (ALGEBRAIC[21] - STATES[26])/ALGEBRAIC[39];
  STATES_NEW[26] = ALGEBRAIC[21] + (STATES[26] - ALGEBRAIC[21]) * exp(-dt/ALGEBRAIC[39]);
  ALGEBRAIC[22] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9] - 19.3000)/15.0000));
  ALGEBRAIC[40] =  0.00100000*0.660000*1.40000*(15.5900/( 1.03700*exp( 0.0900000*(ALGEBRAIC[9]+30.6100))+ 0.369000*exp( - 0.120000*(ALGEBRAIC[9]+23.8400)))+2.98000);
  RATES[27] = (ALGEBRAIC[22] - STATES[27])/ALGEBRAIC[40];
  STATES_NEW[27] = ALGEBRAIC[22] + (STATES[27] - ALGEBRAIC[22]) * exp(-dt/ALGEBRAIC[40]);
  ALGEBRAIC[23] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+10.0144)/7.66070));
  ALGEBRAIC[41] = 0.846554/( 4.20000*exp(ALGEBRAIC[9]/17.0000)+ 0.150000*exp(- ALGEBRAIC[9]/21.6000));
  RATES[28] = (ALGEBRAIC[23] - STATES[28])/ALGEBRAIC[41];
  STATES_NEW[28] = ALGEBRAIC[23] + (STATES[28] - ALGEBRAIC[23]) * exp(-dt/ALGEBRAIC[41]);
  ALGEBRAIC[42] = 1.00000/( 30.0000*exp(ALGEBRAIC[9]/10.0000)+exp(- ALGEBRAIC[9]/12.0000));
  RATES[29] = (ALGEBRAIC[23] - STATES[29])/ALGEBRAIC[42];
  STATES_NEW[29] = ALGEBRAIC[23] + (STATES[29] - ALGEBRAIC[23]) * exp(-dt/ALGEBRAIC[42]);
  ALGEBRAIC[43] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+28.6000)/17.1000));
  ALGEBRAIC[26] = 1.00000/( 100.000*exp(- ALGEBRAIC[9]/54.6450)+ 656.000*exp(ALGEBRAIC[9]/106.157));
  RATES[30] = (ALGEBRAIC[43] - STATES[30])/ALGEBRAIC[26];
  STATES_NEW[30] = ALGEBRAIC[43] + (STATES[30] - ALGEBRAIC[43]) * exp(-dt/ALGEBRAIC[26]);
  ALGEBRAIC[28] =  10.0000*exp( 0.0133000*(ALGEBRAIC[9]+40.0000));
  ALGEBRAIC[45] = CONSTANTS_FS[115]/(CONSTANTS_FS[115]+ALGEBRAIC[28]);
  ALGEBRAIC[51] = 1.00000/(CONSTANTS_FS[115]+ALGEBRAIC[28]);
  RATES[32] = (ALGEBRAIC[45] - STATES[32])/ALGEBRAIC[51];
  STATES_NEW[32] = ALGEBRAIC[45] + (STATES[32] - ALGEBRAIC[45]) * exp(-dt/ALGEBRAIC[51]);
  ALGEBRAIC[12] = 1.00000/(1.00000+exp((ALGEBRAIC[9]+69.8040)/4.45650));
  ALGEBRAIC[31] =  20.0000*exp( - 0.125000*(ALGEBRAIC[9]+75.0000));
  ALGEBRAIC[47] = 2000.00/( 320.000*exp( - 0.100000*(ALGEBRAIC[9]+75.0000))+1.00000);
  ALGEBRAIC[53] = 1.00000/(ALGEBRAIC[31]+ALGEBRAIC[47]);
  RATES[5] = (ALGEBRAIC[12] - STATES[5])/ALGEBRAIC[53];
  STATES_NEW[5] = ALGEBRAIC[12] + (STATES[5] - ALGEBRAIC[12]) * exp(-dt/ALGEBRAIC[53]);
  ALGEBRAIC[27] =  pow((1.00000/(1.00000+exp(- ((ALGEBRAIC[9]+0.638300) - CONSTANTS_FS[114])/10.7071))), 1.0 / 2);
  ALGEBRAIC[44] = 28.0000/(1.00000+exp(- ((ALGEBRAIC[9] - 40.0000) - CONSTANTS_FS[114])/3.00000));
  ALGEBRAIC[50] =  1.00000*exp(- ((ALGEBRAIC[9] - CONSTANTS_FS[114]) - 5.00000)/25.0000);
  ALGEBRAIC[56] = 1.00000/(ALGEBRAIC[44]+ALGEBRAIC[50]);
  RATES[31] = (ALGEBRAIC[27] - STATES[31])/ALGEBRAIC[56];
  STATES_NEW[31] = ALGEBRAIC[27] + (STATES[31] - ALGEBRAIC[27]) * exp(-dt/ALGEBRAIC[56]);
  ALGEBRAIC[11] = 1.00000/(1.00000+exp(- (ALGEBRAIC[9]+42.0504)/8.31060));
  ALGEBRAIC[30] = ALGEBRAIC[9]+41.0000;
  ALGEBRAIC[46] = (fabs(ALGEBRAIC[30])<CONSTANTS_FS[40] ? 2000.00 : ( 200.000*ALGEBRAIC[30])/(1.00000 - exp( - 0.100000*ALGEBRAIC[30])));
  ALGEBRAIC[52] =  8000.00*exp( - 0.0560000*(ALGEBRAIC[9]+66.0000));
  ALGEBRAIC[57] = 1.00000/(ALGEBRAIC[46]+ALGEBRAIC[52]);
  RATES[4] = (ALGEBRAIC[11] - STATES[4])/ALGEBRAIC[57];
  STATES_NEW[4] = ALGEBRAIC[11] + (STATES[4] - ALGEBRAIC[11]) * exp(-dt/ALGEBRAIC[57]);
  ALGEBRAIC[13] = 1.00000/(1.00000+exp(- ((ALGEBRAIC[9] - CONSTANTS_FS[43]) - CONSTANTS_FS[107])/( CONSTANTS_FS[42]*(1.00000+CONSTANTS_FS[108]/100.000))));
  ALGEBRAIC[32] = (ALGEBRAIC[9]==- 41.8000 ? - 41.8000 : ALGEBRAIC[9]==0.00000 ? 0.00000 : ALGEBRAIC[9]==- 6.80000 ? - 6.80001 : ALGEBRAIC[9]);
  ALGEBRAIC[48] = ( - 0.0283900*(ALGEBRAIC[32]+41.8000))/(exp(- (ALGEBRAIC[32]+41.8000)/2.50000) - 1.00000) - ( 0.0849000*(ALGEBRAIC[32]+6.80000))/(exp(- (ALGEBRAIC[32]+6.80000)/4.80000) - 1.00000);
  ALGEBRAIC[54] = (ALGEBRAIC[9]==- 1.80000 ? - 1.80001 : ALGEBRAIC[9]);
  ALGEBRAIC[58] = ( 0.0114300*(ALGEBRAIC[54]+1.80000))/(exp((ALGEBRAIC[54]+1.80000)/2.50000) - 1.00000);
  ALGEBRAIC[60] = 0.00100000/(ALGEBRAIC[48]+ALGEBRAIC[58]);
  RATES[6] = (ALGEBRAIC[13] - STATES[6])/ALGEBRAIC[60];
  STATES_NEW[6] = ALGEBRAIC[13] + (STATES[6] - ALGEBRAIC[13]) * exp(-dt/ALGEBRAIC[60]);
  ALGEBRAIC[18] = STATES[2];
  ALGEBRAIC[36] =  CONSTANTS_FS[91]*log(CONSTANTS_FS[11]/ALGEBRAIC[18]);
  ALGEBRAIC[61] =  CONSTANTS_FS[102]*CONSTANTS_FS[23]*pow(1.00000+pow(CONSTANTS_FS[21]/CONSTANTS_FS[13], 1.20000), - 1.00000)*pow(1.00000+pow(CONSTANTS_FS[22]/ALGEBRAIC[18], 1.30000), - 1.00000)*pow(1.00000+exp(- ((ALGEBRAIC[9] - ALGEBRAIC[36])+110.000)/20.0000), - 1.00000);
  ALGEBRAIC[63] = exp(( - CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91]));
  ALGEBRAIC[69] = 1.00000+ (CONSTANTS_FS[14]/CONSTANTS_FS[36])*(1.00000+exp(( CONSTANTS_FS[27]*ALGEBRAIC[9])/CONSTANTS_FS[91]))+ (CONSTANTS_FS[11]/CONSTANTS_FS[34])*(1.00000+ (CONSTANTS_FS[11]/CONSTANTS_FS[35])*(1.00000+CONSTANTS_FS[11]/CONSTANTS_FS[33]));
  ALGEBRAIC[71] = ( (( (CONSTANTS_FS[11]/CONSTANTS_FS[34])*CONSTANTS_FS[11])/CONSTANTS_FS[35])*(1.00000+CONSTANTS_FS[11]/CONSTANTS_FS[33])*exp(( - CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91])))/ALGEBRAIC[69];
  ALGEBRAIC[70] = ( (CONSTANTS_FS[14]/CONSTANTS_FS[36])*exp(( CONSTANTS_FS[27]*ALGEBRAIC[9])/CONSTANTS_FS[91]))/ALGEBRAIC[69];
  ALGEBRAIC[67] = exp(( CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91]));
  ALGEBRAIC[62] = ALGEBRAIC[18]/(CONSTANTS_FS[28]+ALGEBRAIC[18]);
  ALGEBRAIC[72] =  ALGEBRAIC[63]*CONSTANTS_FS[104]*(ALGEBRAIC[71]+ALGEBRAIC[70])+ ALGEBRAIC[70]*ALGEBRAIC[67]*(ALGEBRAIC[62]+ALGEBRAIC[63]);
  ALGEBRAIC[64] = 1.00000+ (STATES[1]/CONSTANTS_FS[29])*(1.00000+exp(( - CONSTANTS_FS[25]*ALGEBRAIC[9])/CONSTANTS_FS[91])+ALGEBRAIC[18]/CONSTANTS_FS[32])+ (ALGEBRAIC[18]/CONSTANTS_FS[30])*(1.00000+ (ALGEBRAIC[18]/CONSTANTS_FS[31])*(1.00000+ALGEBRAIC[18]/CONSTANTS_FS[28]));
  ALGEBRAIC[65] = ( (STATES[1]/CONSTANTS_FS[29])*exp(( - CONSTANTS_FS[25]*ALGEBRAIC[9])/CONSTANTS_FS[91]))/ALGEBRAIC[64];
  ALGEBRAIC[66] = ( (( (ALGEBRAIC[18]/CONSTANTS_FS[30])*ALGEBRAIC[18])/CONSTANTS_FS[31])*(1.00000+ALGEBRAIC[18]/CONSTANTS_FS[28])*exp(( CONSTANTS_FS[26]*ALGEBRAIC[9])/( 2.00000*CONSTANTS_FS[91])))/ALGEBRAIC[64];
  ALGEBRAIC[68] =  ALGEBRAIC[67]*ALGEBRAIC[62]*(ALGEBRAIC[66]+ALGEBRAIC[65])+ ALGEBRAIC[63]*ALGEBRAIC[65]*(CONSTANTS_FS[104]+ALGEBRAIC[67]);
  ALGEBRAIC[73] =  ALGEBRAIC[66]*ALGEBRAIC[62]*(ALGEBRAIC[71]+ALGEBRAIC[70])+ ALGEBRAIC[65]*ALGEBRAIC[71]*(ALGEBRAIC[62]+ALGEBRAIC[63]);
  ALGEBRAIC[74] =  ALGEBRAIC[71]*CONSTANTS_FS[104]*(ALGEBRAIC[66]+ALGEBRAIC[65])+ ALGEBRAIC[66]*ALGEBRAIC[70]*(CONSTANTS_FS[104]+ALGEBRAIC[67]);
  ALGEBRAIC[75] = ( (1.00000 - CONSTANTS_FS[37])*CONSTANTS_FS[24]*( ALGEBRAIC[68]*ALGEBRAIC[70] -  ALGEBRAIC[72]*ALGEBRAIC[65]))/(ALGEBRAIC[72]+ALGEBRAIC[68]+ALGEBRAIC[73]+ALGEBRAIC[74]);
  ALGEBRAIC[76] =  CONSTANTS_FS[91]*log((CONSTANTS_FS[11]+ 0.120000*CONSTANTS_FS[13])/(ALGEBRAIC[18]+ 0.120000*CONSTANTS_FS[12]));
  ALGEBRAIC[77] =  CONSTANTS_FS[38]*pow(STATES[4], 3.00000)*STATES[5]*(ALGEBRAIC[9] - ALGEBRAIC[76]);
  ALGEBRAIC[78] =  CONSTANTS_FS[39]*pow(STATES[4], 3.00000)*(ALGEBRAIC[9] - ALGEBRAIC[76]);
  ALGEBRAIC[79] = ALGEBRAIC[77]+ALGEBRAIC[78];
  ALGEBRAIC[49] =  STATES[3]*CONSTANTS_FS[103]*(ALGEBRAIC[9] - ALGEBRAIC[36])*(1.00000 - CONSTANTS_FS[19]);
  ALGEBRAIC[82] =  (( 1.85000e-05*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))))*(ALGEBRAIC[18] -  CONSTANTS_FS[11]*exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  RATES[2] = ( (1.00000 - CONSTANTS_FS[15])*- 1.00000*(ALGEBRAIC[79]+ALGEBRAIC[49]+ALGEBRAIC[82]+ 3.00000*ALGEBRAIC[61]+ 3.00000*ALGEBRAIC[75]))/( 1.00000*(CONSTANTS_FS[112]+CONSTANTS_FS[110])*CONSTANTS_FS[2]);
  ALGEBRAIC[90] =  CONSTANTS_FS[71]*STATES[1]*(1.00000 - STATES[20]) -  CONSTANTS_FS[76]*STATES[20];
  RATES[20] = ALGEBRAIC[90];
  ALGEBRAIC[84] =  (( 2.00000*CONSTANTS_FS[48]*ALGEBRAIC[9])/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*ALGEBRAIC[9]*2.00000)/CONSTANTS_FS[91]))))*(STATES[1] -  CONSTANTS_FS[14]*exp(( - 2.00000*ALGEBRAIC[9])/CONSTANTS_FS[91]))*STATES[9]*STATES[10];
  ALGEBRAIC[80] =  (( 2.00000*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000)*2.00000)/CONSTANTS_FS[91]))))*(STATES[1] -  CONSTANTS_FS[14]*exp(( - 2.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  ALGEBRAIC[86] =  CONSTANTS_FS[50]*STATES[12]*(STATES[15] - STATES[1]);
  ALGEBRAIC[88] = (STATES[1] - STATES[17])/CONSTANTS_FS[59];
  RATES[1] = ( ALGEBRAIC[86]*CONSTANTS_FS[111])/CONSTANTS_FS[110] - (((ALGEBRAIC[80]+ALGEBRAIC[84]) -  2.00000*ALGEBRAIC[75])/( 2.00000*CONSTANTS_FS[2]*CONSTANTS_FS[110])+ALGEBRAIC[88]+ CONSTANTS_FS[66]*ALGEBRAIC[90]);
  ALGEBRAIC[93] =  CONSTANTS_FS[68]*STATES[17]*(1.00000 - STATES[21]) -  CONSTANTS_FS[73]*STATES[21];
  RATES[21] = ALGEBRAIC[93];
  ALGEBRAIC[91] = CONSTANTS_FS[98]/(1.00000+exp((- STATES[17]+CONSTANTS_FS[62])/CONSTANTS_FS[63]));
  ALGEBRAIC[94] = (STATES[16] - STATES[15])/CONSTANTS_FS[60];
  RATES[16] = ALGEBRAIC[91] - ( ALGEBRAIC[94]*CONSTANTS_FS[111])/CONSTANTS_FS[113];
  ALGEBRAIC[96] =  CONSTANTS_FS[70]*STATES[17]*(1.00000 - (STATES[22]+STATES[18])) -  CONSTANTS_FS[74]*STATES[22];
  RATES[22] = ALGEBRAIC[96];
  ALGEBRAIC[97] =  CONSTANTS_FS[72]*STATES[15]*(1.00000 - STATES[23]) -  CONSTANTS_FS[77]*STATES[23];
  RATES[23] = ALGEBRAIC[97];
  RATES[15] = ALGEBRAIC[94] - (ALGEBRAIC[86]+ CONSTANTS_FS[67]*ALGEBRAIC[97]);
  ALGEBRAIC[99] =  CONSTANTS_FS[71]*STATES[17]*(1.00000 - STATES[19]) -  CONSTANTS_FS[76]*STATES[19];
  RATES[19] = ALGEBRAIC[99];
  RATES[17] = ( 1.00000*( ALGEBRAIC[88]*CONSTANTS_FS[110] -  ALGEBRAIC[91]*CONSTANTS_FS[113]))/CONSTANTS_FS[112] - ( CONSTANTS_FS[66]*ALGEBRAIC[99]+ CONSTANTS_FS[64]*ALGEBRAIC[93]+ CONSTANTS_FS[65]*ALGEBRAIC[96]);
  ALGEBRAIC[55] =  STATES[3]*CONSTANTS_FS[101]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*(1.00000 - CONSTANTS_FS[19]);
  ALGEBRAIC[59] = ALGEBRAIC[49]+ALGEBRAIC[55];
  ALGEBRAIC[89] =  CONSTANTS_FS[87]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*( 0.900000*STATES[29]+ 0.100000*STATES[28])*STATES[30];
  ALGEBRAIC[92] =  CONSTANTS_FS[91]*log((CONSTANTS_FS[13]+ 0.120000*CONSTANTS_FS[11])/(CONSTANTS_FS[12]+ 0.120000*ALGEBRAIC[18]));
  ALGEBRAIC[95] =  CONSTANTS_FS[94]*(ALGEBRAIC[9] - ALGEBRAIC[92])*pow(STATES[31], 2.00000);
  ALGEBRAIC[87] =  CONSTANTS_FS[86]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*STATES[26]*STATES[27];
  ALGEBRAIC[81] =  (( 0.000365000*CONSTANTS_FS[41]*(ALGEBRAIC[9] - 0.00000))/( CONSTANTS_FS[91]*(1.00000 - exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))))*(CONSTANTS_FS[12] -  CONSTANTS_FS[13]*exp(( - 1.00000*(ALGEBRAIC[9] - 0.00000))/CONSTANTS_FS[91]))*STATES[6]*STATES[7]*STATES[8];
  ALGEBRAIC[83] =  (ALGEBRAIC[80]+ALGEBRAIC[81]+ALGEBRAIC[82])*(1.00000 - CONSTANTS_FS[106])*1.00000*CONSTANTS_FS[105];
  ALGEBRAIC[98] = (CONSTANTS_FS[9]>0.00000 ?  CONSTANTS_FS[90]*CONSTANTS_FS[89]*(ALGEBRAIC[9] - CONSTANTS_FS[96])*(1.00000+exp((ALGEBRAIC[9]+20.0000)/20.0000))*STATES[32] : 0.00000);
  ALGEBRAIC[85] =  CONSTANTS_FS[85]*STATES[24]*STATES[25]*(ALGEBRAIC[9] - CONSTANTS_FS[96]);
  ALGEBRAIC[100] = I_f_factor*ALGEBRAIC[59]+I_Kr_factor*ALGEBRAIC[89]+ALGEBRAIC[95]+ALGEBRAIC[87]+ALGEBRAIC[61]+ALGEBRAIC[75]+ALGEBRAIC[79]+ALGEBRAIC[83]+ALGEBRAIC[84]+ALGEBRAIC[98]+ALGEBRAIC[85];
  RATES[0] = - ALGEBRAIC[100]/CONSTANTS_FS[3];

  for (int i = 0; i < 33; i++)
    if (i != 3 and i != 4 and i != 5 and i != 6 and i != 7 and i != 8 and i != 9 and i != 10 and i != 24 and i != 25 and i != 26 and i != 27 and i != 28 and i != 29 and i != 30 and i != 31 and i != 32)
      STATES_NEW[i] = STATES[i] + dt*RATES[i];

}

__device__ void derivsMaleckar(PDEFIELD_TYPE VOI, PDEFIELD_TYPE* STATES,PDEFIELD_TYPE* RATES, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE I_Na_factor, int id, PDEFIELD_TYPE activation_strength){


  /*
    There are a total of 70 entries in the algebraic variable array.
    There are a total of 30 entries in each of the rate and state variable arrays.
    There are a total of 51 entries in the constant variable array.
  */
  /*
  * VOI is time in component environment (second).
  * STATES[0] is V in component membrane (millivolt).
  * CONSTANTS_M[0] is R in component membrane (millijoule_per_mole_kelvin).
  * CONSTANTS_M[1] is T in component membrane (kelvin).
  * CONSTANTS_M[2] is F in component membrane (coulomb_per_mole).
  * CONSTANTS_M[3] is Cm in component membrane (nanoF).
  * ALGEBRAIC[0] is Q_tot in component membrane (millivolt).
  * ALGEBRAIC[37] is i_Na in component sodium_current (picoA).
  * ALGEBRAIC[41] is i_Ca_L in component L_type_Ca_channel (picoA).
  * ALGEBRAIC[44] is i_t in component Ca_independent_transient_outward_K_current (picoA).
  * ALGEBRAIC[45] is i_Kur in component ultra_rapid_K_current (picoA).
  * ALGEBRAIC[46] is i_K1 in component inward_rectifier (picoA).
  * ALGEBRAIC[49] is i_Kr in component delayed_rectifier_K_currents (picoA).
  * ALGEBRAIC[47] is i_Ks in component delayed_rectifier_K_currents (picoA).
  * ALGEBRAIC[50] is i_B_Na in component background_currents (picoA).
  * ALGEBRAIC[52] is i_B_Ca in component background_currents (picoA).
  * ALGEBRAIC[54] is i_NaK in component sodium_potassium_pump (picoA).
  * ALGEBRAIC[55] is i_CaP in component sarcolemmal_calcium_pump_current (picoA).
  * ALGEBRAIC[56] is i_NaCa in component Na_Ca_ion_exchanger_current (picoA).
  * ALGEBRAIC[57] is i_KACh in component ACh_dependent_K_current (picoA).
  * ALGEBRAIC[59] is I in component membrane (pA_per_nF).
  * ALGEBRAIC[24] is i_Stim in component membrane (pA_per_nF).
  * CONSTANTS_M[4] is stim_offset in component membrane (second).
  * CONSTANTS_M[5] is stim_period in component membrane (second).
  * CONSTANTS_M[6] is stim_duration in component membrane (second).
  * CONSTANTS_M[7] is stim_amplitude in component membrane (pA_per_nF).
  * ALGEBRAIC[1] is past in component membrane (second).
  * ALGEBRAIC[35] is E_Na in component sodium_current (millivolt).
  * CONSTANTS_M[8] is P_Na in component sodium_current (nanolitre_per_second).
  * STATES[1] is Na_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[2] is Na_i in component intracellular_ion_concentrations (millimolar).
  * STATES[3] is m in component sodium_current_m_gate (dimensionless).
  * STATES[4] is h1 in component sodium_current_h1_gate (dimensionless).
  * STATES[5] is h2 in component sodium_current_h2_gate (dimensionless).
  * ALGEBRAIC[14] is m_infinity in component sodium_current_m_gate (dimensionless).
  * ALGEBRAIC[2] is m_factor in component sodium_current_m_gate (dimensionless).
  * ALGEBRAIC[26] is tau_m in component sodium_current_m_gate (second).
  * ALGEBRAIC[3] is h_infinity in component sodium_current_h1_gate (dimensionless).
  * ALGEBRAIC[15] is h_factor in component sodium_current_h1_gate (dimensionless).
  * ALGEBRAIC[27] is tau_h1 in component sodium_current_h1_gate (second).
  * ALGEBRAIC[28] is tau_h2 in component sodium_current_h2_gate (second).
  * CONSTANTS_M[9] is g_Ca_L in component L_type_Ca_channel (nanoS).
  * CONSTANTS_M[10] is E_Ca_app in component L_type_Ca_channel (millivolt).
  * ALGEBRAIC[39] is f_Ca in component L_type_Ca_channel (dimensionless).
  * CONSTANTS_M[11] is k_Ca in component L_type_Ca_channel (millimolar).
  * STATES[6] is Ca_d in component intracellular_ion_concentrations (millimolar).
  * STATES[7] is d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
  * STATES[8] is f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * STATES[9] is f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
  * ALGEBRAIC[4] is d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless).
  * ALGEBRAIC[16] is d_L_factor in component L_type_Ca_channel_d_L_gate (dimensionless).
  * ALGEBRAIC[29] is tau_d_L in component L_type_Ca_channel_d_L_gate (second).
  * ALGEBRAIC[5] is f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * ALGEBRAIC[17] is f_L_factor in component L_type_Ca_channel_f_L1_gate (millivolt).
  * ALGEBRAIC[30] is tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second).
  * ALGEBRAIC[31] is tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second).
  * ALGEBRAIC[43] is E_K in component Ca_independent_transient_outward_K_current (millivolt).
  * CONSTANTS_M[12] is g_t in component Ca_independent_transient_outward_K_current (nanoS).
  * STATES[10] is K_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[11] is K_i in component intracellular_ion_concentrations (millimolar).
  * STATES[12] is r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * STATES[13] is s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * ALGEBRAIC[18] is tau_r in component Ca_independent_transient_outward_K_current_r_gate (second).
  * ALGEBRAIC[6] is r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * ALGEBRAIC[32] is tau_s in component Ca_independent_transient_outward_K_current_s_gate (second).
  * ALGEBRAIC[7] is s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * ALGEBRAIC[19] is s_factor in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * CONSTANTS_M[13] is g_kur in component ultra_rapid_K_current (nanoS).
  * STATES[14] is a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
  * STATES[15] is i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
  * ALGEBRAIC[8] is a_ur_infinity in component ultra_rapid_K_current_aur_gate (dimensionless).
  * ALGEBRAIC[20] is tau_a_ur in component ultra_rapid_K_current_aur_gate (second).
  * ALGEBRAIC[9] is i_ur_infinity in component ultra_rapid_K_current_iur_gate (dimensionless).
  * ALGEBRAIC[21] is tau_i_ur in component ultra_rapid_K_current_iur_gate (second).
  * CONSTANTS_M[14] is g_K1 in component inward_rectifier (nanoS).
  * CONSTANTS_M[15] is g_Ks in component delayed_rectifier_K_currents (nanoS).
  * CONSTANTS_M[16] is g_Kr in component delayed_rectifier_K_currents (nanoS).
  * STATES[16] is n in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * STATES[17] is pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * ALGEBRAIC[48] is pip in component delayed_rectifier_K_currents_pi_gate (dimensionless).
  * ALGEBRAIC[33] is tau_n in component delayed_rectifier_K_currents_n_gate (second).
  * ALGEBRAIC[10] is n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * ALGEBRAIC[22] is n_factor in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * ALGEBRAIC[34] is tau_pa in component delayed_rectifier_K_currents_pa_gate (second).
  * ALGEBRAIC[23] is pa_factor in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * ALGEBRAIC[11] is p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * CONSTANTS_M[17] is g_B_Na in component background_currents (nanoS).
  * CONSTANTS_M[18] is g_B_Ca in component background_currents (nanoS).
  * ALGEBRAIC[51] is E_Ca in component background_currents (millivolt).
  * STATES[18] is Ca_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[19] is Ca_i in component intracellular_ion_concentrations (millimolar).
  * CONSTANTS_M[19] is K_NaK_K in component sodium_potassium_pump (millimolar).
  * CONSTANTS_M[20] is i_NaK_max in component sodium_potassium_pump (picoA).
  * CONSTANTS_M[21] is pow_K_NaK_Na_15 in component sodium_potassium_pump (millimolar15).
  * ALGEBRAIC[53] is pow_Na_i_15 in component sodium_potassium_pump (millimolar15).
  * CONSTANTS_M[22] is i_CaP_max in component sarcolemmal_calcium_pump_current (picoA).
  * CONSTANTS_M[23] is k_CaP in component sarcolemmal_calcium_pump_current (millimolar).
  * CONSTANTS_M[24] is K_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4).
  * CONSTANTS_M[25] is d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4).
  * CONSTANTS_M[26] is gamma_Na in component Na_Ca_ion_exchanger_current (dimensionless).
  * CONSTANTS_M[27] is ACh in component ACh_dependent_K_current (millimolar).
  * CONSTANTS_M[28] is phi_Na_en in component intracellular_ion_concentrations (picoA).
  * CONSTANTS_M[29] is Vol_i in component intracellular_ion_concentrations (nanolitre).
  * CONSTANTS_M[30] is Vol_d in component intracellular_ion_concentrations (nanolitre).
  * ALGEBRAIC[58] is i_di in component intracellular_ion_concentrations (picoA).
  * CONSTANTS_M[31] is tau_di in component intracellular_ion_concentrations (second).
  * ALGEBRAIC[67] is i_up in component Ca_handling_by_the_SR (picoA).
  * ALGEBRAIC[66] is i_rel in component Ca_handling_by_the_SR (picoA).
  * ALGEBRAIC[63] is J_O in component intracellular_Ca_buffering (per_second).
  * STATES[20] is O_C in component intracellular_Ca_buffering (dimensionless).
  * STATES[21] is O_TC in component intracellular_Ca_buffering (dimensionless).
  * STATES[22] is O_TMgC in component intracellular_Ca_buffering (dimensionless).
  * STATES[23] is O_TMgMg in component intracellular_Ca_buffering (dimensionless).
  * STATES[24] is O in component intracellular_Ca_buffering (dimensionless).
  * ALGEBRAIC[60] is J_O_C in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[61] is J_O_TC in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[62] is J_O_TMgC in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[12] is J_O_TMgMg in component intracellular_Ca_buffering (per_second).
  * CONSTANTS_M[32] is Mg_i in component intracellular_Ca_buffering (millimolar).
  * CONSTANTS_M[33] is Vol_c in component cleft_space_ion_concentrations (nanolitre).
  * CONSTANTS_M[34] is tau_Na in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[35] is tau_K in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[36] is tau_Ca in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[37] is Na_b in component cleft_space_ion_concentrations (millimolar).
  * CONSTANTS_M[38] is Ca_b in component cleft_space_ion_concentrations (millimolar).
  * CONSTANTS_M[39] is K_b in component cleft_space_ion_concentrations (millimolar).
  * ALGEBRAIC[68] is i_tr in component Ca_handling_by_the_SR (picoA).
  * CONSTANTS_M[40] is I_up_max in component Ca_handling_by_the_SR (picoA).
  * CONSTANTS_M[41] is k_cyca in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[42] is k_srca in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[43] is k_xcs in component Ca_handling_by_the_SR (dimensionless).
  * CONSTANTS_M[44] is alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar).
  * STATES[25] is Ca_rel in component Ca_handling_by_the_SR (millimolar).
  * STATES[26] is Ca_up in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[45] is Vol_up in component Ca_handling_by_the_SR (nanolitre).
  * CONSTANTS_M[46] is Vol_rel in component Ca_handling_by_the_SR (nanolitre).
  * ALGEBRAIC[40] is r_act in component Ca_handling_by_the_SR (per_second).
  * ALGEBRAIC[42] is r_inact in component Ca_handling_by_the_SR (per_second).
  * CONSTANTS_M[47] is r_recov in component Ca_handling_by_the_SR (per_second).
  * ALGEBRAIC[13] is r_Ca_d_term in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[25] is r_Ca_i_term in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[36] is r_Ca_d_factor in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[38] is r_Ca_i_factor in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[64] is i_rel_f2 in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[65] is i_rel_factor in component Ca_handling_by_the_SR (dimensionless).
  * STATES[27] is O_Calse in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[69] is J_O_Calse in component Ca_handling_by_the_SR (per_second).
  * STATES[28] is F1 in component Ca_handling_by_the_SR (dimensionless).
  * STATES[29] is F2 in component Ca_handling_by_the_SR (dimensionless).
  * CONSTANTS_M[48] is tau_tr in component Ca_handling_by_the_SR (second).
  * CONSTANTS_M[49] is k_rel_i in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[50] is k_rel_d in component Ca_handling_by_the_SR (millimolar).
  * RATES[0] is d/dt V in component membrane (millivolt).
  * RATES[3] is d/dt m in component sodium_current_m_gate (dimensionless).
  * RATES[4] is d/dt h1 in component sodium_current_h1_gate (dimensionless).
  * RATES[5] is d/dt h2 in component sodium_current_h2_gate (dimensionless).
  * RATES[7] is d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
  * RATES[8] is d/dt f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * RATES[9] is d/dt f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
  * RATES[12] is d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * RATES[13] is d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * RATES[14] is d/dt a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
  * RATES[15] is d/dt i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
  * RATES[16] is d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * RATES[17] is d/dt pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * RATES[11] is d/dt K_i in component intracellular_ion_concentrations (millimolar).
  * RATES[2] is d/dt Na_i in component intracellular_ion_concentrations (millimolar).
  * RATES[19] is d/dt Ca_i in component intracellular_ion_concentrations (millimolar).
  * RATES[6] is d/dt Ca_d in component intracellular_ion_concentrations (millimolar).
  * RATES[20] is d/dt O_C in component intracellular_Ca_buffering (dimensionless).
  * RATES[21] is d/dt O_TC in component intracellular_Ca_buffering (dimensionless).
  * RATES[22] is d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless).
  * RATES[23] is d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless).
  * RATES[24] is d/dt O in component intracellular_Ca_buffering (dimensionless).
  * RATES[18] is d/dt Ca_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[10] is d/dt K_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[1] is d/dt Na_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[28] is d/dt F1 in component Ca_handling_by_the_SR (dimensionless).
  * RATES[29] is d/dt F2 in component Ca_handling_by_the_SR (dimensionless).
  * RATES[27] is d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless).
  * RATES[26] is d/dt Ca_up in component Ca_handling_by_the_SR (millimolar).
  * RATES[25] is d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar).
  */

  PDEFIELD_TYPE CONSTANTS_M[116];
  PDEFIELD_TYPE ALGEBRAIC[101];
  //PDEFIELD_TYPE scalarZyantkekorov = 81/50;

  CONSTANTS_M[0] = 8314;
  CONSTANTS_M[1] = 306.15;
  CONSTANTS_M[2] = 96487;
  CONSTANTS_M[3] = 50;
  CONSTANTS_M[4] = 0;
  CONSTANTS_M[5] = 5;
  CONSTANTS_M[6] = 0.006;
  CONSTANTS_M[7] = activation_strength;//-15;
  if (id/410 >= 600000){
    //printf("id = %i, so x = %i", id, id%410);
    CONSTANTS_M[7] = -2000;//-15;
  }
  CONSTANTS_M[8] = 0.0018;
  CONSTANTS_M[9] = 6.75;
  CONSTANTS_M[10] = 60;
  CONSTANTS_M[11] = 0.025;
  CONSTANTS_M[12] = 8.25;
  CONSTANTS_M[13] = 2.25;
  CONSTANTS_M[14] = 3.1;
  CONSTANTS_M[15] = 1;
  CONSTANTS_M[16] = 0.5;
  CONSTANTS_M[17] = 0.060599;
  CONSTANTS_M[18] = 0.078681;
  CONSTANTS_M[19] = 1;
  CONSTANTS_M[20] = 68.55;
  CONSTANTS_M[21] = 36.4829;
  CONSTANTS_M[22] = 4;
  CONSTANTS_M[23] = 0.0002;
  CONSTANTS_M[24] = 0.0374842;
  CONSTANTS_M[25] = 0.0003;
  CONSTANTS_M[26] = 0.45;
  CONSTANTS_M[27] = 1e-24;
  CONSTANTS_M[28] = 0;
  CONSTANTS_M[29] = 0.005884;
  CONSTANTS_M[30] = 0.00011768;
  CONSTANTS_M[31] = 0.01;
  CONSTANTS_M[32] = 2.5;
  CONSTANTS_M[33] = 0.000800224;
  CONSTANTS_M[34] = 14.3;
  CONSTANTS_M[35] = 10;
  CONSTANTS_M[36] = 24.7;
  CONSTANTS_M[37] = 130;
  CONSTANTS_M[38] = 1.8;
  CONSTANTS_M[39] = 5.4;
  CONSTANTS_M[40] = 2800;
  CONSTANTS_M[41] = 0.0003;
  CONSTANTS_M[42] = 0.5;
  CONSTANTS_M[43] = 0.4;
  CONSTANTS_M[44] = 200000;
  CONSTANTS_M[45] = 0.0003969;
  CONSTANTS_M[46] = 0.0000441;
  CONSTANTS_M[47] = 0.815;
  CONSTANTS_M[48] = 0.01;
  CONSTANTS_M[49] = 0.0003;
  CONSTANTS_M[50] = 0.003;

  ALGEBRAIC[12] =  2000.00*CONSTANTS_M[32]*((1.00000 - STATES[22]) - STATES[23]) -  666.000*STATES[23];
  RATES[23] = ALGEBRAIC[12];
  ALGEBRAIC[18] =  0.00350000*exp((( - STATES[0]*STATES[0])/30.0000)/30.0000)+0.00150000;
  ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0] - 1.00000)/- 11.0000));
  RATES[12] = (ALGEBRAIC[6] - STATES[12])/ALGEBRAIC[18];
  ALGEBRAIC[8] = 1.00000/(1.00000+exp(- (STATES[0]+6.00000)/8.60000));
  ALGEBRAIC[20] = 0.00900000/(1.00000+exp((STATES[0]+5.00000)/12.0000))+0.000500000;
  RATES[14] = (ALGEBRAIC[8] - STATES[14])/ALGEBRAIC[20];
  ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+7.50000)/10.0000));
  ALGEBRAIC[21] = 0.590000/(1.00000+exp((STATES[0]+60.0000)/10.0000))+3.05000;
  RATES[15] = (ALGEBRAIC[9] - STATES[15])/ALGEBRAIC[21];
  ALGEBRAIC[14] = 1.00000/(1.00000+exp((STATES[0]+27.1200)/- 8.21000));
  ALGEBRAIC[2] = (STATES[0]+25.5700)/28.8000;
  ALGEBRAIC[26] =  4.20000e-05*exp( - ALGEBRAIC[2]*ALGEBRAIC[2])+2.40000e-05;
  RATES[3] = (ALGEBRAIC[14] - STATES[3])/ALGEBRAIC[26];
  ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+63.6000)/5.30000));
  ALGEBRAIC[15] = 1.00000/(1.00000+exp((STATES[0]+35.1000)/3.20000));
  ALGEBRAIC[27] =  0.0300000*ALGEBRAIC[15]+0.000300000;
  RATES[4] = (ALGEBRAIC[3] - STATES[4])/ALGEBRAIC[27];
  ALGEBRAIC[28] =  0.120000*ALGEBRAIC[15]+0.00300000;
  RATES[5] = (ALGEBRAIC[3] - STATES[5])/ALGEBRAIC[28];
  ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+9.00000)/- 5.80000));
  ALGEBRAIC[16] = (STATES[0]+35.0000)/30.0000;
  ALGEBRAIC[29] =  0.00270000*exp( - ALGEBRAIC[16]*ALGEBRAIC[16])+0.00200000;
  RATES[7] = (ALGEBRAIC[4] - STATES[7])/ALGEBRAIC[29];
  ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0]+27.4000)/7.10000));
  ALGEBRAIC[17] = STATES[0]+40.0000;
  ALGEBRAIC[30] =  0.161000*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.4000)/14.4000)+0.0100000;
  RATES[8] = (ALGEBRAIC[5] - STATES[8])/ALGEBRAIC[30];
  ALGEBRAIC[31] =  1.33230*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.2000)/14.2000)+0.0626000;
  RATES[9] = (ALGEBRAIC[5] - STATES[9])/ALGEBRAIC[31];
  ALGEBRAIC[19] = (STATES[0]+52.4500)/15.8827;
  ALGEBRAIC[32] =  0.0256350*exp( - ALGEBRAIC[19]*ALGEBRAIC[19])+0.0141400;
  ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[0]+40.5000)/11.5000));
  RATES[13] = (ALGEBRAIC[7] - STATES[13])/ALGEBRAIC[32];
  ALGEBRAIC[22] = (STATES[0] - 20.0000)/20.0000;
  ALGEBRAIC[33] = 0.700000+ 0.400000*exp( - ALGEBRAIC[22]*ALGEBRAIC[22]);
  ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0] - 19.9000)/- 12.7000));
  RATES[16] = (ALGEBRAIC[10] - STATES[16])/ALGEBRAIC[33];
  ALGEBRAIC[23] = (STATES[0]+20.1376)/22.1996;
  ALGEBRAIC[34] = 0.0311800+ 0.217180*exp( - ALGEBRAIC[23]*ALGEBRAIC[23]);
  ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+15.0000)/- 6.00000));
  RATES[17] = (ALGEBRAIC[11] - STATES[17])/ALGEBRAIC[34];
  ALGEBRAIC[13] = STATES[6]/(STATES[6]+CONSTANTS_M[50]);
  ALGEBRAIC[36] =  ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13];
  ALGEBRAIC[25] = STATES[19]/(STATES[19]+CONSTANTS_M[49]);
  ALGEBRAIC[38] =  ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25];
  ALGEBRAIC[40] =  203.800*(ALGEBRAIC[38]+ALGEBRAIC[36]);
  RATES[28] =  CONSTANTS_M[47]*((1.00000 - STATES[28]) - STATES[29]) -  ALGEBRAIC[40]*STATES[28];
  ALGEBRAIC[42] = 33.9600+ 339.600*ALGEBRAIC[38];
  RATES[29] =  ALGEBRAIC[40]*STATES[28] -  ALGEBRAIC[42]*STATES[29];
  ALGEBRAIC[43] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/CONSTANTS_M[2])*log(STATES[10]/STATES[11]);
  ALGEBRAIC[44] =  CONSTANTS_M[12]*STATES[12]*STATES[13]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[45] =  CONSTANTS_M[13]*STATES[14]*STATES[15]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[46] = ( CONSTANTS_M[14]*pow(STATES[10]/1.00000, 0.445700)*(STATES[0] - ALGEBRAIC[43]))/(1.00000+exp(( 1.50000*((STATES[0] - ALGEBRAIC[43])+3.60000)*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])));
  ALGEBRAIC[48] = 1.00000/(1.00000+exp((STATES[0]+55.0000)/24.0000));
  ALGEBRAIC[49] =  CONSTANTS_M[16]*STATES[17]*ALGEBRAIC[48]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[47] =  CONSTANTS_M[15]*STATES[16]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[53] = pow(STATES[2], 1.50000);
  ALGEBRAIC[54] = ( (( (( CONSTANTS_M[20]*STATES[10])/(STATES[10]+CONSTANTS_M[19]))*ALGEBRAIC[53])/(ALGEBRAIC[53]+CONSTANTS_M[21]))*(STATES[0]+150.000))/(STATES[0]+200.000);
  ALGEBRAIC[1] =  floor(VOI/CONSTANTS_M[5])*CONSTANTS_M[5];
  ALGEBRAIC[24] = (VOI - ALGEBRAIC[1]>=CONSTANTS_M[4]&&VOI - ALGEBRAIC[1]<=CONSTANTS_M[4]+CONSTANTS_M[6] ? CONSTANTS_M[7] : 0.00000);
  RATES[11] = - (((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])+ ALGEBRAIC[24]*CONSTANTS_M[3])/( CONSTANTS_M[29]*CONSTANTS_M[2]);
  RATES[10] = (CONSTANTS_M[39] - STATES[10])/CONSTANTS_M[35]+((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])/( CONSTANTS_M[33]*CONSTANTS_M[2]);
  ALGEBRAIC[35] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/CONSTANTS_M[2])*log(STATES[1]/STATES[2]);
  ALGEBRAIC[37] = ( (( CONSTANTS_M[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*STATES[0]*CONSTANTS_M[2]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))*(exp(( (STATES[0] - ALGEBRAIC[35])*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000))/(exp(( STATES[0]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000);
  if(!isfinite(ALGEBRAIC[37]))
    ALGEBRAIC[37] = ( (( CONSTANTS_M[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*CONSTANTS_M[2]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))*(exp((-ALGEBRAIC[35]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000))/(CONSTANTS_M[2]/(CONSTANTS_M[0]*CONSTANTS_M[1]) - 1.00000);
  ALGEBRAIC[37] = ALGEBRAIC[37]*I_Na_factor;
  ALGEBRAIC[50] =  CONSTANTS_M[17]*(STATES[0] - ALGEBRAIC[35]);
  ALGEBRAIC[56] = ( CONSTANTS_M[24]*( STATES[2]*STATES[2]*STATES[2]*STATES[18]*exp(( CONSTANTS_M[2]*STATES[0]*CONSTANTS_M[26])/( CONSTANTS_M[0]*CONSTANTS_M[1])) -  STATES[1]*STATES[1]*STATES[1]*STATES[19]*exp(( (CONSTANTS_M[26] - 1.00000)*STATES[0]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))))/(1.00000+ CONSTANTS_M[25]*( STATES[1]*STATES[1]*STATES[1]*STATES[19]+ STATES[2]*STATES[2]*STATES[2]*STATES[18]));
  ALGEBRAIC[39] = STATES[6]/(STATES[6]+CONSTANTS_M[11]);
  ALGEBRAIC[41] =  CONSTANTS_M[9]*STATES[7]*( ALGEBRAIC[39]*STATES[8]+ (1.00000 - ALGEBRAIC[39])*STATES[9])*(STATES[0] - CONSTANTS_M[10]);
  ALGEBRAIC[51] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/( 2.00000*CONSTANTS_M[2]))*log(STATES[18]/STATES[19]);
  ALGEBRAIC[52] =  CONSTANTS_M[18]*(STATES[0] - ALGEBRAIC[51]);
  ALGEBRAIC[55] = ( CONSTANTS_M[22]*STATES[19])/(STATES[19]+CONSTANTS_M[23]);
  RATES[18] = (CONSTANTS_M[38] - STATES[18])/CONSTANTS_M[36]+((ALGEBRAIC[41]+ALGEBRAIC[52]+ALGEBRAIC[55]) -  2.00000*ALGEBRAIC[56])/( 2.00000*CONSTANTS_M[33]*CONSTANTS_M[2]);
  RATES[1] = (CONSTANTS_M[37] - STATES[1])/CONSTANTS_M[34]+(ALGEBRAIC[37]+ALGEBRAIC[50]+ 3.00000*ALGEBRAIC[56]+ 3.00000*ALGEBRAIC[54]+CONSTANTS_M[28])/( CONSTANTS_M[33]*CONSTANTS_M[2]);
  ALGEBRAIC[58] = ( (STATES[6] - STATES[19])*2.00000*CONSTANTS_M[30]*CONSTANTS_M[2])/CONSTANTS_M[31];
  RATES[6] = - (ALGEBRAIC[41]+ALGEBRAIC[58])/( 2.00000*CONSTANTS_M[30]*CONSTANTS_M[2]);
  ALGEBRAIC[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(CONSTANTS_M[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((STATES[0]+59.5300)/17.1800)))*(STATES[0] - ALGEBRAIC[43])*CONSTANTS_M[3];
  ALGEBRAIC[59] = (ALGEBRAIC[37]+ALGEBRAIC[41]+ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[49]+ALGEBRAIC[47]+ALGEBRAIC[50]+ALGEBRAIC[52]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[57])/CONSTANTS_M[3]+ALGEBRAIC[24];
  RATES[0] =  - ALGEBRAIC[59]*1000.00;
  ALGEBRAIC[60] =  200000.*STATES[19]*(1.00000 - STATES[20]) -  476.000*STATES[20];
  RATES[20] = ALGEBRAIC[60];
  ALGEBRAIC[61] =  78400.0*STATES[19]*(1.00000 - STATES[21]) -  392.000*STATES[21];
  RATES[21] = ALGEBRAIC[61];
  ALGEBRAIC[62] =  200000.*STATES[19]*((1.00000 - STATES[22]) - STATES[23]) -  6.60000*STATES[22];
  RATES[22] = ALGEBRAIC[62];
  ALGEBRAIC[63] =  0.0800000*ALGEBRAIC[61]+ 0.160000*ALGEBRAIC[62]+ 0.0450000*ALGEBRAIC[60];
  RATES[24] = ALGEBRAIC[63];
  ALGEBRAIC[67] = ( CONSTANTS_M[40]*(STATES[19]/CONSTANTS_M[41] - ( CONSTANTS_M[43]*CONSTANTS_M[43]*STATES[26])/CONSTANTS_M[42]))/((STATES[19]+CONSTANTS_M[41])/CONSTANTS_M[41]+( CONSTANTS_M[43]*(STATES[26]+CONSTANTS_M[42]))/CONSTANTS_M[42]);
  ALGEBRAIC[64] = STATES[29]/(STATES[29]+0.250000);
  ALGEBRAIC[65] =  ALGEBRAIC[64]*ALGEBRAIC[64];
  ALGEBRAIC[66] =  CONSTANTS_M[44]*ALGEBRAIC[65]*(STATES[25] - STATES[19]);
  RATES[19] = - ((ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[67]) - (ALGEBRAIC[58]+ALGEBRAIC[66]+ 2.00000*ALGEBRAIC[56]))/( 2.00000*CONSTANTS_M[29]*CONSTANTS_M[2]) -  1.00000*ALGEBRAIC[63];
  ALGEBRAIC[68] = ( (STATES[26] - STATES[25])*2.00000*CONSTANTS_M[46]*CONSTANTS_M[2])/CONSTANTS_M[48];
  RATES[26] = (ALGEBRAIC[67] - ALGEBRAIC[68])/( 2.00000*CONSTANTS_M[45]*CONSTANTS_M[2]);
  ALGEBRAIC[69] =  480.000*STATES[25]*(1.00000 - STATES[27]) -  400.000*STATES[27];
  RATES[27] = ALGEBRAIC[69];
  RATES[25] = (ALGEBRAIC[68] - ALGEBRAIC[66])/( 2.00000*CONSTANTS_M[46]*CONSTANTS_M[2]) -  31.0000*ALGEBRAIC[69];

}

void derivsMaleckar_host(PDEFIELD_TYPE VOI, PDEFIELD_TYPE* STATES,PDEFIELD_TYPE* RATES, PDEFIELD_TYPE pacing_interval, int id, PDEFIELD_TYPE activation_strength, PDEFIELD_TYPE t_A){


  /*
    There are a total of 70 entries in the algebraic variable array.
    There are a total of 30 entries in each of the rate and state variable arrays.
    There are a total of 51 entries in the constant variable array.
  */
  /*
  * VOI is time in component environment (second).
  * STATES[0] is V in component membrane (millivolt).
  * CONSTANTS_M[0] is R in component membrane (millijoule_per_mole_kelvin).
  * CONSTANTS_M[1] is T in component membrane (kelvin).
  * CONSTANTS_M[2] is F in component membrane (coulomb_per_mole).
  * CONSTANTS_M[3] is Cm in component membrane (nanoF).
  * ALGEBRAIC[0] is Q_tot in component membrane (millivolt).
  * ALGEBRAIC[37] is i_Na in component sodium_current (picoA).
  * ALGEBRAIC[41] is i_Ca_L in component L_type_Ca_channel (picoA).
  * ALGEBRAIC[44] is i_t in component Ca_independent_transient_outward_K_current (picoA).
  * ALGEBRAIC[45] is i_Kur in component ultra_rapid_K_current (picoA).
  * ALGEBRAIC[46] is i_K1 in component inward_rectifier (picoA).
  * ALGEBRAIC[49] is i_Kr in component delayed_rectifier_K_currents (picoA).
  * ALGEBRAIC[47] is i_Ks in component delayed_rectifier_K_currents (picoA).
  * ALGEBRAIC[50] is i_B_Na in component background_currents (picoA).
  * ALGEBRAIC[52] is i_B_Ca in component background_currents (picoA).
  * ALGEBRAIC[54] is i_NaK in component sodium_potassium_pump (picoA).
  * ALGEBRAIC[55] is i_CaP in component sarcolemmal_calcium_pump_current (picoA).
  * ALGEBRAIC[56] is i_NaCa in component Na_Ca_ion_exchanger_current (picoA).
  * ALGEBRAIC[57] is i_KACh in component ACh_dependent_K_current (picoA).
  * ALGEBRAIC[59] is I in component membrane (pA_per_nF).
  * ALGEBRAIC[24] is i_Stim in component membrane (pA_per_nF).
  * CONSTANTS_M[4] is stim_offset in component membrane (second).
  * CONSTANTS_M[5] is stim_period in component membrane (second).
  * CONSTANTS_M[6] is stim_duration in component membrane (second).
  * CONSTANTS_M[7] is stim_amplitude in component membrane (pA_per_nF).
  * ALGEBRAIC[1] is past in component membrane (second).
  * ALGEBRAIC[35] is E_Na in component sodium_current (millivolt).
  * CONSTANTS_M[8] is P_Na in component sodium_current (nanolitre_per_second).
  * STATES[1] is Na_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[2] is Na_i in component intracellular_ion_concentrations (millimolar).
  * STATES[3] is m in component sodium_current_m_gate (dimensionless).
  * STATES[4] is h1 in component sodium_current_h1_gate (dimensionless).
  * STATES[5] is h2 in component sodium_current_h2_gate (dimensionless).
  * ALGEBRAIC[14] is m_infinity in component sodium_current_m_gate (dimensionless).
  * ALGEBRAIC[2] is m_factor in component sodium_current_m_gate (dimensionless).
  * ALGEBRAIC[26] is tau_m in component sodium_current_m_gate (second).
  * ALGEBRAIC[3] is h_infinity in component sodium_current_h1_gate (dimensionless).
  * ALGEBRAIC[15] is h_factor in component sodium_current_h1_gate (dimensionless).
  * ALGEBRAIC[27] is tau_h1 in component sodium_current_h1_gate (second).
  * ALGEBRAIC[28] is tau_h2 in component sodium_current_h2_gate (second).
  * CONSTANTS_M[9] is g_Ca_L in component L_type_Ca_channel (nanoS).
  * CONSTANTS_M[10] is E_Ca_app in component L_type_Ca_channel (millivolt).
  * ALGEBRAIC[39] is f_Ca in component L_type_Ca_channel (dimensionless).
  * CONSTANTS_M[11] is k_Ca in component L_type_Ca_channel (millimolar).
  * STATES[6] is Ca_d in component intracellular_ion_concentrations (millimolar).
  * STATES[7] is d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
  * STATES[8] is f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * STATES[9] is f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
  * ALGEBRAIC[4] is d_L_infinity in component L_type_Ca_channel_d_L_gate (dimensionless).
  * ALGEBRAIC[16] is d_L_factor in component L_type_Ca_channel_d_L_gate (dimensionless).
  * ALGEBRAIC[29] is tau_d_L in component L_type_Ca_channel_d_L_gate (second).
  * ALGEBRAIC[5] is f_L_infinity in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * ALGEBRAIC[17] is f_L_factor in component L_type_Ca_channel_f_L1_gate (millivolt).
  * ALGEBRAIC[30] is tau_f_L1 in component L_type_Ca_channel_f_L1_gate (second).
  * ALGEBRAIC[31] is tau_f_L2 in component L_type_Ca_channel_f_L2_gate (second).
  * ALGEBRAIC[43] is E_K in component Ca_independent_transient_outward_K_current (millivolt).
  * CONSTANTS_M[12] is g_t in component Ca_independent_transient_outward_K_current (nanoS).
  * STATES[10] is K_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[11] is K_i in component intracellular_ion_concentrations (millimolar).
  * STATES[12] is r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * STATES[13] is s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * ALGEBRAIC[18] is tau_r in component Ca_independent_transient_outward_K_current_r_gate (second).
  * ALGEBRAIC[6] is r_infinity in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * ALGEBRAIC[32] is tau_s in component Ca_independent_transient_outward_K_current_s_gate (second).
  * ALGEBRAIC[7] is s_infinity in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * ALGEBRAIC[19] is s_factor in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * CONSTANTS_M[13] is g_kur in component ultra_rapid_K_current (nanoS).
  * STATES[14] is a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
  * STATES[15] is i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
  * ALGEBRAIC[8] is a_ur_infinity in component ultra_rapid_K_current_aur_gate (dimensionless).
  * ALGEBRAIC[20] is tau_a_ur in component ultra_rapid_K_current_aur_gate (second).
  * ALGEBRAIC[9] is i_ur_infinity in component ultra_rapid_K_current_iur_gate (dimensionless).
  * ALGEBRAIC[21] is tau_i_ur in component ultra_rapid_K_current_iur_gate (second).
  * CONSTANTS_M[14] is g_K1 in component inward_rectifier (nanoS).
  * CONSTANTS_M[15] is g_Ks in component delayed_rectifier_K_currents (nanoS).
  * CONSTANTS_M[16] is g_Kr in component delayed_rectifier_K_currents (nanoS).
  * STATES[16] is n in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * STATES[17] is pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * ALGEBRAIC[48] is pip in component delayed_rectifier_K_currents_pi_gate (dimensionless).
  * ALGEBRAIC[33] is tau_n in component delayed_rectifier_K_currents_n_gate (second).
  * ALGEBRAIC[10] is n_infinity in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * ALGEBRAIC[22] is n_factor in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * ALGEBRAIC[34] is tau_pa in component delayed_rectifier_K_currents_pa_gate (second).
  * ALGEBRAIC[23] is pa_factor in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * ALGEBRAIC[11] is p_a_infinity in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * CONSTANTS_M[17] is g_B_Na in component background_currents (nanoS).
  * CONSTANTS_M[18] is g_B_Ca in component background_currents (nanoS).
  * ALGEBRAIC[51] is E_Ca in component background_currents (millivolt).
  * STATES[18] is Ca_c in component cleft_space_ion_concentrations (millimolar).
  * STATES[19] is Ca_i in component intracellular_ion_concentrations (millimolar).
  * CONSTANTS_M[19] is K_NaK_K in component sodium_potassium_pump (millimolar).
  * CONSTANTS_M[20] is i_NaK_max in component sodium_potassium_pump (picoA).
  * CONSTANTS_M[21] is pow_K_NaK_Na_15 in component sodium_potassium_pump (millimolar15).
  * ALGEBRAIC[53] is pow_Na_i_15 in component sodium_potassium_pump (millimolar15).
  * CONSTANTS_M[22] is i_CaP_max in component sarcolemmal_calcium_pump_current (picoA).
  * CONSTANTS_M[23] is k_CaP in component sarcolemmal_calcium_pump_current (millimolar).
  * CONSTANTS_M[24] is K_NaCa in component Na_Ca_ion_exchanger_current (picoA_per_millimolar_4).
  * CONSTANTS_M[25] is d_NaCa in component Na_Ca_ion_exchanger_current (per_millimolar_4).
  * CONSTANTS_M[26] is gamma_Na in component Na_Ca_ion_exchanger_current (dimensionless).
  * CONSTANTS_M[27] is ACh in component ACh_dependent_K_current (millimolar).
  * CONSTANTS_M[28] is phi_Na_en in component intracellular_ion_concentrations (picoA).
  * CONSTANTS_M[29] is Vol_i in component intracellular_ion_concentrations (nanolitre).
  * CONSTANTS_M[30] is Vol_d in component intracellular_ion_concentrations (nanolitre).
  * ALGEBRAIC[58] is i_di in component intracellular_ion_concentrations (picoA).
  * CONSTANTS_M[31] is tau_di in component intracellular_ion_concentrations (second).
  * ALGEBRAIC[67] is i_up in component Ca_handling_by_the_SR (picoA).
  * ALGEBRAIC[66] is i_rel in component Ca_handling_by_the_SR (picoA).
  * ALGEBRAIC[63] is J_O in component intracellular_Ca_buffering (per_second).
  * STATES[20] is O_C in component intracellular_Ca_buffering (dimensionless).
  * STATES[21] is O_TC in component intracellular_Ca_buffering (dimensionless).
  * STATES[22] is O_TMgC in component intracellular_Ca_buffering (dimensionless).
  * STATES[23] is O_TMgMg in component intracellular_Ca_buffering (dimensionless).
  * STATES[24] is O in component intracellular_Ca_buffering (dimensionless).
  * ALGEBRAIC[60] is J_O_C in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[61] is J_O_TC in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[62] is J_O_TMgC in component intracellular_Ca_buffering (per_second).
  * ALGEBRAIC[12] is J_O_TMgMg in component intracellular_Ca_buffering (per_second).
  * CONSTANTS_M[32] is Mg_i in component intracellular_Ca_buffering (millimolar).
  * CONSTANTS_M[33] is Vol_c in component cleft_space_ion_concentrations (nanolitre).
  * CONSTANTS_M[34] is tau_Na in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[35] is tau_K in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[36] is tau_Ca in component cleft_space_ion_concentrations (second).
  * CONSTANTS_M[37] is Na_b in component cleft_space_ion_concentrations (millimolar).
  * CONSTANTS_M[38] is Ca_b in component cleft_space_ion_concentrations (millimolar).
  * CONSTANTS_M[39] is K_b in component cleft_space_ion_concentrations (millimolar).
  * ALGEBRAIC[68] is i_tr in component Ca_handling_by_the_SR (picoA).
  * CONSTANTS_M[40] is I_up_max in component Ca_handling_by_the_SR (picoA).
  * CONSTANTS_M[41] is k_cyca in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[42] is k_srca in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[43] is k_xcs in component Ca_handling_by_the_SR (dimensionless).
  * CONSTANTS_M[44] is alpha_rel in component Ca_handling_by_the_SR (picoA_per_millimolar).
  * STATES[25] is Ca_rel in component Ca_handling_by_the_SR (millimolar).
  * STATES[26] is Ca_up in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[45] is Vol_up in component Ca_handling_by_the_SR (nanolitre).
  * CONSTANTS_M[46] is Vol_rel in component Ca_handling_by_the_SR (nanolitre).
  * ALGEBRAIC[40] is r_act in component Ca_handling_by_the_SR (per_second).
  * ALGEBRAIC[42] is r_inact in component Ca_handling_by_the_SR (per_second).
  * CONSTANTS_M[47] is r_recov in component Ca_handling_by_the_SR (per_second).
  * ALGEBRAIC[13] is r_Ca_d_term in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[25] is r_Ca_i_term in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[36] is r_Ca_d_factor in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[38] is r_Ca_i_factor in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[64] is i_rel_f2 in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[65] is i_rel_factor in component Ca_handling_by_the_SR (dimensionless).
  * STATES[27] is O_Calse in component Ca_handling_by_the_SR (dimensionless).
  * ALGEBRAIC[69] is J_O_Calse in component Ca_handling_by_the_SR (per_second).
  * STATES[28] is F1 in component Ca_handling_by_the_SR (dimensionless).
  * STATES[29] is F2 in component Ca_handling_by_the_SR (dimensionless).
  * CONSTANTS_M[48] is tau_tr in component Ca_handling_by_the_SR (second).
  * CONSTANTS_M[49] is k_rel_i in component Ca_handling_by_the_SR (millimolar).
  * CONSTANTS_M[50] is k_rel_d in component Ca_handling_by_the_SR (millimolar).
  * RATES[0] is d/dt V in component membrane (millivolt).
  * RATES[3] is d/dt m in component sodium_current_m_gate (dimensionless).
  * RATES[4] is d/dt h1 in component sodium_current_h1_gate (dimensionless).
  * RATES[5] is d/dt h2 in component sodium_current_h2_gate (dimensionless).
  * RATES[7] is d/dt d_L in component L_type_Ca_channel_d_L_gate (dimensionless).
  * RATES[8] is d/dt f_L1 in component L_type_Ca_channel_f_L1_gate (dimensionless).
  * RATES[9] is d/dt f_L2 in component L_type_Ca_channel_f_L2_gate (dimensionless).
  * RATES[12] is d/dt r in component Ca_independent_transient_outward_K_current_r_gate (dimensionless).
  * RATES[13] is d/dt s in component Ca_independent_transient_outward_K_current_s_gate (dimensionless).
  * RATES[14] is d/dt a_ur in component ultra_rapid_K_current_aur_gate (dimensionless).
  * RATES[15] is d/dt i_ur in component ultra_rapid_K_current_iur_gate (dimensionless).
  * RATES[16] is d/dt n in component delayed_rectifier_K_currents_n_gate (dimensionless).
  * RATES[17] is d/dt pa in component delayed_rectifier_K_currents_pa_gate (dimensionless).
  * RATES[11] is d/dt K_i in component intracellular_ion_concentrations (millimolar).
  * RATES[2] is d/dt Na_i in component intracellular_ion_concentrations (millimolar).
  * RATES[19] is d/dt Ca_i in component intracellular_ion_concentrations (millimolar).
  * RATES[6] is d/dt Ca_d in component intracellular_ion_concentrations (millimolar).
  * RATES[20] is d/dt O_C in component intracellular_Ca_buffering (dimensionless).
  * RATES[21] is d/dt O_TC in component intracellular_Ca_buffering (dimensionless).
  * RATES[22] is d/dt O_TMgC in component intracellular_Ca_buffering (dimensionless).
  * RATES[23] is d/dt O_TMgMg in component intracellular_Ca_buffering (dimensionless).
  * RATES[24] is d/dt O in component intracellular_Ca_buffering (dimensionless).
  * RATES[18] is d/dt Ca_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[10] is d/dt K_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[1] is d/dt Na_c in component cleft_space_ion_concentrations (millimolar).
  * RATES[28] is d/dt F1 in component Ca_handling_by_the_SR (dimensionless).
  * RATES[29] is d/dt F2 in component Ca_handling_by_the_SR (dimensionless).
  * RATES[27] is d/dt O_Calse in component Ca_handling_by_the_SR (dimensionless).
  * RATES[26] is d/dt Ca_up in component Ca_handling_by_the_SR (millimolar).
  * RATES[25] is d/dt Ca_rel in component Ca_handling_by_the_SR (millimolar).
  */

  PDEFIELD_TYPE CONSTANTS_M[116];
  PDEFIELD_TYPE ALGEBRAIC[101];
  //PDEFIELD_TYPE scalarZyantkekorov = 81/50;

  CONSTANTS_M[0] = 8314;
  CONSTANTS_M[1] = 306.15;
  CONSTANTS_M[2] = 96487;
  CONSTANTS_M[3] = 50;
  CONSTANTS_M[4] = 0;
  CONSTANTS_M[5] = 5;
  CONSTANTS_M[6] = t_A;
  CONSTANTS_M[7] = activation_strength;//-15;
  /*if (id/410 >= 600000){
    //printf("id = %i, so x = %i", id, id%410);
    CONSTANTS_M[7] = -2000;//-15;
  }*/
  CONSTANTS_M[8] = 0.0018;
  CONSTANTS_M[9] = 6.75;
  CONSTANTS_M[10] = 60;
  CONSTANTS_M[11] = 0.025;
  CONSTANTS_M[12] = 8.25;
  CONSTANTS_M[13] = 2.25;
  CONSTANTS_M[14] = 3.1;
  CONSTANTS_M[15] = 1;
  CONSTANTS_M[16] = 0.5;
  CONSTANTS_M[17] = 0.060599;
  CONSTANTS_M[18] = 0.078681;
  CONSTANTS_M[19] = 1;
  CONSTANTS_M[20] = 68.55;
  CONSTANTS_M[21] = 36.4829;
  CONSTANTS_M[22] = 4;
  CONSTANTS_M[23] = 0.0002;
  CONSTANTS_M[24] = 0.0374842;
  CONSTANTS_M[25] = 0.0003;
  CONSTANTS_M[26] = 0.45;
  CONSTANTS_M[27] = 1e-24;
  CONSTANTS_M[28] = 0;
  CONSTANTS_M[29] = 0.005884;
  CONSTANTS_M[30] = 0.00011768;
  CONSTANTS_M[31] = 0.01;
  CONSTANTS_M[32] = 2.5;
  CONSTANTS_M[33] = 0.000800224;
  CONSTANTS_M[34] = 14.3;
  CONSTANTS_M[35] = 10;
  CONSTANTS_M[36] = 24.7;
  CONSTANTS_M[37] = 130;
  CONSTANTS_M[38] = 1.8;
  CONSTANTS_M[39] = 5.4;
  CONSTANTS_M[40] = 2800;
  CONSTANTS_M[41] = 0.0003;
  CONSTANTS_M[42] = 0.5;
  CONSTANTS_M[43] = 0.4;
  CONSTANTS_M[44] = 200000;
  CONSTANTS_M[45] = 0.0003969;
  CONSTANTS_M[46] = 0.0000441;
  CONSTANTS_M[47] = 0.815;
  CONSTANTS_M[48] = 0.01;
  CONSTANTS_M[49] = 0.0003;
  CONSTANTS_M[50] = 0.003;

  ALGEBRAIC[12] =  2000.00*CONSTANTS_M[32]*((1.00000 - STATES[22]) - STATES[23]) -  666.000*STATES[23];
  RATES[23] = ALGEBRAIC[12];
  ALGEBRAIC[18] =  0.00350000*exp((( - STATES[0]*STATES[0])/30.0000)/30.0000)+0.00150000;
  ALGEBRAIC[6] = 1.00000/(1.00000+exp((STATES[0] - 1.00000)/- 11.0000));
  RATES[12] = (ALGEBRAIC[6] - STATES[12])/ALGEBRAIC[18];
  ALGEBRAIC[8] = 1.00000/(1.00000+exp(- (STATES[0]+6.00000)/8.60000));
  ALGEBRAIC[20] = 0.00900000/(1.00000+exp((STATES[0]+5.00000)/12.0000))+0.000500000;
  RATES[14] = (ALGEBRAIC[8] - STATES[14])/ALGEBRAIC[20];
  ALGEBRAIC[9] = 1.00000/(1.00000+exp((STATES[0]+7.50000)/10.0000));
  ALGEBRAIC[21] = 0.590000/(1.00000+exp((STATES[0]+60.0000)/10.0000))+3.05000;
  RATES[15] = (ALGEBRAIC[9] - STATES[15])/ALGEBRAIC[21];
  ALGEBRAIC[14] = 1.00000/(1.00000+exp((STATES[0]+27.1200)/- 8.21000));
  ALGEBRAIC[2] = (STATES[0]+25.5700)/28.8000;
  ALGEBRAIC[26] =  4.20000e-05*exp( - ALGEBRAIC[2]*ALGEBRAIC[2])+2.40000e-05;
  RATES[3] = (ALGEBRAIC[14] - STATES[3])/ALGEBRAIC[26];
  ALGEBRAIC[3] = 1.00000/(1.00000+exp((STATES[0]+63.6000)/5.30000));
  ALGEBRAIC[15] = 1.00000/(1.00000+exp((STATES[0]+35.1000)/3.20000));
  ALGEBRAIC[27] =  0.0300000*ALGEBRAIC[15]+0.000300000;
  RATES[4] = (ALGEBRAIC[3] - STATES[4])/ALGEBRAIC[27];
  ALGEBRAIC[28] =  0.120000*ALGEBRAIC[15]+0.00300000;
  RATES[5] = (ALGEBRAIC[3] - STATES[5])/ALGEBRAIC[28];
  ALGEBRAIC[4] = 1.00000/(1.00000+exp((STATES[0]+9.00000)/- 5.80000));
  ALGEBRAIC[16] = (STATES[0]+35.0000)/30.0000;
  ALGEBRAIC[29] =  0.00270000*exp( - ALGEBRAIC[16]*ALGEBRAIC[16])+0.00200000;
  RATES[7] = (ALGEBRAIC[4] - STATES[7])/ALGEBRAIC[29];
  ALGEBRAIC[5] = 1.00000/(1.00000+exp((STATES[0]+27.4000)/7.10000));
  ALGEBRAIC[17] = STATES[0]+40.0000;
  ALGEBRAIC[30] =  0.161000*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.4000)/14.4000)+0.0100000;
  RATES[8] = (ALGEBRAIC[5] - STATES[8])/ALGEBRAIC[30];
  ALGEBRAIC[31] =  1.33230*exp((( - ALGEBRAIC[17]*ALGEBRAIC[17])/14.2000)/14.2000)+0.0626000;
  RATES[9] = (ALGEBRAIC[5] - STATES[9])/ALGEBRAIC[31];
  ALGEBRAIC[19] = (STATES[0]+52.4500)/15.8827;
  ALGEBRAIC[32] =  0.0256350*exp( - ALGEBRAIC[19]*ALGEBRAIC[19])+0.0141400;
  ALGEBRAIC[7] = 1.00000/(1.00000+exp((STATES[0]+40.5000)/11.5000));
  RATES[13] = (ALGEBRAIC[7] - STATES[13])/ALGEBRAIC[32];
  ALGEBRAIC[22] = (STATES[0] - 20.0000)/20.0000;
  ALGEBRAIC[33] = 0.700000+ 0.400000*exp( - ALGEBRAIC[22]*ALGEBRAIC[22]);
  ALGEBRAIC[10] = 1.00000/(1.00000+exp((STATES[0] - 19.9000)/- 12.7000));
  RATES[16] = (ALGEBRAIC[10] - STATES[16])/ALGEBRAIC[33];
  ALGEBRAIC[23] = (STATES[0]+20.1376)/22.1996;
  ALGEBRAIC[34] = 0.0311800+ 0.217180*exp( - ALGEBRAIC[23]*ALGEBRAIC[23]);
  ALGEBRAIC[11] = 1.00000/(1.00000+exp((STATES[0]+15.0000)/- 6.00000));
  RATES[17] = (ALGEBRAIC[11] - STATES[17])/ALGEBRAIC[34];
  ALGEBRAIC[13] = STATES[6]/(STATES[6]+CONSTANTS_M[50]);
  ALGEBRAIC[36] =  ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13]*ALGEBRAIC[13];
  ALGEBRAIC[25] = STATES[19]/(STATES[19]+CONSTANTS_M[49]);
  ALGEBRAIC[38] =  ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25]*ALGEBRAIC[25];
  ALGEBRAIC[40] =  203.800*(ALGEBRAIC[38]+ALGEBRAIC[36]);
  RATES[28] =  CONSTANTS_M[47]*((1.00000 - STATES[28]) - STATES[29]) -  ALGEBRAIC[40]*STATES[28];
  ALGEBRAIC[42] = 33.9600+ 339.600*ALGEBRAIC[38];
  RATES[29] =  ALGEBRAIC[40]*STATES[28] -  ALGEBRAIC[42]*STATES[29];
  ALGEBRAIC[43] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/CONSTANTS_M[2])*log(STATES[10]/STATES[11]);
  ALGEBRAIC[44] =  CONSTANTS_M[12]*STATES[12]*STATES[13]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[45] =  CONSTANTS_M[13]*STATES[14]*STATES[15]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[46] = ( CONSTANTS_M[14]*pow(STATES[10]/1.00000, 0.445700)*(STATES[0] - ALGEBRAIC[43]))/(1.00000+exp(( 1.50000*((STATES[0] - ALGEBRAIC[43])+3.60000)*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])));
  ALGEBRAIC[48] = 1.00000/(1.00000+exp((STATES[0]+55.0000)/24.0000));
  ALGEBRAIC[49] =  CONSTANTS_M[16]*STATES[17]*ALGEBRAIC[48]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[47] =  CONSTANTS_M[15]*STATES[16]*(STATES[0] - ALGEBRAIC[43]);
  ALGEBRAIC[53] = pow(STATES[2], 1.50000);
  ALGEBRAIC[54] = ( (( (( CONSTANTS_M[20]*STATES[10])/(STATES[10]+CONSTANTS_M[19]))*ALGEBRAIC[53])/(ALGEBRAIC[53]+CONSTANTS_M[21]))*(STATES[0]+150.000))/(STATES[0]+200.000);
  ALGEBRAIC[1] =  floor(VOI/CONSTANTS_M[5])*CONSTANTS_M[5];
  ALGEBRAIC[24] = (VOI - ALGEBRAIC[1]>=CONSTANTS_M[4]&&VOI - ALGEBRAIC[1]<=CONSTANTS_M[4]+CONSTANTS_M[6] ? CONSTANTS_M[7] : 0.00000);
  RATES[11] = - (((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])+ ALGEBRAIC[24]*CONSTANTS_M[3])/( CONSTANTS_M[29]*CONSTANTS_M[2]);
  RATES[10] = (CONSTANTS_M[39] - STATES[10])/CONSTANTS_M[35]+((ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[47]+ALGEBRAIC[49]) -  2.00000*ALGEBRAIC[54])/( CONSTANTS_M[33]*CONSTANTS_M[2]);
  ALGEBRAIC[35] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/CONSTANTS_M[2])*log(STATES[1]/STATES[2]);
  ALGEBRAIC[37] = ( (( CONSTANTS_M[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*STATES[0]*CONSTANTS_M[2]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))*(exp(( (STATES[0] - ALGEBRAIC[35])*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000))/(exp(( STATES[0]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000);
  if(!isfinite(ALGEBRAIC[37]))
    ALGEBRAIC[37] = ( (( CONSTANTS_M[8]*STATES[3]*STATES[3]*STATES[3]*( 0.900000*STATES[4]+ 0.100000*STATES[5])*STATES[1]*CONSTANTS_M[2]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))*(exp((-ALGEBRAIC[35]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1])) - 1.00000))/(CONSTANTS_M[2]/(CONSTANTS_M[0]*CONSTANTS_M[1]) - 1.00000);
  ALGEBRAIC[50] =  CONSTANTS_M[17]*(STATES[0] - ALGEBRAIC[35]);
  ALGEBRAIC[56] = ( CONSTANTS_M[24]*( STATES[2]*STATES[2]*STATES[2]*STATES[18]*exp(( CONSTANTS_M[2]*STATES[0]*CONSTANTS_M[26])/( CONSTANTS_M[0]*CONSTANTS_M[1])) -  STATES[1]*STATES[1]*STATES[1]*STATES[19]*exp(( (CONSTANTS_M[26] - 1.00000)*STATES[0]*CONSTANTS_M[2])/( CONSTANTS_M[0]*CONSTANTS_M[1]))))/(1.00000+ CONSTANTS_M[25]*( STATES[1]*STATES[1]*STATES[1]*STATES[19]+ STATES[2]*STATES[2]*STATES[2]*STATES[18]));
  ALGEBRAIC[39] = STATES[6]/(STATES[6]+CONSTANTS_M[11]);
  ALGEBRAIC[41] =  CONSTANTS_M[9]*STATES[7]*( ALGEBRAIC[39]*STATES[8]+ (1.00000 - ALGEBRAIC[39])*STATES[9])*(STATES[0] - CONSTANTS_M[10]);
  ALGEBRAIC[51] =  (( CONSTANTS_M[0]*CONSTANTS_M[1])/( 2.00000*CONSTANTS_M[2]))*log(STATES[18]/STATES[19]);
  ALGEBRAIC[52] =  CONSTANTS_M[18]*(STATES[0] - ALGEBRAIC[51]);
  ALGEBRAIC[55] = ( CONSTANTS_M[22]*STATES[19])/(STATES[19]+CONSTANTS_M[23]);
  RATES[18] = (CONSTANTS_M[38] - STATES[18])/CONSTANTS_M[36]+((ALGEBRAIC[41]+ALGEBRAIC[52]+ALGEBRAIC[55]) -  2.00000*ALGEBRAIC[56])/( 2.00000*CONSTANTS_M[33]*CONSTANTS_M[2]);
  RATES[1] = (CONSTANTS_M[37] - STATES[1])/CONSTANTS_M[34]+(ALGEBRAIC[37]+ALGEBRAIC[50]+ 3.00000*ALGEBRAIC[56]+ 3.00000*ALGEBRAIC[54]+CONSTANTS_M[28])/( CONSTANTS_M[33]*CONSTANTS_M[2]);
  ALGEBRAIC[58] = ( (STATES[6] - STATES[19])*2.00000*CONSTANTS_M[30]*CONSTANTS_M[2])/CONSTANTS_M[31];
  RATES[6] = - (ALGEBRAIC[41]+ALGEBRAIC[58])/( 2.00000*CONSTANTS_M[30]*CONSTANTS_M[2]);
  ALGEBRAIC[57] =  (10.0000/(1.00000+( 9.13652*pow(1.00000, 0.477811))/pow(CONSTANTS_M[27], 0.477811)))*(0.0517000+0.451600/(1.00000+exp((STATES[0]+59.5300)/17.1800)))*(STATES[0] - ALGEBRAIC[43])*CONSTANTS_M[3];
  ALGEBRAIC[59] = (ALGEBRAIC[37]+ALGEBRAIC[41]+ALGEBRAIC[44]+ALGEBRAIC[45]+ALGEBRAIC[46]+ALGEBRAIC[49]+ALGEBRAIC[47]+ALGEBRAIC[50]+ALGEBRAIC[52]+ALGEBRAIC[54]+ALGEBRAIC[55]+ALGEBRAIC[56]+ALGEBRAIC[57])/CONSTANTS_M[3]+ALGEBRAIC[24];
  RATES[0] =  - ALGEBRAIC[59]*1000.00;
  ALGEBRAIC[60] =  200000.*STATES[19]*(1.00000 - STATES[20]) -  476.000*STATES[20];
  RATES[20] = ALGEBRAIC[60];
  ALGEBRAIC[61] =  78400.0*STATES[19]*(1.00000 - STATES[21]) -  392.000*STATES[21];
  RATES[21] = ALGEBRAIC[61];
  ALGEBRAIC[62] =  200000.*STATES[19]*((1.00000 - STATES[22]) - STATES[23]) -  6.60000*STATES[22];
  RATES[22] = ALGEBRAIC[62];
  ALGEBRAIC[63] =  0.0800000*ALGEBRAIC[61]+ 0.160000*ALGEBRAIC[62]+ 0.0450000*ALGEBRAIC[60];
  RATES[24] = ALGEBRAIC[63];
  ALGEBRAIC[67] = ( CONSTANTS_M[40]*(STATES[19]/CONSTANTS_M[41] - ( CONSTANTS_M[43]*CONSTANTS_M[43]*STATES[26])/CONSTANTS_M[42]))/((STATES[19]+CONSTANTS_M[41])/CONSTANTS_M[41]+( CONSTANTS_M[43]*(STATES[26]+CONSTANTS_M[42]))/CONSTANTS_M[42]);
  ALGEBRAIC[64] = STATES[29]/(STATES[29]+0.250000);
  ALGEBRAIC[65] =  ALGEBRAIC[64]*ALGEBRAIC[64];
  ALGEBRAIC[66] =  CONSTANTS_M[44]*ALGEBRAIC[65]*(STATES[25] - STATES[19]);
  RATES[19] = - ((ALGEBRAIC[52]+ALGEBRAIC[55]+ALGEBRAIC[67]) - (ALGEBRAIC[58]+ALGEBRAIC[66]+ 2.00000*ALGEBRAIC[56]))/( 2.00000*CONSTANTS_M[29]*CONSTANTS_M[2]) -  1.00000*ALGEBRAIC[63];
  ALGEBRAIC[68] = ( (STATES[26] - STATES[25])*2.00000*CONSTANTS_M[46]*CONSTANTS_M[2])/CONSTANTS_M[48];
  RATES[26] = (ALGEBRAIC[67] - ALGEBRAIC[68])/( 2.00000*CONSTANTS_M[45]*CONSTANTS_M[2]);
  ALGEBRAIC[69] =  480.000*STATES[25]*(1.00000 - STATES[27]) -  400.000*STATES[27];
  RATES[27] = ALGEBRAIC[69];
  RATES[25] = (ALGEBRAIC[68] - ALGEBRAIC[66])/( 2.00000*CONSTANTS_M[46]*CONSTANTS_M[2]) -  31.0000*ALGEBRAIC[69];

}



__global__ void ODEstepRL_Paci(PDEFIELD_TYPE dt, PDEFIELD_TYPE ddt, double thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype, PDEFIELD_TYPE* next_stepsize, PDEFIELD_TYPE stepsize_min, PDEFIELD_TYPE eps, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE pacing_duration, PDEFIELD_TYPE pacing_strength){
  /* Ordinary Differential Equation step Runge Kutta Adaptive
  Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and
  adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n]
  at the starting value of the independent variable x. Also input are the stepsize to be attempted
  htry, the required accuracy eps, and the vector yscal[1..n] against which the error is
  scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was
  actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied
  routine that computes the right-hand side derivatives. */
  
  PDEFIELD_TYPE end_time = dt + thetime;
  int nr_iterations = int(dt/ddt);

  PDEFIELD_TYPE y[ARRAY_SIZE];
  PDEFIELD_TYPE y_new[ARRAY_SIZE];
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  bool celltype2 = false;
  int i;
  
  int id = blockDim.x * blockIdx.x + threadIdx.x; 
  if (id < sizex * sizey){
    if (celltype[id] < 1){
      for (i = 0; i < layers; i++) //fill with current PDE values
        alt_PDEvars[i*sizex*sizey + id]= PDEvars[i*sizex*sizey + id];
    }    
    else{
      for (i=0;i<layers;i++) 
        y[i]=PDEvars[i*sizex*sizey + id];
      for (int it = 0; it < nr_iterations; it ++){
        if (it != 0)
          for (i=0;i<layers;i++)
            y[i] = y_new[i];
        celltype2 = false; 
        if (celltype[id] == 2)
          celltype2 = true;    
        derivsPaci_RL(thetime,ddt,y,y_new,dydt,celltype2,pacing_interval,pacing_duration,pacing_strength, id);


        thetime = thetime + ddt;
      }
      for (i=0;i<layers;i++) 
        alt_PDEvars[i*sizex*sizey+id] = y_new[i];
    }
  }
}

__global__ void ODEstepFE(PDEFIELD_TYPE dt, PDEFIELD_TYPE ddt, double thetime, int layers, int sizex, int sizey, PDEFIELD_TYPE* PDEvars, PDEFIELD_TYPE* alt_PDEvars, int* celltype, int* sigmafield, PDEFIELD_TYPE* next_stepsize, PDEFIELD_TYPE stepsize_min, PDEFIELD_TYPE pacing_interval, PDEFIELD_TYPE I_Na_factor, PDEFIELD_TYPE I_f_factor, PDEFIELD_TYPE I_Kr_factor){
  //PDEFIELD_TYPE ddt = 2e-7; //for couplingcoefficient 1e-4 
  //PDEFIELD_TYPE ddt = 1e-6; //for couplingcoefficient 1e-5
  int nr_of_iterations = round(dt/ddt);
  if (fabs(dt/ddt - nr_of_iterations) > 0.001)
    printf("dt and ddt do not divide properly!");
  PDEFIELD_TYPE begin_time,stepsize_next,stepsize_did,stepsize, end_time;
  PDEFIELD_TYPE yscal[ARRAY_SIZE];
  PDEFIELD_TYPE y[ARRAY_SIZE];
  PDEFIELD_TYPE y_new[ARRAY_SIZE];
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE current_time;
  PDEFIELD_TYPE MaxTimeError = 5e-7;
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
      for (int it = 0; it < nr_of_iterations; it++){

        overshot = false;
        if (celltype2)
          derivsFabbriSeveri(current_time,y,dydt,pacing_interval, I_f_factor, I_Kr_factor, id);
        else{
          derivsMaleckar(current_time,y,dydt,pacing_interval, I_Na_factor, id, 0);
        }
        //derivsFitzHughNagumo(current_time,y,dydt,celltype2, sigmafield, pacing_interval,pacing_duration,pacing_strength, id, FHN_interval_beats, FHN_pulse_duration, FHN_pulse_strength,  a, b, tau, FHN_a, FHN_b, FHN_tau);
        current_time += ddt;
        if (it == nr_of_iterations-1) { //Are we done?
          for (i=0;i<layers;i++) {
            alt_PDEvars[i*sizex*sizey + id] = y[i]+ddt*dydt[i];
          }
        }
        else{  
          for (i=0;i<layers;i++) {
            y[i]=y[i]+ddt*dydt[i];  
          }
        }
      }
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

__global__ void CheckActivationTimes(PDEFIELD_TYPE thetime, PDEFIELD_TYPE* activation_times, int* sigmafield,  PDEFIELD_TYPE* PDEvars, int sizex, int sizey){
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  for (int id = index; id < sizex*sizey; id += stride){
    if (PDEvars[id] > 0 && activation_times[id] == 0 && sigmafield[id] > 0)
      activation_times[id] = thetime;
  }
}

void ComputeQthr(PDEFIELD_TYPE* y, PDEFIELD_TYPE t_A, PDEFIELD_TYPE ddt, PDEFIELD_TYPE &Q_thr){
  cout << "t_A = " << t_A << endl;
  PDEFIELD_TYPE activation_strength = 50;
  PDEFIELD_TYPE stepsize_activation = activation_strength/2;
  PDEFIELD_TYPE elapsed_time = 0;
  PDEFIELD_TYPE dydt[ARRAY_SIZE];
  PDEFIELD_TYPE y_copy[ARRAY_SIZE];


  while(stepsize_activation > 1e-7){
    cout << stepsize_activation << endl;
    for (int i = 0; i < ARRAY_SIZE; i ++){
      y_copy[i] = y[i];
    }
    elapsed_time = 0;
    while (y_copy[0] < 0 && elapsed_time < 1){
      //derivsMaleckar_host(elapsed_time,y_copy,dydt,0, 0, -activation_strength, t_A); //linear
      derivsMaleckar_host(elapsed_time,y_copy,dydt,0, 0, -activation_strength, t_A);  //parabolic
      elapsed_time += ddt;
      for (int i = 0; i < ARRAY_SIZE; i++) {
        y_copy[i]=y_copy[i]+ddt*dydt[i];  
      }
      //cout  << "y_copy[0] = " << y_copy[0] << " after " << elapsed_time << " has elapsed." << endl;

    }
    if (y_copy[0] < 0){
      activation_strength = activation_strength + stepsize_activation;
      cout << "For activation_strength = " << activation_strength << ". Failure after " << elapsed_time << " seconds, y_copy[0] = " << y_copy[0] << endl;
    }
    else{
      activation_strength = activation_strength - stepsize_activation;
      cout << "For activation_strength = " << activation_strength << ". Success after " << elapsed_time << " seconds, y_copy[0] = " << y_copy[0] << endl;
    }
    stepsize_activation = stepsize_activation / 2;
    cout << "activation_strength of " << activation_strength << " gives y_copy[0] = " << y_copy[0] << endl;
  }



  Q_thr = activation_strength *1000*50e-9*t_A;
}


void PDE::cuPDEsteps(CellularPotts * cpm, int repeat){
  if (thetime == 0 && par.SF_all)
    InitializeSFComputation(cpm);
  //copy current couplingcoefficient matrix and celltype matrix from host to device
  couplingcoefficient = cpm->getCouplingCoefficient();
  //couplingcoefficient = cpm->getCouplingCoefficient_Gradient();

  //int** cellnumber = cpm -> getSigma(); 
  cudaError_t errSync;
  cudaError_t errAsync;
  celltype = cpm->getTau();
  sigmafield = cpm->getSigma(); 
  cudaMemcpy(d_couplingcoefficient, couplingcoefficient[0], sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice); 
  cudaMemcpy(d_celltype, celltype[0], sizex*sizey*sizeof(int), cudaMemcpyHostToDevice); 
  cudaMemcpy(d_sigmafield, sigmafield[0], sizex*sizey*sizeof(int), cudaMemcpyHostToDevice); 
  cudaMemcpy(d_PDEvars, PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyHostToDevice); 

  
  int nr_blocks = sizex*sizey/par.threads_per_core + 1;
  PDEFIELD_TYPE Cm_maleckar = 50; //in nF
  PDEFIELD_TYPE I_m;
  bool afterdiffusion;

  for (int iteration = 0; iteration < repeat; iteration++){
    if (par.SF_all){
      cuSFChecker();
    }
      //cout << "Iteration = " << iteration << endl;

      //setup matrices for upperdiagonal, diagonal and lower diagonal for both the horizontal and vertical direction, since these remain the same during once MCS
    InitializeDiagonals<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, lowerH, upperH, diagH, lowerV, upperV, diagV, d_couplingcoefficient);
    cudaDeviceSynchronize();
    errSync  = cudaGetLastError();
    errAsync = cudaDeviceSynchronize();
    if (errSync != cudaSuccess) 
      printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
    if (errAsync != cudaSuccess)
      printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));


    cuODEstep();
    afterdiffusion = false;
    if (par.SF_all){
      cuCopyVoltageForSF(afterdiffusion);
    }
  
    cuHorizontalADIstep();
    afterdiffusion = true;
    if (par.SF_all){
      cuCopyVoltageForSF(afterdiffusion);
    }

    //increase time by dt/2
    thetime = thetime + dt/2;  
    cuODEstep();
    afterdiffusion = false;
    if (par.SF_all){
      cuCopyVoltageForSF(afterdiffusion);
    }


    cuVerticalADIstep();
    afterdiffusion = true;
    if (par.SF_all){
      cuCopyVoltageForSF(afterdiffusion);
    }
      
    //cudaMemcpy(alt_PDEvars, d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    //cout << "After second FE step, alt_PDEvars[23885] = " << alt_PDEvars[23885] << endl;
    if (par.SF_all){
      if (SF_start_one && !SF_end_one && par.SF_one_pixel)
        cuComputeSFOne();
      
      if (SF_in_progress && !SF_all_done && par.SF_all)
        cuWriteSFData();
    }
    

    
    //increase time by dt/2
    thetime = thetime + dt/2; 

    if (par.activation_times){    
      CheckActivationTimes<<<par.number_of_cores, par.threads_per_core>>>(thetime, d_Activation_times_array, d_sigmafield,  d_PDEvars, sizex, sizey);
    }   
  }
  cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  cuPDEVarsToFiles();
  if (par.activation_times){ 
    cudaMemcpy(Activation_times_array, d_Activation_times_array, sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);   
    cuWriteActivationTimes();
  }
}

void PDE::cuODEstep(){
  //Do an ODE step of size dt/2
  cudaError_t errSync;
  cudaError_t errAsync;
  //ODEstepRL_Paci<<<nr_blocks, par.threads_per_core>>>(dt/2, ddt, thetime, layers, sizex, sizey, d_PDEvars, d_alt_PDEvars, d_celltype, next_stepsize, min_stepsize, par.eps, pacing_interval, par.pacing_duration, par.pacing_strength);
  //ODEstepRKA<<<par.number_of_cores, par.threads_per_core>>>(dt/2, thetime, layers, sizex, sizey, d_PDEvars, d_alt_PDEvars, d_celltype, next_stepsize, min_stepsize, par.eps, pacing_interval, par.pacing_duration, par.pacing_strength);
  ODEstepFE<<<par.number_of_cores, par.threads_per_core>>>(dt/2, ddt, thetime, layers, sizex, sizey, d_PDEvars, d_alt_PDEvars, d_celltype, d_sigmafield, next_stepsize, min_stepsize, pacing_interval, par.I_f_factor, par.I_Kr_factor, par.I_Na_factor);
  //CopyOriginalToAltPDEvars<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, layers, d_PDEvars, d_alt_PDEvars);

  //cudaMemcpy(alt_PDEvars, d_alt_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
  //cout << "After second FE step, alt_PDEvars[4305] = " << alt_PDEvars[4305] << endl;
  cuErrorChecker(errSync, errAsync);
}

void PDE::cuHorizontalADIstep(){
  //Do a horizontal ADI sweep of size dt/2
  cudaError_t errSync;
  cudaError_t errAsync;
  InitializeHorizontalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BH, d_couplingcoefficient, d_alt_PDEvars);
  cuErrorChecker(errSync, errAsync);
  #ifdef PDEFIELD_DOUBLE
    statusH = cusparseDgtsvInterleavedBatch(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, pbufferH);
  #else
    statusH = cusparseSgtsvInterleavedBatch(handleH, 0, sizex, lowerH, diagH, upperH, BH, sizey, pbufferH);
  #endif
  if (statusH != CUSPARSE_STATUS_SUCCESS)
  {
    cout << statusH << endl;
  }
  NewPDEfieldH0<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, BH, d_PDEvars);    
  cuErrorChecker(errSync, errAsync);
  NewPDEfieldOthers<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, layers, BV, d_PDEvars, d_alt_PDEvars); //////
  cuErrorChecker(errSync, errAsync);

  //cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
  //cout << "After second FE step, PDEvars[4305] = " << PDEvars[4305] << endl;

}

void PDE::cuVerticalADIstep(){
  //Do a vertical ADI sweep of size dt/2
  cudaError_t errSync;
  cudaError_t errAsync;
  InitializeVerticalVectors<<<par.number_of_cores, par.threads_per_core>>>(sizex, sizey, 2/dt, dx2, BV, d_couplingcoefficient, d_alt_PDEvars);
  cuErrorChecker(errSync, errAsync);
  #ifdef PDEFIELD_DOUBLE
    statusV = cusparseDgtsvInterleavedBatch(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, pbufferV);
  #else
    statusV = cusparseSgtsvInterleavedBatch(handleV, 0, sizey, lowerV, diagV, upperV, BV, sizex, pbufferV);
  #endif
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
  cuErrorChecker(errSync, errAsync);

  //cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
  //cout << "After second FE step, PDEvars[4305] = " << PDEvars[4305] << endl;

}

void PDE::cuErrorChecker(cudaError_t errSync, cudaError_t errAsync){
  errSync  = cudaGetLastError();
  errAsync = cudaDeviceSynchronize();
  if (errSync != cudaSuccess) 
    printf("Sync kernel error: %s\n", cudaGetErrorString(errSync));
  if (errAsync != cudaSuccess)
    printf("Async kernel error: %s\n", cudaGetErrorString(errAsync));
}

void PDE::cuPDEVarsToFiles(){

  
  int number_of_measurements = 10;
  int measure_loc;
  ofstream myfile;
  char fname[200];
  for (int measurement = 1; measurement <= number_of_measurements; measurement++){ 
    //evenly spread out measurement locations between (0.5 + 10)*sizey and (sizex-11 + 0.5)*sizey
    measure_loc = int ((0.5 + int(10 * (number_of_measurements-measurement)/(number_of_measurements-1) + (sizex-11)*(measurement-1)/(number_of_measurements-1)))*sizey);  
    sprintf(fname,"location_%03d.txt",measurement);
    myfile.open(fname, std::ios_base::app);
    myfile << thetime << ",";
      for (int i = 0; i < layers; i++)
        myfile << PDEvars[measure_loc+sizex*sizey*i] << ",";
    myfile << endl;
    myfile.close();  
  }
  /*
  int measure_loc;
  ofstream myfile;
  char fname[200];
  int measurement = 1;

  measure_loc = int ((0.5 + 343)*sizey);  
  sprintf(fname,"location_%03d.txt",measurement);
  myfile.open(fname, std::ios_base::app);
  myfile << thetime << ",";
    for (int i = 0; i < layers; i++)
      myfile << PDEvars[measure_loc+sizex*sizey*i] << ",";
  myfile << endl;
  myfile.close();

  measurement++;

  measure_loc = int ((0.5 + 418)*sizey);
  sprintf(fname,"location_%03d.txt",measurement);
  myfile.open(fname, std::ios_base::app);
  myfile << thetime << ",";
    for (int i = 0; i < layers; i++)
      myfile << PDEvars[measure_loc+sizex*sizey*i] << ",";
  myfile << endl;
  myfile.close();
  if (PDEvars[measure_loc] > 0){
    //cout << "We have all data and can safely stop the simulation now" << endl;
    //exit(0);
  }*/

  cout << "PDEvars["<< int(sizey*10.5)<< "] = " << PDEvars[int(sizey*10.5)] << 
  ", PDEvars["<< sizex*sizey-int(sizey*10.5)<< "] = " << PDEvars[sizex*sizey-int(sizey*10.5)] << " and time = " << thetime << endl;
  if (!(PDEvars[int(sizex/2*sizey + 0.5 * sizey) ]>-1000000000 && PDEvars[int(sizex/2*sizey + 0.5 * sizey)] < 1000000000)){
    cout << "We encountered a NaN error. Abort the program. \n";
    exit(1);
  }
}

void PDE::cuSFChecker(){
  //Check whether or not we should be doing SF computations, for one pixel or all pixels.

  //int SF_locator_one = -1;
  string file_loc_base = string(par.datadir) + "/Q_tot_data_";
  string file_loc;
  ofstream Q_tot_file;
  
  
  //Keeps track of which pixels should START exporting data for SF computations
  if (par.SF_all && !SF_all_done){
    for (int i = 0; i < sizex*sizey; i++)
      if (PDEvars[i] > -70 && thetime > 0.2 &&  !SF_start_array[i]){
        SF_in_progress = true;
        SF_start_array[i] = true;
        file_loc = file_loc_base + to_string(i/par.sizey) + "_" + to_string(i%par.sizey) + ".txt";
        Q_tot_file.open(file_loc, std::ios_base::app);
        Q_tot_file << ddt << endl;
        for (int j = 0; j < layers; j++)
        Q_tot_file << PDEvars[i + j*sizex*sizey] << endl;
        Q_tot_file.close();
      }
  }

  //Keeps track of which pixels should STOP exporting data for SF computations
  if (par.SF_all && SF_in_progress && !SF_all_done){
    SF_all_done = true;
    for (int i = 0; i < sizex*sizey; i++){
      if ((PDEvars[i] < -70) &&  SF_start_array[i] && celltype[0][i] == 1 && !SF_end_array[i]){
        SF_in_progress = true;
        SF_end_array[i] = true;
      }
      if (!SF_end_array[i])
        SF_all_done = false;
    }
  }
  
  //Keeps track of the pixel for exporting data for SF for a single pixel
  if (PDEvars[SF_locator_one] > -70 && thetime > 0.05 &&  !SF_start_one && par.SF_one_pixel){
    SF_start_one = true;
    for (int k = 0; k < ARRAY_SIZE; k++){
      copy_for_SF[k] = PDEvars[SF_locator_one + sizex*sizey*k];
      cout << "PDEvars[" << k << "] = " << copy_for_SF[k] << endl;
    }
  }

}

void PDE::cuCopyVoltageForSF(bool afterdiffusion){
  PDEFIELD_TYPE Q_factor = 50e-9;
  if (afterdiffusion){
    //If SF computations need to happen, copy voltage from d_PDEvars
    if (SF_in_progress && !SF_all_done && par.SF_all) {
      cudaMemcpy(PDEvars, d_PDEvars, sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
      for (int i = 0; i < sizex*sizey; i++)
        if (SF_start_array[i] && !SF_end_array[i])
          SF_Q_tot_array[i] += (PDEvars[i] - alt_PDEvars[i])*Q_factor;

    }
    if (SF_start_one && !SF_end_one && par.SF_one_pixel){
      cudaMemcpy(PDEvars, d_PDEvars, layers*sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
      Q_tot += (PDEvars[SF_locator_one] - alt_PDEvars[SF_locator_one])*Q_factor;
    }
  }
  else{
    //If SF computations need to happen, copy voltage from d_alt_PDEvars
    if (SF_in_progress && !SF_all_done && par.SF_all) 
      cudaMemcpy(alt_PDEvars, d_alt_PDEvars, sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    if (SF_start_one && !SF_end_one && par.SF_one_pixel){
      cudaMemcpy(alt_PDEvars, d_alt_PDEvars, sizex*sizey*sizeof(PDEFIELD_TYPE), cudaMemcpyDeviceToHost);
    }
  }
}

void PDE::cuComputeSFOne(){
  Q_tot_store_one[t_A_counter] = Q_tot;
  t_A_counter++;
  cout << "counter = " << t_A_counter << endl;
  if (PDEvars[SF_locator_one] < -70){
    cout << "final counter = " << t_A_counter << endl;
    int stepper = t_A_counter/4;
    int current_loc = t_A_counter/2;
    double SF_lower;
    double SF_upper;
    while (stepper != 0){
      ComputeQthr(copy_for_SF, current_loc*ddt, ddt, Q_thr);
      SF_lower = Q_tot_store_one[current_loc]/Q_thr;
      ComputeQthr(copy_for_SF, (current_loc+1)*ddt, ddt, Q_thr);
      SF_upper = Q_tot_store_one[current_loc+1]/Q_thr;
      if(SF_lower < SF_upper)
        current_loc += stepper;
      else 
        current_loc -= stepper;
      stepper /= 2;
      if (stepper == 0){
        if (SF_lower < SF_upper)
          cout << "SF = " << SF_upper << endl;
        else
          cout << "SF = " << SF_lower << endl;
      }
    }
  }
}

void PDE::cuWriteActivationTimes(){
  ofstream Activation_times_file;
  string file_loc = string(par.datadir) + "/Activation_times_file.txt";
  Activation_times_file.open(file_loc, std::ios_base::app);
  for (int i = 0; i < par.sizex*par.sizey; i++){
    if (Activation_times_array[i] > 0 && !Activation_times_array_written[i]){ 
      Activation_times_array_written[i] = true;  
      Activation_times_file << i/par.sizey << ", " << i%par.sizey << ", " << Activation_times_array[i] << endl;
    }
  }
  Activation_times_file.close();
}

void PDE::cuWriteSFData(){
  ofstream Q_tot_file;
  string file_loc_base = string(par.datadir) + "/Q_tot_data_";
  string file_loc;
  for (int i = 0; i < par.sizex*par.sizey; i++){
    if (SF_start_array[i] && !SF_end_array[i]){
      file_loc = file_loc_base + to_string(i/par.sizey) + "_" + to_string(i%par.sizey) + ".txt";
      Q_tot_file.open(file_loc, std::ios_base::app);
      Q_tot_file << SF_Q_tot_array[i] << endl;
      Q_tot_file.close();
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
  if (val > -100 && val < 100){
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
