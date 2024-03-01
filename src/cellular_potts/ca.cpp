/*

Copyright 1995-2006 Roeland Merks, Nick Savill

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

/* CA.cpp: implementation of Glazier & Graner's Cellular Potts Model */

// This code derives from a Cellular Potts implementation written around 1995
// by Nick Savill


#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

#include "sticky.hpp"
#include "random.hpp"
#include "ca.hpp"
#include "parameter.hpp"
#include "dish.hpp"
#include "sqr.hpp"
#include "crash.hpp"
#include "hull.hpp"
#include "graph.hpp"

#define ZYGFILE(Z) <Z.xpm>
#define XPM(Z) Z##_xpm
#define ZYGXPM(Z) XPM(Z)

#define PI 3.14159265


/* define default zygote */
/* NOTE: ZYGOTE is normally defined in Makefile!!!!!! */
#ifndef ZYGOTE
#define ZYGOTE init
#include "xpm/1.xpm"
#else
#include ZYGFILE(ZYGOTE)
#endif

/* STATIC DATA MEMBER INITIALISATION */
double copyprob[BOLTZMANN];

const int CellularPotts::nx[21] = {0, 0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 1, 2, 2, 1, -1, -2, -2, -1};
const int CellularPotts::ny[21] = {0, -1, 0, 1, 0, -1, 1, 1, -1, -2, 0, 2, 0, -2, -1, 1, 2, 2, 1, -1, -2};

const int CellularPotts::nbh_level[4] = {0, 4, 8, 20};
int CellularPotts::shuffleindex[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

extern Parameter par;

/** PRIVATE **/

using namespace std;
void CellularPotts::BaseInitialisation(vector<Cell> *cells)
{
  CopyProb(par.T);
  cell = cells;
  if (par.neighbours >= 1 && par.neighbours <= 4)
    n_nb = nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4]).";
}

CellularPotts::CellularPotts(vector<Cell> *cells,
                             const int sx, const int sy)
{
  sigma = 0;
  frozen = false;
  thetime = 0;
  zygote_area = 0;

  edgelist = nullptr;
  orderedgelist = nullptr;

  BaseInitialisation(cells);
  sizex = sx;
  sizey = sy;

  AllocateSigma(sx, sy);
  AllocateTau(sx, sy);



  if (par.micropatternmask != "None")
  {
    StoreMask();
    for (int x = 0; x < sizex; x++)
    {
      for (int y = 0; y < sizey; y++)
      {
        if (!mask[x][y])
        {
          sigma[x][y] = -1;
          tau[x][y] = -1;
        }
      }
    }
  }
  else
  {
    // fill borders with special border state
    for (int x = 0; x < sizex; x++)
    {
      sigma[x][0] = -1;
      sigma[x][sizey - 1] = -1;
      tau[x][0] = -1;
      tau[x][sizey - 1] = -1;
    }
    for (int y = 0; y < sizey; y++)
    {
      sigma[0][y] = -1;
      sigma[sizex - 1][y] = -1;
      tau[0][y] = -1;
      tau[sizex - 1][y] = -1;
    }
  }

  if (par.neighbours >= 1 && par.neighbours <= 4)
    n_nb = nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

CellularPotts::CellularPotts(void)
{

  sigma = 0;
  sizex = 0;
  sizey = 0;
  frozen = false;
  thetime = 0;
  zygote_area = 0;

  edgelist = nullptr;
  orderedgelist = nullptr;

  CopyProb(par.T);

  if (par.micropatternmask != string("None"))
  {
    for (int x = 0; x < sizex; x++)
    {
      for (int y = 0; y < sizey; y++)
      {
        if (!mask[x][y])
        {
          sigma[x][y] = -1;
          tau[x][y] = -1;
        }
      }
    }
  }

  else
  {
    // fill borders with special border state
    for (int x = 0; x < sizex; x++)
    {
      sigma[x][0] = -1;
      sigma[x][sizey - 1] = -1;
      tau[x][0] = -1;
      tau[x][sizey - 1] = -1;
    }
    for (int y = 0; y < sizey; y++)
    {
      sigma[0][y] = -1;
      sigma[sizex - 1][y] = -1;
      tau[0][y] = -1;
      tau[sizex - 1][y] = -1;
    }
  }
  if (par.neighbours >= 1 && par.neighbours <= 4)
    n_nb = nbh_level[par.neighbours];
  else
    throw "Panic in CellularPotts: parameter neighbours invalid (choose [1-4])";
}

// destructor (virtual)
CellularPotts::~CellularPotts(void)
{
  if (sigma)
  {
    free(sigma[0]);
    free(sigma);
    sigma = 0;
  }

  if (tau)
  {
    free(tau[0]);
    free(tau);
    tau = 0;
  }

  if (edgelist)
  {
    delete edgelist;
  }

  if (orderedgelist)
  {
    delete orderedgelist;
  }
}

void CellularPotts::AllocateSigma(int sx, int sy)
{

  sizex = sx;
  sizey = sy;

  sigma = (int **)malloc(sizex * sizeof(int *));
  if (sigma == NULL)
    MemoryWarning();

  sigma[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (sigma[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      sigma[i] = sigma[i - 1] + sizey;
  }

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      sigma[0][i] = 0;
  }
}

void CellularPotts::AllocateTau(int sx, int sy)
{

  sizex = sx;
  sizey = sy;

  tau = (int **)malloc(sizex * sizeof(int *));
  if (tau == NULL)
    MemoryWarning();

  tau[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (tau[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      tau[i] = tau[i - 1] + sizey;
  }

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      tau[0][i] = 0;
  }
}

void CellularPotts::AllocateNumberOfEdges(int sx, int sy)
{

  sizex = sx;
  sizey = sy;

  numberofedges = (int **)malloc(sizex * sizeof(int *));
  if (numberofedges == NULL)
    MemoryWarning();

  numberofedges[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (numberofedges[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      numberofedges[i] = numberofedges[i - 1] + sizey;
  }

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      numberofedges[0][i] = 0;
  }
}

void CellularPotts::AllocateCouplingCoefficient(int sx, int sy)
{
  couplingcoeffcient_allocated = true;
  sizex = sx;
  sizey = sy;

  couplingcoefficient = (PDEFIELD_TYPE **)malloc(sizex * sizeof(PDEFIELD_TYPE *));
  if (couplingcoefficient == NULL)
    MemoryWarning();

  couplingcoefficient[0] = (PDEFIELD_TYPE *)malloc(sizex * sizey * sizeof(PDEFIELD_TYPE));
  if (couplingcoefficient[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      couplingcoefficient[i] = couplingcoefficient[i - 1] + sizey;
  }

  /* Clear coupling coeff plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      couplingcoefficient[0][i] = 0;
  }
}

void CellularPotts::AllocateCouplingCoefficient_Gradient(int sx, int sy)
{

  sizex = sx;
  sizey = sy;

  couplingcoefficient_gradient = (PDEFIELD_TYPE **)malloc(sizex * sizeof(PDEFIELD_TYPE *));
  if (couplingcoefficient_gradient == NULL)
    MemoryWarning();

  couplingcoefficient_gradient[0] = (PDEFIELD_TYPE *)malloc(sizex * sizey * sizeof(PDEFIELD_TYPE));
  if (couplingcoefficient_gradient[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      couplingcoefficient_gradient[i] = couplingcoefficient_gradient[i - 1] + sizey;
  }

  /* Clear coupling coeff plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      couplingcoefficient_gradient[0][i] = 0;
  }
}

void CellularPotts::AllocateMask(int sx, int sy)
{
  sizex = sx;
  sizey = sy;

  mask = (int **)malloc(sizex * sizeof(int *));
  if (mask == NULL)
    MemoryWarning();

  mask[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (mask[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      mask[i] = mask[i - 1] + sizey;
  }

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      mask[0][i] = 0;
  }
}

void CellularPotts::InitializeMatrix(Dish &beast)
{
  // sizex; sizey=sy;

  matrix = (int **)malloc(sizex * sizeof(int *));
  if (matrix == NULL)
    MemoryWarning();

  matrix[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (matrix[0] == NULL)
    MemoryWarning();

  {
    for (int i = 1; i < sizex; i++)
      matrix[i] = matrix[i - 1] + sizey;
  }

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      matrix[0][i] = 0;
  }
}

void CellularPotts::IndexShuffle()
{
  int i;
  int temp;
  int index1, index2;

  for (i = 0; i < 9; i++)
  {
    index1 = RandomNumber(8);
    index2 = RandomNumber(8);

    temp = shuffleindex[index1];
    shuffleindex[index1] = shuffleindex[index2];
    shuffleindex[index2] = temp;
  }
}

void CellularPotts::InitializeEdgeList(void)
{
  AllocateNumberOfEdges(par.sizex, par.sizey);
  edgelist = new int[(par.sizex - 2) * (par.sizey - 2) * nbh_level[par.neighbours]];
  orderedgelist = new int[(par.sizex - 2) * (par.sizey - 2) * nbh_level[par.neighbours]];
  sizeedgelist = 0;
  int pixel;
  int neighbour;
  int x, y;
  int xp, yp;
  int c, cp;
  //Initialize both edgelist and orderedgelist to have -1 everywhere
  for (int k = 0; k < (par.sizex - 2) * (par.sizey - 2) * nbh_level[par.neighbours]; k++)
  {
    edgelist[k] = -1;
    orderedgelist[k] = -1;
  }
  for (int k = 0; k < (par.sizex - 2) * (par.sizey - 2) * nbh_level[par.neighbours]; k++)
  {
    //Loop over all edges
    //Outermost loop is over the y-coordinate
    //Middle loop is over the x-coordinate
    //Innermost loop is over the neighbours.
    pixel = k / nbh_level[par.neighbours];
    neighbour = k % nbh_level[par.neighbours] + 1;
    x = pixel % (sizex - 2) + 1;
    y = pixel / (sizex - 2) + 1;
    c = sigma[x][y];
    xp = nx[neighbour] + x; //which cell is the neighbour?
    yp = ny[neighbour] + y; //which cell is the neighbour?

    if (par.micropatternmask != string("None"))
    {
      if (!mask[x][y])
        c = -1;
      if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
        cp = -1;
      else if (!mask[xp][yp])
        cp = -1;
      else
        cp = sigma[xp][yp];
    }



    else
    {
      if (par.periodic_boundaries)
      {
        // since we are asynchronic, we cannot just copy the borders once
        // every MCS
        if (xp <= 0)
          xp = sizex - 2 + xp;
        if (yp <= 0)
          yp = sizey - 2 + yp;
        if (xp >= sizex - 1)
          xp = xp - sizex + 2;
        if (yp >= sizey - 1)
          yp = yp - sizey + 2;
        cp = sigma[xp][yp];
      }
      else if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
        cp = -1;
      else
        cp = sigma[xp][yp];
    }
    if (cp != c && cp != -1 && c != -1)
    { 
      if(!par.second_layer || tau[x][y] == tau[xp][yp]){ //If we do not use a second layer, or if both cell types are equal)
        //if a pixel and its neighbour have a different sigma, add a unique 
        //interger to edgelist       
        edgelist[k] = sizeedgelist;
        //also add a unique integer to the end of orderedgelist, making a bijection between the lists
        orderedgelist[sizeedgelist] = k;
        sizeedgelist ++;
        numberofedges[x][y]++;
      }
    }
  }
}


void CellularPotts::InitializeCouplingCoefficient(void)
{
  cout << "Initialize Coupling coeff\n";
  AllocateCouplingCoefficient(par.sizex, par.sizey);
  for (int x = 0; x < par.sizex; x++)
  {
    for (int y = 0; y < par.sizey; y++)
    {
      if (sigma[x][y] == -1)
        couplingcoefficient[x][y] = par.couplingoffmask;
      else if (sigma[x][y] == 0)
        couplingcoefficient[x][y] = par.couplingmedium;
      else if (sigma[x][y] > 0)
        couplingcoefficient[x][y] = par.couplingcell;
      if (numberofedges[x][y] > 0){
        if (tau[x][y] == 1){
          couplingcoefficient[x][y] = par.couplingAtrialAtrial;
          for (int i = 1; i < n_nb; i++){
            if (tau[x+nx[i]][y+ny[i]] == 2){
              couplingcoefficient[x][y] = par.couplingAtrialPM;
              couplingcoefficient[x+nx[i]][y+ny[i]]= par.couplingAtrialPM;
            }
          }
        }
        else if (tau[x][y] == 2){  
          couplingcoefficient[x][y] = par.couplingPMPM;
          for (int i = 1; i < n_nb; i++){
            if (tau[x+nx[i]][y+ny[i]] == 1){
              couplingcoefficient[x][y] = par.couplingAtrialPM;
              couplingcoefficient[x+nx[i]][y+ny[i]]= par.couplingAtrialPM;
            }
          }
        }
        if(par.second_layer){
          for (int i = 1; i < n_nb; i++){
            if ((tau[x][y] == 1 && tau[x+nx[i]][y+ny[i]] == 2) || (tau[x][y] == 2 && tau[x+nx[i]][y+ny[i]] == 1)){
              couplingcoefficient[x][y] = par.couplingAtrialPM;
              couplingcoefficient[x+nx[i]][y+ny[i]]= par.couplingAtrialPM;
            } 
          }
        }
      }

    }
  }
}

void CellularPotts::InitializeCouplingCoefficientNoCellularDetail(void)
{
  AllocateCouplingCoefficient(par.sizex, par.sizey);
  for (int x = 0; x < par.sizex; x++)
  {
    for (int y = 0; y < par.sizey; y++)
    {
      if (sigma[x][y] == -1)
        couplingcoefficient[x][y] = par.couplingoffmask;
      else if (sigma[x][y] == 0)
        couplingcoefficient[x][y] = par.couplingmedium;
      else if (sigma[x][y] > 0){
        if (tau[x][y] == 1)
          couplingcoefficient[x][y] = par.couplingAtrialAtrial;
        else if (tau[x][y] == 2)
          couplingcoefficient[x][y] = par.couplingPMPM;
        if (numberofedges[x][y] > 0)
          couplingcoefficient[x][y] = par.couplingAtrialPM;
      }
    }
  }
}

void CellularPotts::InitializeCouplingCoefficient_Gradient(void)
{
  DetectSidesIsthmus();
  bool written = false;
  AllocateCouplingCoefficient_Gradient(par.sizex, par.sizey);
  for (int x = 0; x < par.sizex; x++)
  {
    written = false;
    for (int y = 0; y < par.sizey; y++)
    {
      if (mask[x][y]){
        if (x < left_side_isthmus)
          couplingcoefficient_gradient[x][y] = par.couplingPMPM;
        else if (x > right_side_isthmus) 
          couplingcoefficient_gradient[x][y] = par.couplingAtrialAtrial;
        else{
          int isthmus_length = right_side_isthmus - left_side_isthmus;
          //couplingcoefficient_gradient[x][y] = (1.0-float(x - left_side_isthmus)/isthmus_length) * par.couplingPMPM + float(x - left_side_isthmus)/isthmus_length * par.couplingAtrialAtrial;
          couplingcoefficient_gradient[x][y] = exp((1.0-float(x - left_side_isthmus)/isthmus_length) * log(par.couplingPMPM) + float(x - left_side_isthmus)/isthmus_length * log(par.couplingAtrialAtrial));
          if (!written){
            written = true;
            cout << "couplingcoefficient_gradient[" << x << "][" << y << "] = "<< couplingcoefficient_gradient[x][y] << endl;
          }
        }
      }
    }  
  }
}
double sat(double x)
{
  return x / (par.saturation * x + 1.);
  // return x;
}

int CellularPotts::IsingDeltaH(int x, int y, PDE *PDEfield)
{
  int DH = 0, H_before = 0, H_after = 0;
  int i, sxy;
  int neigh_sxy;
  int J = par.lambda;

  /* Compute energydifference *IF* the flip were to occur */
  sxy = sigma[x][y];

  /* DH due to spin alignment */
#ifdef DBG_KAWASAKI
  std::cerr << "[ x = {" << x << ", " << y << "}, xp = {" << xp << ", " << yp << "}, ";
#endif
  for (i = 1; i <= n_nb; i++)
  {
    int xn, yn;
    xn = x + nx[i];
    yn = y + ny[i];

    if (par.periodic_boundaries)
    {

      // since we are asynchronic, we cannot just copy the borders once
      // every MCS

      if (xn <= 0)
        xn = sizex - 2 + xn;
      if (yn <= 0)
        yn = sizey - 2 + yn;
      if (xn >= sizex - 1)
        xn = xn - sizex + 2;
      if (yn >= sizey - 1)
        yn = yn - sizey + 2;

      neigh_sxy = sigma[xn][yn];

    } // periodic boundaries
    else
    { // closed boundaries

      if (xn <= 0 || yn <= 0 || xn >= sizex - 1 || yn >= sizey - 1)
        neigh_sxy = -1;
      else
        neigh_sxy = sigma[xn][yn];
    }

    if (neigh_sxy == -1)
    { // border
      cerr << "Only periodic boundaries implemented for Kawasaki dynamics sofar.\n";
      exit(1);
      //  DH += (sxyp==0?0:par.border_energy)-
      //  (sxy==0?0:par.border_energy);
    }
    else
    {
      H_before += -J * (sxy == 0 ? -1 : 1) * (neigh_sxy == 0 ? -1 : 1);
      H_after += -J * (sxy == 0 ? 1 : -1) * (neigh_sxy == 0 ? -1 : 1);
    }
  }

  DH = H_after - H_before;

  return DH;
}

int CellularPotts::PottsDeltaH(int x, int y, int new_state)
{
  int DH = 0, H_before = 0, H_after = 0;
  int i, sxy;
  int neigh_sxy;
  int J = par.lambda;

  /* Compute energydifference *IF* the flip were to occur */
  sxy = sigma[x][y];

  /* DH due to spin alignment */

  for (i = 1; i <= n_nb; i++)
  {
    int xn, yn;
    xn = x + nx[i];
    yn = y + ny[i];

    if (par.periodic_boundaries)
    {

      // since we are asynchronic, we cannot just copy the borders once
      // every MCS

      if (xn <= 0)
        xn = sizex - 2 + xn;
      if (yn <= 0)
        yn = sizey - 2 + yn;
      if (xn >= sizex - 1)
        xn = xn - sizex + 2;
      if (yn >= sizey - 1)
        yn = yn - sizey + 2;

      neigh_sxy = sigma[xn][yn];

    } // periodic boundaries
    else
    { // closed boundaries

      if (xn <= 0 || yn <= 0 || xn >= sizex - 1 || yn >= sizey - 1)
        neigh_sxy = -1;
      else
        neigh_sxy = sigma[xn][yn];
    }

    if (neigh_sxy == -1)
    { // border
      cerr << "Only periodic boundaries implemented for Potts dynamics sofar.\n";
      exit(1);
      //  DH += (sxyp==0?0:par.border_energy)-
      //  (sxy==0?0:par.border_energy);
    }
    else
    {
      /*
      H_before += -J*(sxy==0?-1:1)*(neigh_sxy==0?-1:1);
      H_after += -J*(sxy==0?1:-1)*(neigh_sxy==0?-1:1);*/
      H_before += J * ((sxy != neigh_sxy) ? 1 : 0);
      H_after += J * ((new_state != neigh_sxy) ? 1 : 0);
    }
  }
  DH = H_after - H_before;
  return DH;
}

int CellularPotts::KawasakiDeltaH(int x, int y, int xp, int yp, PDE *PDEfield)
{
  int DH = 0, H_before = 0, H_after = 0;
  int i, sxy, sxyp;
  int neigh_sxy, neigh_sxyp;

  /* Compute energydifference *IF* the copying were to occur */
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];

  /* DH due to cell adhesion */
#ifdef DBG_KAWASAKI
  std::cerr << "[ x = {" << x << ", " << y << "}, xp = {" << xp << ", " << yp << "}, ";
#endif
  for (i = 1; i <= n_nb; i++)
  {
    int xn, yn;
    xn = x + nx[i];
    yn = y + ny[i];

    int xpn, ypn;
    xpn = xp + nx[i];
    ypn = yp + ny[i];

    if (par.periodic_boundaries)
    {

      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (xn <= 0)
        xn = sizex - 2 + xn;
      if (yn <= 0)
        yn = sizey - 2 + yn;
      if (xn >= sizex - 1)
        xn = xn - sizex + 2;
      if (yn >= sizey - 1)
        yn = yn - sizey + 2;

      neigh_sxy = sigma[xn][yn];

      if (xpn <= 0)
        xpn = sizex - 2 + xpn;
      if (ypn <= 0)
        ypn = sizey - 2 + ypn;
      if (xpn >= sizex - 1)
        xpn = xpn - sizex + 2;
      if (ypn >= sizey - 1)
        ypn = ypn - sizey + 2;

      neigh_sxyp = sigma[xpn][ypn];

    } // periodic boundaries
    else
    { // closed boundaries

      if (xn <= 0 || yn <= 0 || xn >= sizex - 1 || yn >= sizey - 1)
        neigh_sxy = -1;
      else
        neigh_sxy = sigma[xn][yn];

      if (xpn <= 0 || ypn <= 0 || xpn >= sizex - 1 || ypn >= sizey - 1)
        neigh_sxyp = -1;
      else
        neigh_sxyp = sigma[xpn][ypn];
    }

    if (neigh_sxy == -1)
    { // border
      cerr << "Only periodic boundaries implemented for Kawasaki dynamics sofar.\n";
      exit(1);
      //  DH += (sxyp==0?0:par.border_energy)-
      //  (sxy==0?0:par.border_energy);
    }
    else
    {
      H_before += (*cell)[sxy].EnergyDifference((*cell)[neigh_sxy]) +
                  (*cell)[sxyp].EnergyDifference((*cell)[neigh_sxyp]);
      int aft = (*cell)[sxyp].EnergyDifference((*cell)[neigh_sxy]) +
                (*cell)[sxy].EnergyDifference((*cell)[neigh_sxyp]);
#ifdef DBG_KAWASAKI
      cerr << aft << ", ";
#endif
      H_after += aft;
    }
  }
  H_after += 2 * (*cell)[sxy].EnergyDifference((*cell)[sxyp]);
#ifdef DBG_KAWASAKI
  cerr << H_after << ", " << H_before << " ]";
#endif
  DH = H_after - H_before;
  // the rest we will do later - in any case no volume constraint :-)
  return DH;
}

int CellularPotts::DeltaH(int x, int y, int xp, int yp, PDE *PDEfield)
{
  int DH = 0;
  int i, sxy, sxyp;
  int neighsite;

  /* Compute energydifference *IF* the copying were to occur */
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];

  /* DH due to cell adhesion */
  for (i = 1; i <= n_nb; i++)
  {
    int xp2, yp2;
    xp2 = x + nx[i];
    yp2 = y + ny[i];

    if (par.periodic_boundaries)
    {
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (xp2 <= 0)
        xp2 = sizex - 2 + xp2;
      if (yp2 <= 0)
        yp2 = sizey - 2 + yp2;
      if (xp2 >= sizex - 1)
        xp2 = xp2 - sizex + 2;
      if (yp2 >= sizey - 1)
        yp2 = yp2 - sizey + 2;
      neighsite = sigma[xp2][yp2];
    }
    else
    {
      if (xp2 <= 0 || yp2 <= 0 || xp2 >= sizex - 1 || yp2 >= sizey - 1)
        neighsite = -1;
      else
        neighsite = sigma[xp2][yp2];
    }
    if (par.micropatternmask != string("None"))
    {
      if (!mask[xp2][yp2])
        neighsite = -1;
    }
    if (neighsite == -1)
    {
      // border
      DH += (sxyp == 0 ? 0 : par.border_energy) - (sxy == 0 ? 0 : par.border_energy);
    }
    else
    {
      DH += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
    }
  }
  // lambda is determined by chemical 0
  // cerr << "[" << lambda << "]";
  if (sxyp == MEDIUM)
  {
    DH += (int)(par.lambda * (1. - 2. *
                                       (double)((*cell)[sxy].Area() - (*cell)[sxy].TargetArea())));
  }
  else if (sxy == MEDIUM)
  {
    DH += (int)((par.lambda * (1. + 2. *
                                        (double)((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()))));
  }
  else
    DH += (int)((par.lambda * (2. + 2. * (double)((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea() - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea()))));

  /* Chemotaxis */ /*
  if (PDEfield && (par.vecadherinknockout || (sxyp==0 || sxy==0))) {
    // copying from (xp, yp) into (x,y)
    // If par.extensiononly == true, apply CompuCell's method, i.e.
    // only chemotactic extensions contribute to energy change
    if (!( par.extensiononly && sxyp==0)) {
      int DDH=(int)(par.chemotaxis*(sat(PDEfield->PDEVARS(0,x,y))-sat(PDEfield->PDEVARS(0,xp,yp))));
      DH-=DDH;
    }
  }*/

  const double lambda2 = par.lambda2;
  /* Length constraint */
  // sp is expanding cell, s is retracting cell
  if (sxyp == MEDIUM)
  {
    DH -= (int)(lambda2 * (DSQR((*cell)[sxy].Length() - (*cell)[sxy].TargetLength()) - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x, y) -
                                                                                            (*cell)[sxy].TargetLength())));
  }
  else if (sxy == MEDIUM)
  {
    DH -= (int)(lambda2 * (DSQR((*cell)[sxyp].Length() - (*cell)[sxyp].TargetLength()) - DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x, y) - (*cell)[sxyp].TargetLength())));
  }
  else
  {
    DH -= (int)(lambda2 * ((DSQR((*cell)[sxyp].Length() - (*cell)[sxyp].TargetLength()) - DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x, y) - (*cell)[sxyp].TargetLength())) +
                           (DSQR((*cell)[sxy].Length() - (*cell)[sxy].TargetLength()) - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x, y) -
                                                                                             (*cell)[sxy].TargetLength()))));
  }


  if (par.conv_ext && x > 305 && x < 405 && y > 95 && y < 315){
    vector<Cell>::iterator c = cell->begin();
	  ++c;
    double pulling_addition = 0;
    bool cell1_increase; // Indicate whether the first / second cells increase / decrease from this copy attempt
    bool cell2_increase;
    bool cell1_decrease;
    bool cell2_decrease;
    
    if (sxyp == MEDIUM) // source pixel belongs to medium
    {
      for (; c != cell->end(); c++)
      {
        if (c->Links[(*cell)[sxy].sigma] && c->area != 0 && (*cell)[sxy].area != 0) // link from c to cell with sxy
        {																			// cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = false;
          cell2_increase = false;
          cell1_decrease = true;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
        if ((*cell)[sxy].Links[c->sigma] && c->area != 0 && (*cell)[sxy].area != 0) // link from cell with sxy to c
        {																			// cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = false;
          cell2_increase = false;
          cell1_decrease = true;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
      }
    }
    else if (sxy == MEDIUM) // target pixel belongs to medium
    {
      for (; c != cell->end(); c++)
      {
        if (c->Links[(*cell)[sxyp].sigma] && c->area != 0 && (*cell)[sxyp].area != 0) // link from c to cell with sxyp
        {																			  // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = true;
          cell2_increase = false;
          cell1_decrease = false;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxyp], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
        if ((*cell)[sxyp].Links[c->sigma] && c->area != 0 && (*cell)[sxyp].area != 0) // link from cell with sxyp to c
        {																			  // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = true;
          cell2_increase = false;
          cell1_decrease = false;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxyp], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
      }
    }
    else // neither pixel belongs to medium
    {
      for (; c != cell->end(); c++)
      {
        if (c->Links[(*cell)[sxy].sigma] && c->area != 0 && (*cell)[sxy].area != 0 && c->sigma != (*cell)[sxyp].sigma) // link from c to cell with sxy
        {																											   // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = false;
          cell2_increase = false;
          cell1_decrease = true;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }

        if ((*cell)[sxy].Links[c->sigma] && c->area != 0 && (*cell)[sxy].area != 0 && c->sigma != (*cell)[sxyp].sigma) // link from cell with sxy to c
        {																											   // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = false;
          cell2_increase = false;
          cell1_decrease = true;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }

        if (c->Links[(*cell)[sxyp].sigma] && c->area != 0 && (*cell)[sxyp].area != 0 && c->sigma != (*cell)[sxy].sigma) // link from c to cell with sxyp
        {																												// cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = true;
          cell2_increase = false;
          cell1_decrease = false;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxyp], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
        if ((*cell)[sxyp].Links[c->sigma] && c->area != 0 && (*cell)[sxyp].area != 0 && c->sigma != (*cell)[sxy].sigma) // link from cell with sxyp to c
        {																												// cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
          cell1_increase = true;
          cell2_increase = false;
          cell1_decrease = false;
          cell2_decrease = false;
          pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxyp], &(*c), x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
        }
      }
      if ((*cell)[sxy].Links[(*cell)[sxyp].sigma]) // link from cell with sxy to cell with sxy
      {											 // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
        cell1_increase = true;
        cell2_increase = false;
        cell1_decrease = false;
        cell2_decrease = true;
        pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*cell)[sxyp], x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
      }
      if ((*cell)[sxyp].Links[(*cell)[sxy].sigma]) // link from cell with sxyp to cell with sxy
      {											 // cells do not always get apoptosed correctly. Do not try to pull non-exsistent cells.
        cell1_increase = true;
        cell2_increase = false;
        cell1_decrease = false;
        cell2_decrease = true;
        pulling_addition -= (int)(par.lambda_force * LengthDifference(&(*cell)[sxy], &(*cell)[sxyp], x, y, cell1_increase, cell2_increase, cell1_decrease, cell2_decrease));
      }
    }
    DH += int(pulling_addition);
  }
  return DH;
}

double CellularPotts::LengthDifference(Cell *cell1, Cell *cell2, int x, int y, bool cell1_increase, bool cell2_increase, bool cell1_decrease, bool cell2_decrease)
{

	if ((cell1_increase && cell2_increase) || (cell1_decrease && cell2_decrease) || (cell1_decrease && cell1_increase) || (cell2_decrease && cell2_increase))
	{
		cout << "Something went wrong in function LengthDifference, impossible situation" << endl;
		exit(1);
	}

	double length_before; // length before a copy attempt
	double length_after;  // length after a copy attempt

	double x1_before, x2_before, y1_before, y2_before, x1_after, x2_after, y1_after, y2_after;

	x1_before = cell1->getCenterX();
	x2_before = cell2->getCenterX();
	y1_before = cell1->getCenterY();
	y2_before = cell2->getCenterY();

	length_before = sqrt(pow(x1_before - x2_before, 2) + pow(y1_before - y2_before, 2));

	if (cell1_increase)
	{
		x1_after = cell1->getCenterXIfXYWereAdded(x, y);
		y1_after = cell1->getCenterYIfXYWereAdded(x, y);
	}
	else if (cell1_decrease)
	{
		x1_after = cell1->getCenterXIfXYWereRemoved(x, y);
		y1_after = cell1->getCenterYIfXYWereRemoved(x, y);
	}
	else
	{
		x1_after = x1_before;
		y1_after = y1_before;
	}

	if (cell2_increase)
	{
		x2_after = cell2->getCenterXIfXYWereAdded(x, y);
		y2_after = cell2->getCenterYIfXYWereAdded(x, y);
	}
	else if (cell2_decrease)
	{
		x2_after = cell2->getCenterXIfXYWereRemoved(x, y);
		y2_after = cell2->getCenterYIfXYWereRemoved(x, y);
	}
	else
	{
		x2_after = x2_before;
		y2_after = y2_before;
	}

	length_after = sqrt(pow(x1_after - x2_after, 2) + pow(y1_after - y2_after, 2));

	return length_before - length_after;
}

int CellularPotts::Act_AmoebaeMove(PDE *PDEfield)
{
  int loop, p;
  thetime++;
  int SumDH = 0;
  if (frozen)
    return 0;
  loop = (sizex - 2) * (sizey - 2);
  for (int i = 0; i < loop; i++)
  {
    // take a random site
    int xy = (int)(RANDOM() * (sizex - 2) * (sizey - 2));
    int x = xy % (sizex - 2) + 1;
    int y = xy / (sizex - 2) + 1;
    // take a random neighbour
    int xyp = (int)(n_nb * RANDOM() + 1);
    int xp = nx[xyp] + x;
    int yp = ny[xyp] + y;
    int k = sigma[x][y];
    int kp;
    if (par.periodic_boundaries)
    {
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (xp <= 0)
        xp = sizex - 2 + xp;
      if (yp <= 0)
        yp = sizey - 2 + yp;
      if (xp >= sizex - 1)
        xp = xp - sizex + 2;
      if (yp >= sizey - 1)
        yp = yp - sizey + 2;
      kp = sigma[xp][yp];
    }
    else
    {
      if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
      {
        kp = -1;
      }
      else
      {
        kp = sigma[xp][yp];
      }
    }
    // test for border state (relevant only if we do not use
    // periodic boundaries)
    if (kp != -1)
    {
      // Don't even think of copying the special border state into you!
      if (k >= 0 && k != kp)
      {
        if (par.cluster_connectivity == false || ConnectivityPreservedPCluster(x, y))
        {
          /* Try to copy if sites do not belong to the same cell */
          // connectivity dissipation:
          int H_diss = 0;
          if (!ConnectivityPreservedP(x, y))
            H_diss = par.conn_diss;
          int D_H = Act_DeltaH(x, y, xp, yp, PDEfield);
          if ((p = CopyvProb(D_H, H_diss, false)) > 0)
          {
            ConvertSpin(x, y, xp, yp);
            SumDH += D_H;
            if (par.lambda_Act > 0)
            {
              // Update actin field
              if (sigma[x][y] > 0)
              {
                actPixels[{x, y}] = par.max_Act;
                std::unordered_set<std::array<int, 2>>::const_iterator it = (alivePixels.find({x, y}));
                if (it == alivePixels.end())
                {
                  alivePixels.insert({x, y});
                }
              }
              else
              {
                std::unordered_set<std::array<int, 2>>::const_iterator it = (alivePixels.find({x, y}));
                if (it != alivePixels.end())
                {
                  alivePixels.erase({x, y});
                }
                std::unordered_map<std::array<int, 2>, int>::const_iterator ap = (actPixels.find({x, y}));
                if (ap != actPixels.end())
                {
                  actPixels.erase({x, y});
                }
              }
            }
            // Update adhesive areas
            if (kp == 0)
            {
              getCell(k).DecrementAdhesiveArea(GetMatrixLevel(x, y));
            }
            else if (k == 0)
            {
              // getCell(kp).IncrementAdhesiveArea(1);
            }
            else
            {
              getCell(k).DecrementAdhesiveArea(GetMatrixLevel(x, y));
              // getCell(kp).IncrementAdhesiveArea(1);
            }
            if (par.lambda_matrix > 0)
            {
              // Update matrix interaction field
              if (sigma[x][y] > 0)
              {
                // matrixPixels[{x,y}]=0;
                matrix[x][y] = 0;
              }
              else
              {
                // matrixPixels.erase({x,y});
                matrix[x][y] = 0;
              }
            }
          }
        }
      }
    }
  }
  return SumDH;
}

int CellularPotts::Act_DeltaH(int x, int y, int xp, int yp, PDE *PDEfield)
{
  int DH = 0;
  int i, sxy, sxyp;
  int neighsite;

  /* Compute energydifference *IF* the copying were to occur */
  int DH_adhesive_energy = 0;
  sxy = sigma[x][y];
  sxyp = sigma[xp][yp];
  if (sxyp <= -1)
  {
    sxyp = 0;
  } // allow for medium to copy in from box border or pillars

  /* DH due to cell adhesion */
  // also compute changes in neighbours for alignment with newneighbours
  std::vector<int> xy_neighbour_changes(cell->size(), 0);
  std::vector<int> xyp_neighbour_changes(cell->size(), 0);
  for (i = 1; i <= n_nb; i++)
  {
    int xp2, yp2;
    xp2 = x + nx[i];
    yp2 = y + ny[i];
    if (par.periodic_boundaries)
    {

      // since we are asynchronic, we cannot just copy the borders once
      // every MCS

      if (xp2 <= 0)
        xp2 = sizex - 2 + xp2;
      if (yp2 <= 0)
        yp2 = sizey - 2 + yp2;
      if (xp2 >= sizex - 1)
        xp2 = xp2 - sizex + 2;
      if (yp2 >= sizey - 1)
        yp2 = yp2 - sizey + 2;

      neighsite = sigma[xp2][yp2];
    }
    else
    {

      if (xp2 <= 0 || yp2 <= 0 || xp2 >= sizex - 1 || yp2 >= sizey - 1)
        neighsite = -1;
      else
        neighsite = sigma[xp2][yp2];
    }
    if (neighsite == -1)
    { // border
      DH_adhesive_energy += (sxyp == 0 ? 0 : par.border_energy) -
                            (sxy == 0 ? 0 : par.border_energy);
    }
    else
    {
      DH_adhesive_energy += (*cell)[sxyp].EnergyDifference((*cell)[neighsite]) - (*cell)[sxy].EnergyDifference((*cell)[neighsite]);
      if ((i <= 4) | par.extended_neighbour_border)
      {
        if (neighsite != sxy)
          xy_neighbour_changes[neighsite] -= 1;
        if (neighsite != sxyp)
          xyp_neighbour_changes[neighsite] += 1;
      }
    }
  }

  DH += DH_adhesive_energy;

  // lambda is determined by chemical 0
  int DH_area = 0;
  if (par.area_constraint_type == 0)
  {
    // cerr << "[" << lambda << "]";
    if (sxyp == MEDIUM)
    {
      DH_area += (int)(par.lambda * (1. - 2. *
                                              (double)((*cell)[sxy].Area() - (*cell)[sxy].TargetArea())));
    }
    else if (sxy == MEDIUM)
    {
      DH_area += (int)((par.lambda * (1. + 2. *
                                               (double)((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea()))));
    }
    else
      DH_area += (int)(par.lambda * (2. + 2. * (double)((*cell)[sxyp].Area() - (*cell)[sxyp].TargetArea() - (*cell)[sxy].Area() + (*cell)[sxy].TargetArea())));
  }

  DH += DH_area;

  //! Perimeter constraint, only available when par.area_constraint_type==0, in other words,
  // when the target area constraint is in place
  int DH_perimeter = 0;
  if (par.area_constraint_type == 0)
  {
    if (sxyp == MEDIUM)
    {

      DH_perimeter -= par.lambda_perimeter * (DSQR((*cell)[sxy].Perimeter() - (*cell)[sxy].TargetPerimeter()) - DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y) - (*cell)[sxy].TargetPerimeter()));
    }
    else if (sxy == MEDIUM)
    {

      DH_perimeter -= par.lambda_perimeter * (DSQR((*cell)[sxyp].Perimeter() - (*cell)[sxyp].TargetPerimeter()) - DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y) - (*cell)[sxyp].TargetPerimeter()));
    }
    // they're both cells
    else
    {

      DH_perimeter -= par.lambda_perimeter * ((DSQR((*cell)[sxyp].Perimeter() - (*cell)[sxyp].TargetPerimeter()) - DSQR(GetNewPerimeterIfXYWereAdded(sxyp, x, y) - (*cell)[sxyp].TargetPerimeter())));

      DH_perimeter -= par.lambda_perimeter * (DSQR((*cell)[sxy].Perimeter() - (*cell)[sxy].TargetPerimeter()) - DSQR(GetNewPerimeterIfXYWereRemoved(sxy, x, y) - (*cell)[sxy].TargetPerimeter()));
    }
  }
  DH += DH_perimeter;
  /* Chemotaxis */
  int DDH = 0;
  if (PDEfield && (par.vecadherinknockout || (sxyp == 0 || sxy == 0)))
  {

    // copying from (xp, yp) into (x,y)
    // If par.extensiononly == true, apply CompuCell's method, i.e.
    // only chemotactic extensions contribute to energy change
    if (!(par.extensiononly && sxyp == 0))
    {
      DDH = (int)(par.chemotaxis * (sat(PDEfield->PDEVARS(0, x, y)) - sat(PDEfield->PDEVARS(0, xp, yp))));

      DH -= DDH;
    }
  }

  const double lambda2 = par.lambda2;

  /* Length constraint */
  // sp is expanding cell, s is retracting cell

  int DH_length = 0;
  if (sxyp == MEDIUM)
  {
    DH_length -= (int)(lambda2 * (DSQR((*cell)[sxy].Length() - (*cell)[sxy].TargetLength()) - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x, y) -
                                                                                                   (*cell)[sxy].TargetLength())));
  }
  else if (sxy == MEDIUM)
  {
    DH_length -= (int)(lambda2 * (DSQR((*cell)[sxyp].Length() - (*cell)[sxyp].TargetLength()) - DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x, y) - (*cell)[sxyp].TargetLength())));
  }
  else
  {
    DH_length -= (int)(lambda2 * ((DSQR((*cell)[sxyp].Length() - (*cell)[sxyp].TargetLength()) - DSQR((*cell)[sxyp].GetNewLengthIfXYWereAdded(x, y) - (*cell)[sxyp].TargetLength())) +
                                  (DSQR((*cell)[sxy].Length() - (*cell)[sxy].TargetLength()) - DSQR((*cell)[sxy].GetNewLengthIfXYWereRemoved(x, y) -
                                                                                                    (*cell)[sxy].TargetLength()))));
  }
  DH += DH_length;

  /************************The Act model****************/
  // let the cell extend with
  int DH_act = 0;
  if (par.lambda_Act)
  {
    double Act_expanding = 1, Act_retracting = 1;
    int nxp = 0, nret = 0;

    for (int i1 = -1; i1 <= 1; i1++)
      for (int i2 = -1; i2 <= 1; i2++)
      {

        if (sigma[xp + i1][yp + i2] >= 0 && sigma[xp + i1][yp + i2] == sigma[xp][yp])
        {
          Act_expanding *= GetActLevel(xp + i1, yp + i2);
          nxp++;
        }

        if (sigma[x + i1][y + i2] >= 0 && sigma[x + i1][y + i2] == sigma[x][y])
        {
          Act_retracting *= GetActLevel(x + i1, y + i2);
          nret++;
        }
      }

    // apply the smoothing
    // Act_expanding*= pow( PDEfield->Sigma(2,xp,yp), w-1);
    // Act_retracting*= pow(PDEfield->Sigma(2,x,y), w-1);
    // nxp += w-1;
    // nret += w-1;

    Act_expanding = pow(Act_expanding, 1. / nxp);
    Act_retracting = pow(Act_retracting, 1. / nret);

    // Act model activation dependent on total adhesion area of the cell
    // If adhesion area exceeds threshold, Act model is fully functional,
    // otherwise, it starts at base*lambda_Act and for increasing adhesion areas
    // it increases linearly to lambda_Act.
    double threshold = par.threshold;
    double base = par.start_level;
    double strength;
    if ((*cell)[sxyp].sigma > 0)
    {
      double adhesion_fraction = (double)(*cell)[sxyp].GetAdhesiveArea() / (double)(*cell)[sxyp].area;
      if (adhesion_fraction >= threshold)
      {
        strength = 1;
      }
      else
      {
        strength = base + ((1 - base) / threshold) * adhesion_fraction;
      }
      DH_act -= (par.lambda_Act * strength) / par.max_Act * Act_expanding;
    }

    if ((*cell)[sxy].sigma > 0)
    {
      double adhesion_fraction = (double)(*cell)[sxy].GetAdhesiveArea() / (double)(*cell)[sxy].area;
      if (adhesion_fraction >= threshold)
      {
        strength = 1;
      }
      else
      {
        strength = base + ((1 - base) / threshold) * adhesion_fraction;
      }
      DH_act += (par.lambda_Act * strength) / par.max_Act * Act_retracting;
    }
  }
  DH += DH_act;

  /****** Matrix interaction retraction yield energy ****/
  // Retractiong of lattice sites that contain an adhesion is penalized with
  // lambda_matrix
  int DH_matrix_interaction = 0;
  if (sxyp == MEDIUM && par.lambda_matrix)
  { // should be done for all retractions, I assume.
    DH_matrix_interaction += par.lambda_matrix * (GetMatrixLevel(x, y));
  }
  DH += DH_matrix_interaction;
  // std::cout << "DH_matrix: " << DH_matrix_interaction << std::endl;
  return DH;
}

bool CellularPotts::Probability(int DH)
{
  if (DH > BOLTZMANN - 1)
    return false;
  else if (DH < 0 || RANDOM() < copyprob[DH])
    return true;
  return false;
}

void CellularPotts::ConvertSpin(int x, int y, int xp, int yp)
{
  int tmpcell;
  if ((tmpcell = sigma[x][y]))
  { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x, y);
    (*cell)[tmpcell].SetPerimeter(GetNewPerimeterIfXYWereRemoved(tmpcell, x, y));
    if (!(*cell)[tmpcell].Area())
    {
      (*cell)[tmpcell].Apoptose();
    }
  }

  if ((tmpcell = sigma[xp][yp]))
  { // if tmpcell is not MEDIUM
    (*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x, y);
    (*cell)[tmpcell].SetPerimeter(GetNewPerimeterIfXYWereAdded(tmpcell, x, y));
  }
  sigma[x][y] = sigma[xp][yp];
  tau[x][y] = tau[xp][yp];
}

void CellularPotts::ExchangeSpin(int x, int y, int xp, int yp)
{
  int tmpcell;
  if ((tmpcell = sigma[x][y]))
  { // if tmpcell is not MEDIUM
    //(*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x, y);
  }

  if ((tmpcell = sigma[xp][yp]))
  { // if tmpcell is not MEDIUM
    //(*cell)[tmpcell].DecrementArea();
    (*cell)[tmpcell].RemoveSiteFromMoments(x, y);
  }

  if ((tmpcell = sigma[x][y]))
  { // if tmpcell is not MEDIUM
    //(*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x, y);
  }

  if ((tmpcell = sigma[xp][yp]))
  { // if tmpcell is not MEDIUM
    //(*cell)[tmpcell].IncrementArea();
    (*cell)[tmpcell].AddSiteToMoments(x, y);
  }

  // Exchange spins
  tmpcell = sigma[x][y];
  sigma[x][y] = sigma[xp][yp];
  sigma[xp][yp] = tmpcell;
  tmpcell = tau[x][y];
  tau[x][y] = tau[xp][yp];
  tau[xp][yp] = tmpcell;
}

/** PUBLIC **/
int CellularPotts::CopyvProb(int DH, double stiff, bool anneal)
{
  if (stiff == par.conn_diss && par.conn_diss > 0){
    return 0;
  }
  double dd;
  int s;
  s = (int)stiff;
  if (DH <= -s)
    return 2;
  if (anneal)
    return 0;
  // if DH becomes extremely large, calculate probability on-the-fly
  if (DH + s > BOLTZMANN - 1)
    dd = exp(-((double)(DH + s) / par.T));
  else
    dd = copyprob[DH + s];

  if (RANDOM() < dd)
    return 1;
  else
    return 0;
}

void CellularPotts::CopyProb(double T)
{
  int i;
  for (i = 0; i < BOLTZMANN; i++)
    copyprob[i] = exp(-((double)(i) / T));
}

void CellularPotts::FreezeAmoebae(void)
{
  if (frozen)
    frozen = FALSE;
  else
    frozen = TRUE;
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::AmoebaeMove(PDE *PDEfield, bool anneal)
{
  int p;
  float loop;
  thetime++;
  int SumDH = 0;

  int positionedge;
  int targetedge;
  int targetsite;
  int targetneighbour;
  int x, y;
  int xp, yp;

  int H_diss;
  int D_H;

  int edgeadjusting;
  int xn, yn; // neighbour cells

  for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
    
  

  if (frozen)
    return 0;

  loop = static_cast<float>(sizeedgelist) / static_cast<float>(n_nb);
  for (int i = 0; i < loop; i++){
    // take a random entry of the edgelist
    positionedge = (int)(RANDOM()*sizeedgelist); 
    // find the corresponding edge
    targetedge = orderedgelist[positionedge];
    // find the lattice site corresponding to this edge
    targetsite = targetedge / n_nb;
     // find the neighbour corresponding to this edge
    targetneighbour = (targetedge % n_nb)+1;
    
    // find the x and y coordinate corresponding to the target site
    x = targetsite%(sizex-2)+1;
    y = targetsite/(sizex-2)+1; 
    
    // find the neighbouring site corresponding to this edge
    xp = nx[targetneighbour]+x;
    yp = ny[targetneighbour]+y;
    if (par.periodic_boundaries) {
       // since we are asynchronic, we cannot just copy the borders once 
       // every MCS
      if (xp<=0)
        xp=sizex-2+xp;
      if (yp<=0)
        yp=sizey-2+yp;
      if (xp>=sizex-1)
        xp=xp-sizex+2;
      if (yp>=sizey-1)
        yp=yp-sizey+2;
    }


    D_H = DeltaH(x, y, xp, yp, PDEfield);
    if ((p = CopyvProb(D_H, H_diss, anneal)) > 0 && LocalConnectedness(x,y,sigma[x][y]) && LocalConnectedness(x,y,sigma[xp][yp]))
    {
      ConvertSpin(x, y, xp, yp); // sigma(x,y) will get the same value as sigma(xp,yp)
      if (par.n_chem > 0)
        CopyPDEvars(x, y, xp, yp, PDEfield);
      for (int j = 1; j <= n_nb; j++)
      {
        xn = nx[j] + x; //Update the edgelist for all neighbours of x,y
        yn = ny[j] + y;
        edgeadjusting = targetsite * n_nb + j - 1;

        if (par.periodic_boundaries)
        {
          // since we are asynchronic, we cannot just copy the borders once
          // every MCS
          if (xn <= 0)
            xn = sizex - 2 + xn;
          if (yn <= 0)
            yn = sizey - 2 + yn;
          if (xn >= sizex - 1)
            xn = xn - sizex + 2;
          if (yn >= sizey - 1)
            yn = yn - sizey + 2;
        }
        if (par.micropatternmask != string("None"))
        {
          if (xn > 0 && yn > 0 && xn < sizex - 1 && yn < sizey - 1)
          { // if the neighbour site is within the lattice
            if (mask[xn][yn])
            { // and if it is within the micropattern
              if (edgelist[edgeadjusting] == -1 && sigma[xn][yn] != sigma[x][y])
                if(!par.second_layer || tau[x][y] == tau[xn][yn]) //We must either have no second layer, or if we do, both cell types need to be of the same type
                { // if we should add the edge to the edgelist, add it
                  AddEdgeToEdgelist(edgeadjusting);
                  loop += (double)2 / n_nb;
                  numberofedges[x][y]++;
                  numberofedges[xn][yn]++;
                  if ((tau[xn][yn] == 1 && tau[x][y] == 2) || (tau[xn][yn] == 2 && tau[x][y] == 1)){
                    couplingcoefficient[x][y] = par.couplingAtrialPM;
                    couplingcoefficient[xn][yn] = par.couplingAtrialPM;
                  }
                  else if (tau[xn][yn] == 0 || tau[x][y] == 0){
                    couplingcoefficient[x][y] = par.couplingmedium;
                    couplingcoefficient[xn][yn] = par.couplingmedium;
                  }
                  else if (tau[xn][yn] == 1){
                    couplingcoefficient[x][y] = par.couplingAtrialAtrial;
                    if (couplingcoefficient[xn][yn] != par.couplingAtrialPM)
                      couplingcoefficient[xn][yn] = par.couplingAtrialAtrial;
                  }
                  else if (tau[xn][yn] == 2){
                    couplingcoefficient[x][y] = par.couplingPMPM;
                    if (couplingcoefficient[xn][yn] != par.couplingAtrialPM)
                      couplingcoefficient[xn][yn] = par.couplingPMPM;
                  }
                  
                }
              if (edgelist[edgeadjusting] != -1 && ((sigma[xn][yn] == sigma[x][y]) ||  (par.second_layer && tau[xn][yn] != tau[x][y])))
              { // if the sites have the same cellnumber and they have an edge, remove it
                //also remove the edge if we do have a second mask layer and the cells are of different celltype
                RemoveEdgeFromEdgelist(edgeadjusting);
                loop -= (double)2 / n_nb;
                numberofedges[x][y]--;
                numberofedges[xn][yn]--;
                if (numberofedges[x][y] == 0){
                  if (sigma[x][y] == 0)
                    couplingcoefficient[x][y] = par.couplingmedium;
                  else 
                    couplingcoefficient[x][y] = par.couplingcell;
                }
                if (numberofedges[xn][yn] == 0){
                  if (sigma[x][y] == 0)
                    couplingcoefficient[xn][yn] = par.couplingmedium;
                  else
                    couplingcoefficient[xn][yn] = par.couplingcell;
                }
              }
            }
          }
        }
        else
        {
          if (xn > 0 && yn > 0 && xn < sizex - 1 && yn < sizey - 1)
          { // if the neighbour site is within the lattice
            if (edgelist[edgeadjusting] == -1 && sigma[xn][yn] != sigma[x][y])
            { 
              //if there should be an edge between (x,y) and (xn,yn) and it is not there yet, add it
              AddEdgeToEdgelist(edgeadjusting);
              //adjust loop because two edges were added
              loop += 2.0 / n_nb;
              numberofedges[x][y]++;
              numberofedges[xn][yn]++;
              if (couplingcoeffcient_allocated){
                if ((tau[xn][yn] == 1 && tau[x][y] == 2) || (tau[xn][yn] == 2 && tau[x][y] == 1)){
                  couplingcoefficient[x][y] = par.couplingAtrialPM;
                  couplingcoefficient[xn][yn] = par.couplingAtrialPM;
                }
                else if (tau[xn][yn] == 0 || tau[x][y] == 0){
                  couplingcoefficient[x][y] = par.couplingmedium;
                  couplingcoefficient[xn][yn] = par.couplingmedium;
                }
                else if (tau[xn][yn] == 1){
                  couplingcoefficient[x][y] = par.couplingAtrialAtrial;
                  if (couplingcoefficient[xn][yn] != par.couplingAtrialPM)
                    couplingcoefficient[xn][yn] = par.couplingAtrialAtrial;
                }
                else if (tau[xn][yn] == 2){
                  couplingcoefficient[x][y] = par.couplingPMPM;
                  if (couplingcoefficient[xn][yn] != par.couplingAtrialPM)
                    couplingcoefficient[xn][yn] = par.couplingPMPM;
                }
              }    
            }
            if (edgelist[edgeadjusting] != -1 && sigma[xn][yn] == sigma[x][y])
            { 
              //if there should be no edge between (x,y) and (xn,yn), but there is an edge remove it 
              RemoveEdgeFromEdgelist(edgeadjusting);
              //adjust loop because two edges were removed
              loop -= 2.0 / n_nb;
              numberofedges[x][y]--;
              numberofedges[xn][yn]--;
              if (couplingcoeffcient_allocated){
                if (numberofedges[x][y] == 0)
                {
                  if (sigma[x][y] == 0)
                    couplingcoefficient[x][y] = par.couplingmedium;
                  else if (sigma[x][y] == 0)
                    couplingcoefficient[x][y] = par.couplingcell;
                }
                if (numberofedges[xn][yn] == 0)
                {
                  if (sigma[x][y] == 0)
                    couplingcoefficient[xn][yn] = par.couplingmedium;
                  else if (sigma[xn][yn] == 0)
                    couplingcoefficient[xn][yn] = par.couplingcell;
                }
              }
            }
          }
        }
      }
      SumDH += D_H;
    }
  }
  return SumDH;
}

void CellularPotts::AddEdgeToEdgelist(int edge) {//add an edge to the end of edgelist
  int counteredge = CounterEdge(edge);  

  //assign a unique integer to position 'edge' in the edgelist
  edgelist[edge] = sizeedgelist;
  //assign a unique integer at the end of orderedgelist, maintaining the bijection between the lists
  orderedgelist[sizeedgelist] = edge;
  //Increase the size of the array
  sizeedgelist++;

  //Repeat for the counteredge
  edgelist[counteredge] = sizeedgelist;
  orderedgelist[sizeedgelist] = counteredge;
  sizeedgelist++;
}


void CellularPotts::RemoveEdgeFromEdgelist(int edge) { //remove an edge from the edgelist
  int counteredge = CounterEdge(edge);  

  if(edgelist[edge] != sizeedgelist-1){ //if edge is not the last edge in orderedgelist
    // move the edge in the last position to the position of the edge that must be deleted
    orderedgelist[edgelist[edge]] = orderedgelist[sizeedgelist-1];
    edgelist[orderedgelist[sizeedgelist-1]] = edgelist[edge];
  }
  // remove the edge from the edgelist
  edgelist[edge] = -1;
  // free the last position of orderedgelist
  orderedgelist[sizeedgelist-1] = -1;
  //decrease the size of the edgelist
  sizeedgelist--;
  
  //Repeat for counteredge
  if(edgelist[counteredge] != sizeedgelist-1){ 
    orderedgelist[edgelist[counteredge]] = orderedgelist[sizeedgelist-1];
    edgelist[orderedgelist[sizeedgelist-1]] = edgelist[counteredge];
  }
  edgelist[counteredge] = -1;
  orderedgelist[sizeedgelist - 1] = -1;
  sizeedgelist--;
}

int CellularPotts::CounterEdge(int edge){
  // For an edge from (x,y) to (xn,yn), this function returns the edge from (xn,yn) to (x,y)
  
  // find the corresponding lattice site and neighbour for the edge.
  int which_site = edge / n_nb;
  int which_neighbour = edge % n_nb + 1;
  int counterneighbour = 0;

  // find the x and y coordinate corresponding to the lattice site
  int x = which_site%(sizex-2)+1;
  int y = which_site/(sizex-2)+1; 

  // find the x and y coordinate corresponding at the other end of the edge
  int xp = nx[which_neighbour]+x; 
  int yp = ny[which_neighbour]+y; 

  // correct for periodic boundaries of necessary
  if (par.periodic_boundaries) {  
    // since we are asynchronic, we cannot just copy the borders once 
    // every MCS
    if (xp <= 0)
      xp = sizex - 2 + xp;
    if (yp <= 0)
      yp = sizey - 2 + yp;
    if (xp >= sizex - 1)
      xp = xp - sizex + 2;
    if (yp >= sizey - 1)
      yp = yp - sizey + 2;
  }




  // lattice site corresponding to other site of the edge
  int neighbourlocation = xp-1 + (yp-1)*(par.sizex-2);
  
  // find the neighbour pointing the other direction
  const int counterneighbourlist[20] = {3, 4, 1, 2, 7, 8, 5, 6, 11, 12, 9, 10, 17, 18, 19, 20, 13, 14, 15, 16};
  counterneighbour = counterneighbourlist[ which_neighbour - 1 ];
  // compute the final counteredge
  int counteredge = neighbourlocation * n_nb + counterneighbour-1;
  return counteredge;
}

void CellularPotts::CopyPDEvars(int x, int y, int xp, int yp, PDE *PDEFIELD){
  PDEFIELD_TYPE* xy = PDEFIELD->PDE_pointer(x,y);
  PDEFIELD_TYPE* xyp = PDEFIELD->PDE_pointer(xp,yp);
  for (int layer = 0; layer < par.n_chem; layer++)
    xy[layer*sizex*sizey] = xyp[layer*sizex*sizey];
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::KawasakiMove(PDE *PDEfield)
{
  int loop, p;
  // int updated=0;
  thetime++;
  int SumDH = 0;

  if (frozen)
    return 0;

  loop = (sizex - 2) * (sizey - 2);
  for (int i = 0; i < loop; i++)
  {
    // take a random site
    int xy = (int)(RANDOM() * (sizex - 2) * (sizey - 2));
    int x = xy % (sizex - 2) + 1;
    int y = xy / (sizex - 2) + 1;

    // take a random neighbour
    int xyp = (int)(n_nb * RANDOM() + 1);
    int xp = nx[xyp] + x;
    int yp = ny[xyp] + y;

    int k = sigma[x][y];

    int kp;
    if (par.periodic_boundaries)
    {
      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (xp <= 0)
        xp = sizex - 2 + xp;
      if (yp <= 0)
        yp = sizey - 2 + yp;
      if (xp >= sizex - 1)
        xp = xp - sizex + 2;
      if (yp >= sizey - 1)
        yp = yp - sizey + 2;
      kp = sigma[xp][yp];
    }
    else
    {
      if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
        kp = -1;
      else
        kp = sigma[xp][yp];
    }
    // test for border state (relevant only if we do not use
    // periodic boundaries)
    if (kp != -1)
    {
      // Don't even think of copying the special border state into you!
      if (k != kp)
      {
        /* Try to exchange sites if sites do not belong to the same cell */
        // connectivity dissipation:
        int H_diss = 0;
        // if (!ConnectivityPreservedP(x,y)) H_diss=par.conn_diss;
        int D_H = KawasakiDeltaH(x, y, xp, yp, PDEfield);
        if (D_H != 0 && (p = CopyvProb(D_H, H_diss, false)) > 0)
        {
          ExchangeSpin(x, y, xp, yp);
          SumDH += D_H;
        }
        // std::cerr << "[ " << D_H << ", p = " << p << " ]";
      }
    }
  }
  return SumDH;
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::IsingMove(PDE *PDEfield)
{
  int loop, p;
  // int updated=0;
  thetime++;
  int SumDH = 0;

  loop = (sizex - 2) * (sizey - 2);

  for (int i = 0; i < loop; i++)
  {

    // take a random site
    int xy = (int)(RANDOM() * (sizex - 2) * (sizey - 2));
    int x = xy % (sizex - 2) + 1;
    int y = xy / (sizex - 2) + 1;

    int D_H = IsingDeltaH(x, y, PDEfield);

    if (D_H != 0 && (p = CopyvProb(D_H, 0, false) > 0))
    {

      sigma[x][y] = sigma[x][y] == 0 ? 1 : 0;
      // tau[x][y]=(*cell)[sigma[x][y]].getTau();
      SumDH += D_H;
    }
    // std::cerr << "[ " << D_H << ", p = " << p << " ]";
  }

  return SumDH;
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::PottsMove(PDE *PDEfield)
{
  int loop, p;
  // int updated=0;
  thetime++;
  int SumDH = 0;
  loop = (sizex - 2) * (sizey - 2);

  for (int i = 0; i < loop; i++)
  {
    // take a random site
    int xy = (int)(RANDOM() * (sizex - 2) * (sizey - 2));
    int x = xy % (sizex - 2) + 1;
    int y = xy / (sizex - 2) + 1;

    int new_state = (int)(RANDOM() * par.n_init_cells);
    int D_H = PottsDeltaH(x, y, new_state);
    // cerr << "D_H = " << D_H << endl;
    if (D_H < 0 || (p = CopyvProb(D_H, 0, false) > 0))
    {
      sigma[x][y] = new_state;
      tau[x][y] = (*cell)[sigma[x][y]].getTau();
      // cerr << "[ " << x << ", " << y << "]";
      SumDH += D_H;
    }
    // std::cerr << "[ " << D_H << ", p = " << p << " ]";
  }
  return SumDH;
}

//! Monte Carlo Step. Returns summed energy change
int CellularPotts::PottsNeighborMove(PDE *PDEfield)
{
  int loop, p;
  // int updated=0;
  thetime++;
  int SumDH = 0;

  loop = (sizex - 2) * (sizey - 2);

  for (int i = 0; i < loop; i++)
  {
    // take a random site
    int xy = (int)(RANDOM() * (sizex - 2) * (sizey - 2));
    int x = xy % (sizex - 2) + 1;
    int y = xy / (sizex - 2) + 1;
    // take a random neighbour
    int xyp = (int)(n_nb * RANDOM() + 1);
    int xp = nx[xyp] + x;
    int yp = ny[xyp] + y;

    int kp;
    if (par.periodic_boundaries)
    {

      // since we are asynchronic, we cannot just copy the borders once
      // every MCS
      if (xp <= 0)
        xp = sizex - 2 + xp;
      if (yp <= 0)
        yp = sizey - 2 + yp;
      if (xp >= sizex - 1)
        xp = xp - sizex + 2;
      if (yp >= sizey - 1)
        yp = yp - sizey + 2;

      kp = sigma[xp][yp];
    }
    else
    {
      if (xp <= 0 || yp <= 0 || xp >= sizex - 1 || yp >= sizey - 1)
        kp = -1;
      else
        kp = sigma[xp][yp];
    }

    int D_H = PottsDeltaH(x, y, kp);
    // cerr << "D_H = " << D_H << endl;
    if (D_H < 0 || (p = CopyvProb(D_H, 0, false) > 0))
    {
      sigma[x][y] = kp;
      tau[x][y] = (*cell)[sigma[x][y]].getTau();
      // cerr << "[ " << x << ", " << y << "]";
      SumDH += D_H;
    }
    // std::cerr << "[ " << D_H << ", p = " << p << " ]";
  }
  return SumDH;
}

/** A simple method to plot all sigma's in window
    without the black lines */
void CellularPotts::PlotSigma(Graphics *g, int mag)
{
  for (int x = 1; x < sizex - 1; x++)
    for (int y = 1; y < sizey - 1; y++)
    {
      for (int xm = 0; xm < mag; xm++)
        for (int ym = 0; ym < mag; ym++)
          g->Point(sigma[x][y], mag * x + xm, mag * y + ym);
    }
}

/** Plot in black & white for the Ising model **/
void CellularPotts::PlotIsing(Graphics *g, int mag)
{
  for (int x = 1; x < sizex - 1; x++)
    for (int y = 1; y < sizey - 1; y++)
    {
      for (int xm = 0; xm < mag; xm++)
        for (int ym = 0; ym < mag; ym++)
          g->Point(sigma[x][y] == 0 ? 0 : 1, mag * x + xm, mag * y + ym);
    }
}

int **CellularPotts::SearchNandPlot(Graphics *g, bool get_neighbours)
{
  int i, j, q;
  int **neighbours = 0;

  /* Allocate neighbour matrix */
  if (get_neighbours)
  {
    neighbours = (int **)malloc((cell->size() + 1) * sizeof(int *));
    if (neighbours == NULL)
      MemoryWarning();

    neighbours[0] = (int *)malloc((cell->size() + 1) * (cell->size() + 1) * sizeof(int));
    if (neighbours[0] == NULL)
      MemoryWarning();

    for (i = 1; i < (int)cell->size() + 1; i++)
      neighbours[i] = neighbours[i - 1] + (cell->size() + 1);

    /* Clear this matrix */
    for (i = 0; i < ((int)cell->size() + 1) * ((int)cell->size() + 1); i++)
      neighbours[0][i] = EMPTY;
  }

  for (i = 0; i < sizex - 1; i++)
    for (j = 0; j < sizey - 1; j++)
    {
      int colour;
      if (sigma[i][j] <= 0)
      {
        colour = 0;
      }
      else
      {
        colour = (*cell)[sigma[i][j]].Colour();
        // colour = sigma[i][j];
      }

      if (g && sigma[i][j] > 0) /* if draw */
        g->Point(colour, i, j);

      if (sigma[i][j] != sigma[i + 1][j]) /* if cellborder */ /* etc. etc. */
      {
        if (g)
          g->Point(1, i + 1, j);
        if (get_neighbours)
        {
          if (sigma[i][j] > 0)
          {
            for (q = 0; q < (int)cell->size(); q++)
              if (neighbours[sigma[i][j]][q] == EMPTY)
              {
                neighbours[sigma[i][j]][q] = sigma[i + 1][j];
                break;
              }
              else if (neighbours[sigma[i][j]][q] == sigma[i + 1][j])
                break;
          }
          if (sigma[i + 1][j] > 0)
          {
            for (q = 0; q < (int)cell->size(); q++)
              if (neighbours[sigma[i + 1][j]][q] == EMPTY)
              {
                neighbours[sigma[i + 1][j]][q] = sigma[i][j];
                break;
              }
              else if (neighbours[sigma[i + 1][j]][q] == sigma[i][j])
                break;
          }
        }
      }
      else if (g && sigma[i][j] > 0)
        g->Point(colour, i + 1, j);

      if (sigma[i][j] != sigma[i][j + 1])
      {

        if (g)
          g->Point(1, i, j + 1);

        if (get_neighbours)
        {
          if (sigma[i][j] > 0)
          {
            for (q = 0; q < (int)cell->size(); q++)
              if (neighbours[sigma[i][j]][q] == EMPTY)
              {
                neighbours[sigma[i][j]][q] = sigma[i][j + 1];
                break;
              }
              else if (neighbours[sigma[i][j]][q] == sigma[i][j + 1])
                break;
          }

          if (sigma[i][j + 1] > 0)
          {

            for (q = 0; q < (int)cell->size(); q++)
              if (neighbours[sigma[i][j + 1]][q] == EMPTY)
              {
                neighbours[sigma[i][j + 1]][q] = sigma[i][j];
                break;
              }
              else if (neighbours[sigma[i][j + 1]][q] == sigma[i][j])
                break;
          }
        }
      }
      else if (g && sigma[i][j] > 0)
        g->Point(colour, i, j + 1);

      /* Cells that touch eachother's corners are NO neighbours */

      if (sigma[i][j] != sigma[i + 1][j + 1] || sigma[i + 1][j] != sigma[i][j + 1])
      {
        if (g)
          g->Point(1, i + 1, j + 1);
      }
      else if (g && sigma[i][j] > 0)
        g->Point(colour, i + 1, j + 1);
    }

  bool plotLinks = false;
	if (plotLinks)
	{
		vector<Cell>::iterator c = cell->begin();
		++c;
		for (; c != cell->end(); c++)
		{
			if (!(c->sigma % 10))
			{
				vector<Cell>::iterator z = cell->begin();
				++z;
				for (; z != cell->end(); z++)
				{
					if (c->Links[z->sigma])
						g->Line(2 * c->getCenterX(), 2 * c->getCenterY(), 2 * z->getCenterX(), 2 * z->getCenterY(), 7);
				}
			}
		}
	}


  if (get_neighbours)
    return neighbours;
  else
    return 0;
}

void CellularPotts::SearchNandPlotClear(Graphics *g)
{
  for (int i = 0; i < sizex - 1; i++)
  {
    for (int j = 0; j < sizey - 1; j++)
    {
      /* if cellborder */ /* etc. etc. */
      if (sigma[i][j] != sigma[i + 1][j])
      {
        if (g)
          g->Point(1, i + 1, j);
      }
      if (sigma[i][j] != sigma[i][j + 1])
      {
        if (g)
          g->Point(1, i, j + 1);
      }
      /* Cells that touch eachother's corners are NO neighbours */
      if (sigma[i][j] != sigma[i + 1][j + 1] || sigma[i + 1][j] != sigma[i][j + 1])
      {
        if (g)
          g->Point(1, i + 1, j + 1);
      }
    }
  }
}

int **CellularPotts::SearchNeighboursMatrix()
{
  int i, j;
  int **neighbours = new int *[cell->size() + 1];
  for (int i = 0; i < (int)cell->size() + 1; i++)
  {
    neighbours[i] = new int[cell->size() + 1];
  }
  for (i = 0; i < ((int)cell->size() + 1); i++)
  {
    for (j = 0; j < ((int)cell->size() + 1); j++)
    {
      neighbours[i][j] = 0;
    }
  }
  for (i = 0; i < sizex - 1; i++)
  {
    for (j = 0; j < sizey - 1; j++)
    {
      int iplus = i + 1;
      int jplus = j + 1;
      if (par.periodic_boundaries)
      {
        if (iplus <= 0)
        {
          iplus = sizex - 2 + iplus;
        }
        if (jplus <= 0)
        {
          jplus = sizey - 2 + jplus;
        }
        if (iplus >= sizex - 1)
        {
          iplus = iplus - sizex + 2;
        }
        if (jplus >= sizey - 1)
        {
          jplus = jplus - sizey + 2;
        }

        /* if cellborder */ /* etc. etc. */
        if (sigma[i][j] != sigma[i + 1][j])
        {
          neighbours[sigma[i][j]][sigma[iplus][j]] += 1;
          neighbours[sigma[iplus][j]][sigma[i][j]] += 1;
        }
        if (sigma[i][j] != sigma[i][j + 1])
        {
          neighbours[sigma[i][j]][sigma[i][jplus]] += 1;
          neighbours[sigma[i][jplus]][sigma[i][j]] += 1;
        }
        // if extended_neighbour_border is true, also count cells touching by a corner.
        if (par.extended_neighbour_border)
        {
          if (sigma[i][j] != sigma[i + 1][j + 1])
          {
            neighbours[sigma[i][j]][sigma[iplus][jplus]] += 1;
          }
          if (sigma[i + 1][j] != sigma[i][j + 1])
          {
            neighbours[sigma[iplus][j]][sigma[i][jplus]] += 1;
          }
        }
      }
    }
  }
  return neighbours;
}

int CellularPotts::GetNewPerimeterIfXYWereAdded(int sxyp, int x, int y)
{

  /*int n_nb;

   if (par.neighbours>=1 && par.neighbours<=4)
     n_nb=nbh_level[par.neighbours];
  */
  int perim = (*cell)[sxyp].Perimeter();

  /* the cell with sigma sxyp wants to extend by adding lattice site (x, y).
 This means that the sxyp neighbours of (x,y) will not be borders anymore,so they can be
 subtracted from the perimeter of sxyp.
*/
  for (int i = 1; i <= n_nb; i++)
  {

    int xp2, yp2;

    xp2 = x + nx[i];
    yp2 = y + ny[i];

    if (par.periodic_boundaries)
    {

      if (xp2 <= 0)
        xp2 = sizex - 2 + xp2;
      if (yp2 <= 0)
        yp2 = sizey - 2 + yp2;
      if (xp2 >= sizex - 1)
        xp2 = xp2 - sizex + 2;
      if (yp2 >= sizey - 1)
        yp2 = yp2 - sizey + 2;
    }
    if (sigma[xp2][yp2] == sxyp)
    {
      perim--;
    }
    else
    {
      perim++;
    }
  }
  return perim;
}

int CellularPotts::GetActLevel(int x, int y)
{
  if (sigma[x][y] > 0)
    return (actPixels[{x, y}]);
  else
    return (0);
}
// matrix array implementation
int CellularPotts::GetMatrixLevel(int x, int y)
{
  if (matrix[x][y] > 0)
  {
    return (matrix[x][y]);
  }
  else
  {
    return (0);
  }
}

int CellularPotts::GetNewPerimeterIfXYWereRemoved(int sxy, int x, int y)
{
  /*int n_nb;
   if (par.neighbours>=1 && par.neighbours<=4)
    int n_nb=nbh_level[par.neighbours];
  */
  int perim = (*cell)[sxy].Perimeter();
  /* the cell with sigma sxy loses xy
   */
  for (int i = 1; i <= n_nb; i++)
  {

    int xp2, yp2;
    xp2 = x + nx[i];
    yp2 = y + ny[i];
    if (par.periodic_boundaries)
    {

      if (xp2 <= 0)
        xp2 = sizex - 2 + xp2;
      if (yp2 <= 0)
        yp2 = sizey - 2 + yp2;
      if (xp2 >= sizex - 1)
        xp2 = xp2 - sizex + 2;
      if (yp2 >= sizey - 1)
        yp2 = yp2 - sizey + 2;
    }
    if (sigma[xp2][yp2] == sxy)
    {
      perim++;
    }
    else
    {
      perim--;
    }
  }
  return perim;
}
void CellularPotts::ReadZygotePicture(void)
{
  int pix, cells, i, j, c, p, checkx, checky;
  char **pixelmap;
  char pixel[3];

  sscanf(ZYGXPM(ZYGOTE)[0], "%d %d %d %d", &checkx, &checky, &cells, &pix);

  if ((checkx > sizex) || (checky > sizey))
  {
    std::cerr << "ReadZygote: The included xpm picture is smaller than the grid!\n";
    std::cerr << "\n Please adjust either the grid size or the picture size.\n";
    std::cerr << sizex << "," << sizey << "," << checkx << "," << checky << "\n";
    exit(1);
  }
  pixelmap = (char **)malloc(cells * sizeof(char *));
  if (pixelmap == NULL)
    MemoryWarning();

  pixelmap[0] = (char *)malloc(cells * 3 * sizeof(char));
  if (pixelmap[0] == NULL)
    MemoryWarning();

  for (i = 1; i < cells; i++)
    pixelmap[i] = pixelmap[i - 1] + 3;

  for (i = 0; i < cells; i++)
  {
    for (j = 0; j < pix; j++)
      pixelmap[i][j] = ZYGXPM(ZYGOTE)[i + 1][j];
    pixelmap[i][pix] = '\0';
  }

  for (i = 0; i < sizex * sizey; i++)
  {
    sigma[0][i] = 0;
    tau[0][i] = 0;
  }
  fprintf(stderr, "[%d %d]\n", checkx, checky);

  int offs_x, offs_y;
  offs_x = (sizex - checkx) / 2;
  offs_y = (sizey - checky) / 2;

  for (i = 0; i < checkx; i++)
    for (j = 0; j < checky; j++)
    {
      for (p = 0; p < pix; p++)
        pixel[p] = ZYGXPM(ZYGOTE)[cells + 1 + j][i * pix + p];

      pixel[pix] = '\0';

      for (c = 0; c < cells; c++)
      {
        if (!(strcmp(pixelmap[c], pixel)))
        {
          if ((sigma[offs_x + i][offs_y + j] = c))
          {

            // if c is _NOT_ medium (then c=0)
            // assign pixel values from "sigmamax"
            sigma[offs_x + i][offs_y + j] += (Cell::MaxSigma() - 1);
            tau[offs_x + i][offs_y + j] = (*cell)[sigma[offs_x + i][offs_y + i]].getTau();
          }
        }
      }
    }
  free(pixelmap[0]);
  free(pixelmap);
}

void CellularPotts::ConstructInitCells(Dish &beast)
{

  // Get the maximum cell ID (mostly equal to the cell number)
  int loop = sizex * sizey;
  int cells = 0;
  for (int i = 0; i < loop; i++)
  {
    if (cells < sigma[0][i])
      cells = sigma[0][i];
  }

  cerr << "[ cells = " << cells << "]\n";

  // construct enough cells for the zygote.  "cells", contains the
  // number of colours (excluding background).
  {
    for (int i = 0; i < cells; i++)
    {
      cell->push_back(Cell(beast));
    }
  }

  // Set the area and target area of the cell
  // makes use of the pointer to the Cell pointer of Dish
  // which is a member of CellularPotts
  MeasureCellSizes();

  // set zygote_area to mean cell area.
  int mean_area = 0;
  for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
  {
    mean_area += c->Area();
  }
  if (cells != 0)
    mean_area /= cells;

  zygote_area = mean_area;

  cout << "mean_area = " << mean_area << "\n";
  // set all cell areas to the mean area
  {
    for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
    {
      if (par.target_area >= 0)
      { 
        c->SetTargetArea(par.target_area);
      }
      else
      {
        c->SetTargetArea(mean_area);
      }
    }
  }
}

void CellularPotts::ConstructInitCellGrid(Dish &beast)
{  

  // Get the maximum cell ID (mostly equal to the cell number)
  int loop = sizex * sizey;
  int cells = 0;
  int cells1 = 0;
  int cells2 = 0;
  for (int i = 0; i < loop; i++)
  {
    if (cells < sigma[0][i])
      cells = sigma[0][i];
  }


  cerr << "[ cells = " << cells << "]\n";

  // construct enough cells for the zygote.  "cells", contains the
  // number of colours (excluding background).
  {
    for (int i = 0; i < cells; i++)
    {
      cell->push_back(Cell(beast));
    }
  }

  // Set the area and target area of the cell
  // makes use of the pointer to the Cell pointer of Dish
  // which is a member of CellularPotts
  MeasureCellSizes();
  if (par.second_layer)
    SetTypesWithDoubleMask();
  else
    SetTypesWithMask();
  // set zygote_area to mean cell area.

  int mean_area_1 = 0;
  int mean_area_2 = 0;
  for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
  {
    if (c->getTau() == 1){
      mean_area_1 += c->Area();
      cells1++;
    }
    else if (c->getTau() == 2){
      mean_area_2 += c->Area();
      cells2++;
      
    }
  }
  if (cells1 != 0)
    mean_area_1 /= cells1;
  if (cells2 != 0) 
    mean_area_2 /= cells2;
  // set all cell areas to the mean area
  {
    for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
    {
      if (c->getTau() == 1)
        c->SetTargetArea(mean_area_1);
      if (c->getTau() == 2)
        c->SetTargetArea(mean_area_2);    
    }
  }
  ResetTargetLengths();
}


void CellularPotts::GrowCellGrid(Dish &beast){
  DetectSidesIsthmus();
  int delta_x1 = par.celltype1_length;
  int delta_y1 = par.celltype1_width;
  int delta_x2 = par.celltype2_length;
  int delta_y2 = par.celltype2_width;

  cout << "delta_x1 = " << delta_x1 << endl;
  cout << "delta_y1 = " << delta_y1 << endl;
  cout << "delta_x2 = " << delta_x2 << endl;
  cout << "delta_y2 = " << delta_y2 << endl;

  int cell_number = 1;
  bool added_cell = false;
  for (int grid_x = 0; grid_x < right_side_isthmus; grid_x += delta_x2)
    for (int grid_y = 0; grid_y < sizey; grid_y += delta_y2){
      for (int x = 0; x < delta_x2; x++)
        for (int y = 0; y < delta_y2; y++){
          if (grid_x+x < right_side_isthmus && grid_y+y < sizey)
            if (mask[grid_x+x][grid_y+y]){
              added_cell = true;
              sigma[grid_x+x][grid_y+y] = cell_number;
            }        
        }
      if(added_cell)
        cell_number++;
      added_cell = false;
    }
  for (int grid_x = right_side_isthmus; grid_x < sizex; grid_x += delta_x1)
    for (int grid_y = 0; grid_y < sizey; grid_y += delta_y1){
      for (int x = 0; x < delta_x1; x++)
        for (int y = 0; y < delta_y1; y++){
          if (grid_x+x < sizex && grid_y+y < sizey)
            if (mask[grid_x+x][grid_y+y]){
              added_cell = true;
              sigma[grid_x+x][grid_y+y] = cell_number;
            }        
        }
      if(added_cell)
        cell_number++;
      added_cell = false;
    }
}



void CellularPotts::GrowCellGridOnLayers(Dish &beast){
  int delta_x1 = par.celltype1_length;
  int delta_y1 = par.celltype1_width;
  int delta_x2 = par.celltype2_length;
  int delta_y2 = par.celltype2_width;

  cout << "delta_x1 = " << delta_x1 << endl;
  cout << "delta_y1 = " << delta_y1 << endl;
  cout << "delta_x2 = " << delta_x2 << endl;
  cout << "delta_y2 = " << delta_y2 << endl;

  int cell_number = 1;
  bool added_cell = false;
  for (int grid_x = 0; grid_x < sizex; grid_x += delta_x1)
    for (int grid_y = 0; grid_y < sizey; grid_y += delta_y1){
      for (int x = 0; x < delta_x1; x++)
        for (int y = 0; y < delta_y1; y++){
          if(grid_x+x < sizex && grid_y+y < sizey)
            if (mask[grid_x+x][grid_y+y] == 1){
                added_cell = true;
                sigma[grid_x+x][grid_y+y] = cell_number;
              }        
        }
      if(added_cell)
        cell_number++;
      added_cell = false;
    }
  for (int grid_x = 0; grid_x < sizex; grid_x += delta_x2)
    for (int grid_y = 0; grid_y < sizey; grid_y += delta_y2){
      for (int x = 0; x < delta_x2; x++)
        for (int y = 0; y < delta_y2; y++){
          if(grid_x+x < sizex && grid_y+y < sizey)
            if (mask[grid_x+x][grid_y+y] == 2){
                added_cell = true;
                sigma[grid_x+x][grid_y+y] = cell_number;
              }        
        }
      if(added_cell)
        cell_number++;
      added_cell = false;
    }  
}

void CellularPotts::RefreshLinks(void)
{
	vector<Cell>::iterator c = cell->begin();
	++c;
	vector<Cell>::iterator z;
	int p;
	int number_of_links;
	int number_of_cells = 0;
	double centX_1, centY_1, centX_2, centY_2, deltaX, deltaY;
	double polX, polY;
	c = cell->begin();
	++c;
	for (; c != cell->end(); c++)
		number_of_cells++;
	c = cell->begin();
	++c;
	for (; c != cell->end(); c++)
	{
		z = cell->begin();
		++z;
		number_of_links = 0;

		centX_1 = c->getCenterX();
		centY_1 = c->getCenterY();
		polX = cos(c->polarization);
		polY = sin(c->polarization);

		for (; z != cell->end(); z++)
		{
			switch (par.pulling_method)
			{
			case 1:
				if (z->tau == 1)
				{
					centX_2 = z->getCenterX();
					centY_2 = z->getCenterY();
					deltaX = centX_2 - centX_1;
					deltaY = centY_2 - centY_1;
					c->Links[z->Sigma()] = false;
					if ((sqrt(deltaX * deltaX + deltaY * deltaY) < par.r_max) && (c->tau == 1) &&
						((acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180) ||
						 (PI - acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180)))
					{
						number_of_links++;
						c->Links[z->Sigma()] = true;
					}
				}
				if (z->tau == 2)
				{
					centX_2 = z->getCenterX();
					centY_2 = z->getCenterY();
					deltaX = centX_2 - centX_1;
					deltaY = centY_2 - centY_1;
					c->Links[z->Sigma()] = false;
					if ((sqrt(deltaX * deltaX + deltaY * deltaY) < par.r_max) && (c->tau == 2) &&
						((acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180) ||
						 (PI - acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180)))
					{
						number_of_links++;
						c->Links[z->Sigma()] = true;
					}
				}
				break;
			case 2:
				centX_2 = z->getCenterX();
				centY_2 = z->getCenterY();
				deltaX = centX_2 - centX_1;
				deltaY = centY_2 - centY_1;
				c->Links[z->Sigma()] = false;
				if ((sqrt(deltaX * deltaX + deltaY * deltaY) < par.r_max) &&
					((acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180) ||
					 (PI - acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180)))
				{
					if (!c->sigma || !z->sigma)
						cout << "c->Sigma() = " << c->sigma << ", z-Sigma() = " << z->sigma << " and  sqrt(deltaX * deltaX + deltaY * deltaY) = " << sqrt(deltaX * deltaX + deltaY * deltaY) << endl;
					number_of_links++;
					c->Links[z->Sigma()] = true;
				}
				break;
			case 3:
				centX_2 = z->getCenterX();
				centY_2 = z->getCenterY();
				deltaX = centX_2 - centX_1;
				deltaY = centY_2 - centY_1;
				c->Links[z->Sigma()] = false;
				if ((sqrt(deltaX * deltaX + deltaY * deltaY) < par.r_max) && (c->tau == 2) &&
					((acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180) ||
					 (PI - acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180)))
				{
					number_of_links++;
					c->Links[z->Sigma()] = true;
				}
				break;
			case 4:
				if (z->tau == 2)
				{
					centX_2 = z->getCenterX();
					centY_2 = z->getCenterY();
					deltaX = centX_2 - centX_1;
					deltaY = centY_2 - centY_1;
					c->Links[z->Sigma()] = false;
					if ((sqrt(deltaX * deltaX + deltaY * deltaY) < par.r_max) && (c->tau == 2) &&
						((acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180) ||
						 (PI - acos((deltaX * polX + deltaY * polY) / (sqrt(deltaX * deltaX + deltaY * deltaY) * sqrt(polX * polX + polY * polY))) < par.theta_max * PI / 180)))
					{
						number_of_links++;
						c->Links[z->Sigma()] = true;
					}
				}
			}
		}
		while (number_of_links > par.max_links)
		{
			p = rand() % number_of_cells + 1;
			if (c->Links[p])
			{
				number_of_links--;
				c->Links[p] = false;
			}
			// randomly remove links until you only have n_max left
		}
		// Adjust polarizations according to your neighbours.
		if (number_of_links)
		{
			bool belmonte = false; //Use the repolarization method from Belmonte et al.?


			
			z = cell->begin();
			++z;
			double Avgpolarization = 0;
			double AvgX = 0;
			double AvgY = 0;
			double distance = 0;
			if (belmonte){

				for (; z != cell->end(); z++){	
					if (c->Links[z->Sigma()]){
						distance = sqrt(pow(c->getCenterX() - z->getCenterX(),2) + pow(c->getCenterY() - z->getCenterY(),2));
						if (c->tau == 1 && z->tau == 1)
						{
							AvgX += (z->getCenterX()-c->getCenterX()) / distance;
							AvgY += (z->getCenterY()-c->getCenterY()) / distance;
						}
						else if (c->tau == 1 && z->tau == 2)
						{
							AvgX += (z->getCenterX()-c->getCenterX()) / distance;
							AvgY += (z->getCenterY()-c->getCenterY()) / distance;
						}
						else if (c->tau == 2 && z->tau == 1)
						{
							AvgX += (z->getCenterX()-c->getCenterX()) / distance;
							AvgY += (z->getCenterY()-c->getCenterY()) / distance;
						}
						else if (c->tau == 2 && z->tau == 2)
						{
							AvgX += (z->getCenterX()-c->getCenterX()) / distance;
							AvgY += (z->getCenterY()-c->getCenterY()) / distance;
						}
					}
				}
			}

			else{
				for (; z != cell->end(); z++)
				{
					if (c->Links[z->Sigma()] && c->tau == 1 && z->tau == 1)
					{
						AvgX += cos(z->polarization);
						AvgY += sin(z->polarization);
					}
					else if (c->Links[z->Sigma()] && c->tau == 1 && z->tau == 2)
					{
						AvgX += cos(z->polarization);
						AvgY += sin(z->polarization);
					}
					else if (c->Links[z->Sigma()] && c->tau == 2 && z->tau == 1)
					{
						AvgX += cos(z->polarization);
						AvgY += sin(z->polarization);
					}
					else if (c->Links[z->Sigma()] && c->tau == 2 && z->tau == 2)
					{
						AvgX += cos(z->polarization);
						AvgY += sin(z->polarization);
					}
				}
			}

				if (AvgX > 0 && AvgY > 0)
					Avgpolarization = atan(AvgY / AvgX);
				else if (AvgX < 0 && AvgY > 0)
					Avgpolarization = atan(-AvgX / AvgY) + PI / 2;
				else if (AvgX < 0 && AvgY < 0)
					Avgpolarization = atan(AvgY / AvgX) + PI;
				else
					Avgpolarization = atan(-AvgX / AvgY) + 3 * PI / 2;

				double NewX, NewY;
				NewX = par.memory * cos(c->polarization) + (1 - par.memory) * cos(Avgpolarization);
				NewY = par.memory * sin(c->polarization) + (1 - par.memory) * sin(Avgpolarization);

				if (NewX > 0 && NewY > 0)
					c->polarization = atan(NewY / NewX);
				else if (NewX < 0 && NewY > 0)
					c->polarization = atan(-NewX / NewY) + PI / 2;
				else if (NewX < 0 && NewY < 0)
					c->polarization = atan(NewY / NewX) + PI;
				else
					c->polarization = atan(-NewX / NewY) + 3 * PI / 2;

			
		}
	}
}

void CellularPotts::MeasureCellSizes(void)
{
  // Clean areas of all cells, including medium
  for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
  {
    c->SetTargetArea(0);
    c->area = 0;
  }

  // calculate the area of the cells
  for (int x = 1; x < sizex - 1; x++)
  {
    for (int y = 1; y < sizey - 1; y++)
    {
      if (sigma[x][y] > 0)
      {
        (*cell)[sigma[x][y]].IncrementTargetArea();
        (*cell)[sigma[x][y]].IncrementArea();
        (*cell)[sigma[x][y]].AddSiteToMoments(x, y);
      }
    }
  }

  // set the actual area to the target area
  for (vector<Cell>::iterator c = cell->begin(); c != cell->end(); c++)
  {
    c->SetAreaToTarget();
  }
}

void CellularPotts::MeasureCellSize(Cell &c)
{
  c.CleanMoments();
  // calculate the area of the cell
  for (int x = 1; x < sizex - 1; x++)
  {
    for (int y = 1; y < sizey - 1; y++)
    {
      if (sigma[x][y] == c.sigma)
      {
        (*cell)[sigma[x][y]].IncrementTargetArea();
        (*cell)[sigma[x][y]].IncrementArea();
        (*cell)[sigma[x][y]].AddSiteToMoments(x, y);
      }
    }
  }
}

void CellularPotts::MeasureCellPerimeters()
{
  for (int x = 1; x < sizex - 1; x++)
  {
    for (int y = 1; y < sizey - 1; y++)
    {
      if (sigma[x][y] > 0)
      {
        for (int i = 1; i <= n_nb; i++)
        {
          int xp2, yp2;
          xp2 = x + nx[i];
          yp2 = y + ny[i];
          if (par.periodic_boundaries)
          {
            if (xp2 <= 0)
              xp2 = sizex - 2 + xp2;
            if (yp2 <= 0)
              yp2 = sizey - 2 + yp2;
            if (xp2 >= sizex - 1)
              xp2 = xp2 - sizex + 2;
            if (yp2 >= sizey - 1)
              yp2 = yp2 - sizey + 2;
          }
          // did we find a border?
          if (sigma[xp2][yp2] != sigma[x][y])
          {
            // add to the perimeter of the cell
            (*cell)[sigma[x][y]].IncrementTargetPerimeter();
            (*cell)[sigma[x][y]].IncrementPerimeter();
          }
        }
      }
    }
  }
}

Dir *CellularPotts::FindCellDirections(void) const
{
  double *sumx = 0, *sumy = 0;
  double *sumxx = 0, *sumxy = 0, *sumyy = 0;
  double *n = 0;
  double xmean = 0, ymean = 0, sxx = 0, sxy = 0, syy = 0;
  double D, lb1 = 0, lb2 = 0;

  Dir *celldir;

  /* Allocation of sufficient memory space */
  if ((sumx = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();
  else if ((sumy = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();
  else if ((sumxx = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();
  else if ((sumxy = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();
  else if ((sumyy = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();
  else if ((n = (double *)malloc(cell->size() * sizeof(double))) == NULL)
    MemoryWarning();

  if (!(celldir = new Dir[cell->size()]))
    MemoryWarning();

  /* Initialization of the variables */
  for (int i = 0; i < (int)cell->size(); i++)
  {
    sumx[i] = 0.;
    sumy[i] = 0.;
    sumxx[i] = 0.;
    sumxy[i] = 0.;
    sumyy[i] = 0.;
    n[i] = 0L;
  }

  /* Find sumx, sumy, sumxx and sumxy for all cells */
  for (int x = 0; x < sizex; x++)
    for (int y = 0; y < sizey; y++)
      if (sigma[x][y] > 0)
      {
        sumx[0] += (double)x;
        sumy[0] += (double)y;
        sumxx[0] += (double)x * x;
        sumxy[0] += (double)x * y;
        sumyy[0] += (double)y * y;

        n[0]++;

        sumx[sigma[x][y]] += (double)x;
        sumy[sigma[x][y]] += (double)y;

        sumxx[sigma[x][y]] += (double)x * x;
        sumxy[sigma[x][y]] += (double)x * y;
        sumyy[sigma[x][y]] += (double)y * y;

        n[sigma[x][y]]++;
      }

  /* Compute the principal axes for all cells */
  for (int i = 0; i < (int)cell->size(); i++)
  {
    if (n[i] > 10)
    {
      xmean = ((double)sumx[i]) / ((double)n[i]);
      ymean = ((double)sumy[i]) / ((double)n[i]);

      sxx = (double)(sumxx[i]) - ((double)(sumx[i] * sumx[i])) / (double)n[i];
      sxx = sxx / (double)(n[i] - 1);

      sxy = (double)(sumxy[i]) - ((double)(sumx[i] * sumy[i])) / (double)n[i];
      sxy = sxy / (double)(n[i] - 1);

      syy = (double)(sumyy[i]) - ((double)(sumy[i] * sumy[i])) / (double)n[i];
      syy = syy / (double)(n[i] - 1);

      D = sqrt((sxx + syy) * (sxx + syy) - 4. * (sxx * syy - sxy * sxy));
      lb1 = (sxx + syy + D) / 2.;
      lb2 = (sxx + syy - D) / 2.;
      celldir[i].lb1 = lb1;
      celldir[i].lb2 = lb2;
    }
    if (sxy == 0.0)
      celldir[i].bb1 = 1.;
    else
      celldir[i].bb1 = sxy / (lb1 - syy);

    if (fabs(celldir[i].bb1) < .00001)
    {
      if (celldir[i].bb1 > 0.)
        celldir[i].bb1 = .00001;
      else
        celldir[i].bb1 = -.00001;
    }
    celldir[i].aa1 = ymean - xmean * celldir[i].bb1;
    celldir[i].bb2 = (-1.) / celldir[i].bb1;
    celldir[i].aa2 = ymean - celldir[i].bb2 * xmean;
  }

  /* bevrijd gealloceerd geheugen */
  free(sumx);
  free(sumy);
  free(sumxx);
  free(sumxy);
  free(sumyy);
  free(n);

  return celldir;
}

void CellularPotts::ShowDirections(Graphics &g, const Dir *celldir) const
{
  int i;
  if (cell->size() > 1)
    for (i = 1; i < (int)cell->size(); i++)
      g.Line(0, (int)(2 * celldir[i].aa1), sizex * 2, (int)((celldir[i].aa1 + celldir[i].bb1 * sizey) * 2), 2);
}

void CellularPotts::DivideCells(vector<bool> which_cells)
{

  // for the cell directions
  Dir *celldir = 0;

  /* Allocate space for divisionflags */
  int *divflags = (int *)malloc((cell->size() * 2 + 5) * sizeof(int));

  /* Clear divisionflags */
  for (int i = 0; i < (int)(cell->size() * 2 + 5); i++)
    divflags[i] = 0;

  if (!(which_cells.size() == 0 || which_cells.size() >= cell->size()))
  {
    throw "In CellularPotts::DivideCells, Too few elements in vector<int> which_cells.";
  }

  /* division */
  for (int i = 0; i < sizex; i++)
  {
    for (int j = 0; j < sizey; j++)
      if (sigma[i][j] > 0)
      { // i.e. not medium and not border state (-1)
        // Pointer to mother. Warning: Renew pointer after a new
        // cell is added (push_back). Then, the array *cell is relocated and
        // the pointer will be lost...

        Cell *motherp = &((*cell)[sigma[i][j]]);
        Cell *daughterp;

        /* Divide if NOT medium and if DIV bit set or divide_always is set */
        // if which_cells is given, divide only if the cell
        // is marked in which_cells.
        if (!which_cells.size() || which_cells[motherp->sigma])
        {
          if (!(divflags[motherp->Sigma()]))
          {
            // add daughter cell, copying states of mother
            daughterp = new Cell(*(motherp->owner));
            daughterp->CellBirth(*motherp);
            cell->push_back(*daughterp);

            // renew pointer to mother
            motherp = &((*cell)[sigma[i][j]]);

            divflags[motherp->Sigma()] = daughterp->Sigma();
            delete daughterp;

            // array may be relocated after "push_back"

            // renew daughter pointers
            daughterp = &(cell->back());

            /* administration on the onset of mitosis */

            /* Ancestry is taken care of in copy constructor of Cell
               see cell.hh: Cell(const Cell &src, bool newcellP=false) : Cytoplasm(src) {} */

            /* inherit  polarity of mother */
            // All that needs to be copied is copied in the copy constructor
            // of Cell and in the default copy constr. of its base class Cytoplasm
            // note: also the celltype is inherited
          }
          else
          {
            daughterp = &((*cell)[divflags[motherp->Sigma()]]);
          }

          /* Now the actual division takes place */

          /* If celldirections where not yet computed: do it now */
          if (!celldir)
            celldir = FindCellDirections();

          /* if site is below the minor axis of the cell: sigma of new cell */
          if (j > ((int)(celldir[motherp->sigma].aa2 +
                         celldir[motherp->sigma].bb2 * (double)i)))
          {
            motherp->DecrementArea();
            motherp->DecrementTargetArea();
            motherp->RemoveSiteFromMoments(i, j);
            sigma[i][j] = daughterp->Sigma();
            tau[i][j] = (*cell)[sigma[i][j]].getTau();
            daughterp->AddSiteToMoments(i, j);
            daughterp->IncrementArea();
            daughterp->IncrementTargetArea();
          }
        }
      }
  }
  if (celldir)
    delete[](celldir);

  if (divflags)
    free(divflags);
}

/**! Fill the plane with initial cells
 \return actual amount of cells (some are not draw due to overlap) */
int CellularPotts::ThrowInCells(int n, int cellsize)
{

  //  int gapx=(sizex-nx*cellsize)/(nx+1);
  // int gapy=(sizey-ny*cellsize)/(ny+1);

  int cellnum = 1;

  for (int i = 0; i < n; i++)
  {
    // draw a circle at x0, y0
    int x0 = RandomNumber(sizex);
    int y0 = RandomNumber(sizey);

    bool overlap = false;

    // check overlap
    for (int x = 0; x < cellsize; x++)
      for (int y = 0; y < cellsize; y++)
        if ((((x - cellsize / 2) * (x - cellsize / 2) + (y - cellsize / 2) * (y - cellsize / 2)) <
             ((cellsize / 2) * (cellsize / 2))) &&
            (x0 + x < sizex && y0 + y < sizey))
          if (sigma[x0 + x][y0 + y])
          {
            overlap = true;
            break;
          }

    if (!overlap)
    {
      for (int x = 0; x < cellsize; x++)
        for (int y = 0; y < cellsize; y++)
          if ((((x - cellsize / 2) * (x - cellsize / 2) + (y - cellsize / 2) * (y - cellsize / 2)) <
               ((cellsize / 2) * (cellsize / 2))) &&
              (x0 + x < sizex && y0 + y < sizey))
          {
            sigma[x0 + x][y0 + y] = cellnum;
            tau[x0 + x][y0 + y] = (*cell)[cellnum].getTau();
          }
      cellnum++;
    }
  }
  cerr << "[ cellnum = " << cellnum << "]";

  if (par.micropatternmask != string("None"))
  {
    for (int x = 0; x < sizex; x++)
    {
      for (int y = 0; y < sizey; y++)
      {
        if (!mask[x][y])
        {
          sigma[x][y] = -1;
          tau[x][y] = -1;
        }
        else
        {
          sigma[x][y] = 0;
          tau[x][y] = 0;
        }
      }
    }
  }

  else
  {
    // repair borders
    // fill borders with special border state
    for (int x = 0; x < sizex - 1; x++)
    {
      sigma[x][0] = -1;
      sigma[x][sizey - 1] = -1;
      tau[x][0] = -1;
      tau[x][sizey - 1] = -1;
    }
    for (int y = 0; y < sizey - 1; y++)
    {
      sigma[0][y] = -1;
      sigma[sizex - 1][y] = -1;
      tau[0][y] = -1;
      tau[sizex - 1][y] = -1;
    }
    for (int x = 1; x < sizex - 2; x++)
    {
      sigma[x][1] = 0;
      sigma[x][sizey - 2] = 0;
      tau[x][1] = 0;
      tau[x][sizey - 2] = 0;
    }
    for (int y = 1; y < sizey - 2; y++)
    {
      sigma[1][y] = 0;
      sigma[sizex - 2][y] = 0;
      tau[1][y] = 0;
      tau[sizex - 2][y] = 0;
    }
  }
  return cellnum;
}

int CellularPotts::GrowInCellsInMicropattern(int n_cells, int cell_size)
{

  // make initial cells using Eden Growth

  int **new_sigma = (int **)malloc(sizex * sizeof(int *));
  if (new_sigma == NULL)
    MemoryWarning();

  new_sigma[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (new_sigma[0] == NULL)
    MemoryWarning();

  for (int i = 1; i < sizex; i++)
    new_sigma[i] = new_sigma[i - 1] + sizey;

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      new_sigma[0][i] = 0;
  }

  // scatter initial points, or place a cell in the middle
  // if only one cell is desired
  int cellnum = cell->size() - 1;

  bool placed;
  int xposition;
  int yposition;
  for (int i = 0; i < n_cells; i++)
  {
    placed = false;
    while (placed == false)
    {
      xposition = 1 + RandomNumber(sizex - 2);
      yposition = 1 + RandomNumber(sizey - 2);
      if (mask[xposition][yposition])
      {
        sigma[xposition][yposition] = ++cellnum;
        // cout << "tau[xposition][yposition] = " << tau[xposition][yposition] << endl;
        // cout << "(*cell)[sigma[xposition][yposition]].getTau()= " << (*cell)[sigma[xposition][yposition]].getTau() << endl;
        // tau[xposition][yposition]=(*cell)[sigma[xposition][yposition]].getTau();
        placed = true;
      }
    }
  }

  // Do Eden growth for a number of time steps
  {
    for (int i = 0; i < cell_size; i++)
    {
      for (int x = 1; x < sizex - 1; x++)
        for (int y = 1; y < sizey - 1; y++)
        {

          if (sigma[x][y] == 0)
          {
            // take a random neighbour
            int xyp = (int)(8 * RANDOM() + 1);
            int xp = nx[xyp] + x;
            int yp = ny[xyp] + y;
            int kp;
            //  NB removing this border test yields interesting effects :-)
            // You get a ragged border, which you may like!
            if ((kp = sigma[xp][yp]) != -1)
              if (kp > (cellnum - n_cells))
                new_sigma[x][y] = kp;
              else
                new_sigma[x][y] = 0;
            else
              new_sigma[x][y] = 0;
          }
          else
          {
            new_sigma[x][y] = sigma[x][y];
          }
        }

      // copy sigma to new_sigma, but do not touch the border!
      {
        for (int x = 1; x < sizex - 1; x++)
        {
          for (int y = 1; y < sizey - 1; y++)
          {
            sigma[x][y] = new_sigma[x][y];
          }
        }
      }
    }
  }

  free(new_sigma[0]);
  free(new_sigma);
  return cellnum;
}

void CellularPotts::DetectSidesIsthmus(){
  bool interrupted = false;
  for (int y = sizey/2; y < sizey; y++){
    interrupted = false;
    for (int x = sizex -1; x > 0; x--){
      if(mask[x][y] && !mask[x-1][y] && !interrupted){
        interrupted = true;
        right_side_isthmus = x-1;
      }
      else if(!mask[x][y] && mask[x-1][y] && interrupted){
        left_side_isthmus = x;
        return;
      }    
    }      
  }

}

int CellularPotts::GrowInCells(int n_cells, int cell_size, double subfield, int posx, int posy)
{
  int sx = (int)((sizex - 2) / subfield);
  int sy = (int)((sizey - 2) / subfield);

  int offset_x = (sizex - 2 - sx) / 2;
  int offset_y = (sizey - 2 - sy) / 2;

  if (n_cells == 1)
  {
    if (posx < 0)
      posx = sizex / 2;
    if (posy < 0)
      posy = sizey / 2;
    return GrowInCells(1, cell_size, posx, posy, 0, 0);
  }
  else
  {
    return GrowInCells(n_cells, cell_size, sx, sy, offset_x, offset_y);
  }
}

void CellularPotts::RandomSpins(double prob)
{
  for (int x = 1; x <= sizex - 2; x++)
  {
    for (int y = 1; y < sizey - 2; y++)
    {
      sigma[x][y] = (RANDOM() < prob) ? 0 : 1;
      tau[x][y] = (*cell)[sigma[x][y]].getTau();
    }
  }
  cerr << "RandomSpins done" << endl;
}

int CellularPotts::GrowInCells(int n_cells, int cell_size, int sx, int sy, int offset_x, int offset_y)
{
  int randomx;
  int randomy;
  // make initial cells using Eden Growth

  int **new_sigma = (int **)malloc(sizex * sizeof(int *));
  if (new_sigma == NULL)
    MemoryWarning();

  new_sigma[0] = (int *)malloc(sizex * sizey * sizeof(int));
  if (new_sigma[0] == NULL)
    MemoryWarning();

  for (int i = 1; i < sizex; i++)
    new_sigma[i] = new_sigma[i - 1] + sizey;

  /* Clear CA plane */
  {
    for (int i = 0; i < sizex * sizey; i++)
      new_sigma[0][i] = 0;
  }

  // scatter initial points, or place a cell in the middle
  // if only one cell is desired
  int cellnum = cell->size() - 1;

  if (n_cells > 1)
  {

    {
      for (int i = 0; i < n_cells; i++)
      {
        randomx = RandomNumber(sx);
        randomy = RandomNumber(sy);
        sigma[randomx + offset_x][randomy + offset_y] = ++cellnum;
        tau[randomx + offset_x][randomy + offset_y] = (*cell)[sigma[randomx + offset_x][randomy + offset_y]].getTau();
      }
    }
  }
  else
  {
    sigma[sx][sy] = ++cellnum;
    tau[sx][sy] = (*cell)[sigma[sx][sy]].getTau();
  }

  // Do Eden growth for a number of time steps
  {
    for (int i = 0; i < cell_size; i++)
    {
      for (int x = 1; x < sizex - 1; x++)
        for (int y = 1; y < sizey - 1; y++)
        {

          if (sigma[x][y] == 0)
          {
            // take a random neighbour
            int xyp = (int)(8 * RANDOM() + 1);
            int xp = nx[xyp] + x;
            int yp = ny[xyp] + y;
            int kp;
            //  NB removing this border test yields interesting effects :-)
            // You get a ragged border, which you may like!
            if ((kp = sigma[xp][yp]) != -1)
              if (kp > (cellnum - n_cells))
                new_sigma[x][y] = kp;
              else
                new_sigma[x][y] = 0;
            else
              new_sigma[x][y] = 0;
          }
          else
          {
            new_sigma[x][y] = sigma[x][y];
          }
        }

      // copy sigma to new_sigma, but do not touch the border!
      {
        for (int x = 1; x < sizex - 1; x++)
        {
          for (int y = 1; y < sizey - 1; y++)
          {
            sigma[x][y] = new_sigma[x][y];
            tau[x][y] = (*cell)[sigma[x][y]].getTau();
          }
        }
      }
    }
  }
  free(new_sigma[0]);
  free(new_sigma);

  return cellnum;
}

/** Draw a square cell in at (cx,cy) */
int CellularPotts::SquareCell(int sig, int cx, int cy, int size)
{
  int xmin, xmax;
  xmin = cx - size / 2;
  if (xmin < 1)
    xmin = 1;
  xmax = cx + size / 2;
  if (xmax > sizex - 1)
    xmax = sizex - 1;

  int ymin, ymax;
  ymin = cy - size / 2;
  if (ymin < 1)
    ymin = 1;
  ymax = cy + size / 2;
  if (ymax > sizey - 1)
    ymax = sizey - 1;

  for (int x = xmin; x <= xmax; x++)
  {
    for (int y = ymin; y <= ymax; y++)
    {
      sigma[x][y] = sig;
      tau[x][y] = (*cell)[sig].getTau();
    }
  }
  return 1;
}

bool CellularPotts::LocalConnectedness(int x, int y, int s){

  //Algorithm from Durand, M., & Guesnet, E. (2016). An efficient Cellular Potts Model algorithm that forbids cell fragmentation. Computer Physics Communications, 208, 54-63.
  //Checks if cell sigma is locally connected at lattice point (x,y)
   // Use local nx and ny in a cyclic order (starts at upper left corner)
  const int cyc_nx[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
  const int cyc_ny[8] = {0, -1, -1, -1, 0, 1, 1, 1};
  bool connected_component = false; 
  //Currently in a connected component
  int nr_connected_components = 0;
  //Total number of conncected components around x,y
  for (int i = 0; i <= 7; i++){
    int s_nb = sigma[x + cyc_nx[i]][y + cyc_ny[i]];
    if (s_nb == s && !connected_component){
      //start of a connected component
      connected_component = true;
      nr_connected_components++;
    }
    else if (s_nb != s && connected_component){
      //end of a conencted component
      connected_component = false;
    }
  }
  bool looped = false;
  if (sigma[x + cyc_nx[0]][y + cyc_ny[0]] == s && sigma[x + cyc_nx[7]][y + cyc_ny[7]] == s)
    looped = true;
  //Check if the first and last element are connected
  if ((nr_connected_components >= 2 && !looped) || nr_connected_components >= 3 && looped)
  //permit one more component when the first and last element are connected
    return false;
  else
    return true;
}

// Predicate returns true when connectivity is locally preserved
// if the value of the central site would be changed
bool CellularPotts::ConnectivityPreservedP(int x, int y)
{ 
  // Use local nx and ny in a cyclic order (starts at upper left corner)
  // first site is repeated, for easier looping
  const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1};
  const int cyc_ny[10] = {0, -1, -1, -1, 0, 1, 1, 1, 0, -1};

  int sxy = sigma[x][y]; // the central site
  if (sxy == 0)
    return true;

  int n_borders = 0; // to count the amount of sites in state sxy bordering a site !=sxy

  static int stack[8]; // stack to count number of different surrounding cells
  int stackp = -1;
  bool one_of_neighbors_medium = false;
  for (int i = 1; i <= 8; i++)
  {
    int s_nb = sigma[x + cyc_nx[i]][y + cyc_ny[i]];
    int s_next_nb = sigma[x + cyc_nx[i + 1]][y + cyc_ny[i + 1]];

    if ((s_nb == sxy || s_next_nb == sxy) && (s_nb != s_next_nb))
    {
      // check whether s_nb is adjacent to non-identical site,
      // count it
      n_borders++;
    }
    int j;
    bool on_stack_p = false;

    // we need the next heuristic to prevent stalling at
    // cell-cell borders
    // do not enforce constraint at two cell interface(no medium)
    if (s_nb)
    {
      for (j = stackp; j >= 0; j--)
      {
        if (s_nb == stack[j])
        {
          on_stack_p = true;
          break;
        }
      }
      if (!on_stack_p)
      {
        if (stackp > 6)
        {
          cerr << "Stack overflow, stackp=" << stackp << "\n";
        }
        stack[++stackp] = s_nb;
      }
    }
    else
    {
      one_of_neighbors_medium = true;
    }
  }

  // number of different neighbours is stackp+1;
  if (n_borders > 2 && ((stackp + 1) > 2 || one_of_neighbors_medium))
  {
    return false;
  }
  else
    return true;
}

// Predicate returns true when cluster connectivity is locally preserved
// if the value of the central site would be changed
bool CellularPotts::ConnectivityPreservedPCluster(int x, int y)
{

  // Use local nx and ny in a cyclic order (starts at upper left corner)
  // first site is repeated, for easier looping
  const int cyc_nx[10] = {-1, -1, 0, 1, 1, 1, 0, -1, -1, -1};
  const int cyc_ny[10] = {0, -1, -1, -1, 0, 1, 1, 1, 0, -1};

  int sxy = sigma[x][y]; // the central site
  if (sxy == 0)
    return true;

  int n_borders = 0; // to count the amount of sites in state sxy bordering a site !=sxy

  static int stack[8]; // stack to count number of different surrounding cells
  int stackp = -1;
  bool one_of_neighbors_medium = false;

  for (int i = 1; i <= 8; i++)
  {
    int xcn = x + cyc_nx[i];
    int ycn = y + cyc_ny[i];
    int xncn = x + cyc_nx[i + 1];
    int yncn = y + cyc_ny[i + 1];

    if (par.periodic_boundaries)
    {
      if (xcn <= 0)
        xcn = sizex - 2 + xcn;
      if (ycn <= 0)
        ycn = sizey - 2 + ycn;
      if (xcn >= sizex - 1)
        xcn = xcn - sizex + 2;
      if (ycn >= sizey - 1)
        ycn = ycn - sizey + 2;
      if (xncn <= 0)
        xncn = sizex - 2 + xncn;
      if (yncn <= 0)
        yncn = sizey - 2 + yncn;
      if (xncn >= sizex - 1)
        xncn = xncn - sizex + 2;
      if (yncn >= sizey - 1)
        yncn = yncn - sizey + 2;
    }

    int s_nb = sigma[xcn][ycn];
    int s_next_nb = sigma[xncn][yncn];

    if ((s_nb > 0 || s_next_nb > 0) && (s_nb == 0 || s_next_nb == 0))
    {

      // check whether s_nb is adjacent to non-identical site,
      // count it
      n_borders++;
    }
    int j;
    bool on_stack_p = false;
    // we need the next heuristic to prevent stalling at
    // cell-cell borders
    // do not enforce constraint at two cell interface(no medium)
    if (s_nb)
    {
      for (j = stackp; j >= 0; j--)
      {
        if (s_nb == stack[j])
        {
          on_stack_p = true;
          break;
        }
      }
      if (!on_stack_p)
      {
        if (stackp > 6)
        {
          cerr << "Stack overflow, stackp=" << stackp << "\n";
        }
        stack[++stackp] = s_nb;
      }
    }
    else
    {
      one_of_neighbors_medium = true;
    }
  }

  // number of different neighbours is stackp+1;
  if (n_borders > 2 && ((stackp + 1) > 2 || one_of_neighbors_medium))
  {
    return false;
  }
  else
    return true;
}

double CellularPotts::CellDensity(void) const
{
  // return the density of cells
  int sum = 0;
  for (int i = 0; i < sizex * sizey; i++)
  {
    if (sigma[0][i])
    {
      sum++;
    }
  }
  return (double)sum / (double)(sizex * sizey);
}

double CellularPotts::MeanCellArea(void) const
{
  int sum_area = 0, n = 0;
  double sum_length = 0.;
  vector<Cell>::iterator c = cell->begin();
  ++c;

  for (;
       c != cell->end();
       c++)
  {
    sum_area += c->Area();
    sum_length += c->Length();
    n++;
  }
  cerr << "Mean cell length is " << sum_length / ((double)n) << endl;
  return (double)sum_area / (double)n;
}

void CellularPotts::ResetTargetLengths(void)
{
  vector<Cell>::iterator c = cell->begin();
  ++c;
  for (; c != cell->end(); c++)
  {
    if (c->getTau() == 1){
      c->SetTargetLength(par.celltype1_length);
    }
    else if (c->getTau() == 2){
      c->SetTargetLength(par.celltype2_length);
    }
  }
}

void CellularPotts::SetRandomTypes(void)
{
  // each cell gets a random type 1..maxtau
  vector<Cell>::iterator c = cell->begin();
  ++c;
  for (; c != cell->end(); c++)
  {
    int celltype = RandomNumber(Cell::maxtau);
    // cerr << "Setting celltype " << celltype << endl;
    c->setTau(celltype);
  }
}

void CellularPotts::SetTypesWithMask(void)
{ 
  DetectSidesIsthmus();
  int xcoord;
  int celltype;

  vector<Cell>::iterator c = cell->begin();
  ++c;
  for (;
       c != cell->end();
       c++)
  { 

    int xcoord = c->getCenterX();
    if (xcoord > right_side_isthmus)
      celltype = 1;
    else
      celltype = 2;
    c->setTau(celltype);
  }
}

void CellularPotts::SetTypesWithDoubleMask(void)
{ 
  for (int x = 0; x < sizex; x++)
    for (int y = 0; y < sizey; y++){
      if (mask[x][y] == 1 || mask[x][y] == 2)
        (*cell)[sigma[x][y]].setTau(mask[x][y]);}
}


void CellularPotts::SetUpTauMatrix(int sizex, int sizey)
{
  for (int x = 0; x < sizex; x++)
  {
    for (int y = 0; y < sizey; y++)
    {
      if (sigma[x][y] == -1)
        tau[x][y] = -1;
      else
        tau[x][y] = (*cell)[sigma[x][y]].getTau();
    }
  }
}

void CellularPotts::GrowAndDivideCells(int growth_rate)
{
  vector<Cell>::iterator c = cell->begin();
  ++c;
  vector<bool> which_cells(cell->size());
  for (; c != cell->end(); c++)
  {
    // only tumor cells grow and divide
    if (c->getTau() == 2)
    {
      c->SetTargetArea(c->TargetArea() + growth_rate);

      if (c->Area() > par.target_area)
      {
        which_cells[c->Sigma()] = true;
      }
      else
      {
        which_cells[c->Sigma()] = false;
      }
      if (c->chem[1] < 0.9)
      { // arbitrary oxygen threshold for the moment
        c->setTau(3);
      }
    }
    else
    {
      which_cells[c->Sigma()] = false;
    }
  }
  DivideCells(which_cells);
}

double CellularPotts::DrawConvexHull(Graphics *g, int color)
{
  // Draw the convex hull of the cells
  // using Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1. Prepare data for 2D hull code

  // count number of points to determine size of array
  int np = 0;
  for (int x = 1; x < sizex - 1; x++)
    for (int y = 1; y < sizey - 1; y++)
    {
      if (sigma[x][y])
      {
        np++;
      }
    }

  Point *p = new Point[np];
  int pc = 0;
  for (int x = 1; x < sizex - 1; x++)
  {
    for (int y = 1; y < sizey - 1; y++)
    {
      if (sigma[x][y])
      {
        p[pc++] = Point(x, y);
      }
    }
  }
  // Step 2: call 2D Hull code
  Point *hull = new Point[np];
  int nph = chainHull_2D(p, np, hull);

  // Step 3: draw it
  for (int i = 0; i < nph - 1; i++)
  {
    g->Line(hull[i].x, hull[i].y, hull[i + 1].x, hull[i + 1].y, color);
  }

  // Step 4: calculate area of convex hull
  double hull_area = 0.;
  for (int i = 0; i < nph - 1; i++)
  {
    hull_area += hull[i].x * hull[i + 1].y - hull[i + 1].x * hull[i].y;
  }
  hull_area /= 2.;
  // cerr << "Area = " << hull_area << "\n";

  delete[] p;
  delete[] hull;
  return hull_area;
}

void CellularPotts::CropSurface(int* bounds){
  DetectSidesIsthmus();
  int top = 0;
  int bottom = sizey;
  int left = sizex;
  int right = 0;
  for (int y = 0; y < sizey; y++){
    if (mask[right_side_isthmus-1][y]){
      if (y > top)
        top = y;
      if (y < bottom)
        bottom = y;
    }
    for (int x = 0; x < sizex; x++){
      if (tau[x][y] == 2 && x > right)
        right = x;
      if (tau[x][y] == 1 && x < left)
        left = x;
    }
  }
  bounds[0] = bottom;
  bounds[1] = top;
  bounds[2] = left;
  bounds[3] = right;
}

double CellularPotts::Convexity(void){
  //Compute the convexity (or concavity) at the surface between the two cell types
  int bounds[4];
  bounds[1] = 1;
  CropSurface(bounds);
  //for (int i = 0; i < 4; i++)
  //  cout << "bounds[" << i << "] = " << bounds[i] << endl;
  double compactness_1 = Compactness(bounds, 1);
  double compactness_2 = Compactness(bounds, 2);
  //cout << "Compactness celltype 1 = " << compactness_1 << endl;
  //cout << "Compactness celltype 2 = " << compactness_2 << endl;
  double convexity = compactness_2-compactness_1;
  //cout << "Convexity = " << convexity << endl;
  return convexity;
}


double CellularPotts::Compactness(int *bounds, int celltype)
{
  // Calculate compactness using the convex hull of the cells, including the corner points of pixels
  // We use Andrew's Monotone Chain Algorithm (see hull.cpp)

  // Step 1: calculate total cell area

  double cell_area = 0;
  for (int x = bounds[2]; x < bounds[3]+1; x++) //count only within the box
    for (int y = bounds[0]; y < bounds[1]+1; y++)
    {
      if (tau[x][y] == celltype) //Only consider one celltype
      {
       cell_area++;
      }
    }

  int np = 0;
  // Step 2. Prepare data for 2D hull code

  // Step 2a. Count number of corner points to determine size of array
  
  //First consider the left-most column, a corner point if there is a pixel (or to the bottom of it)
  if (tau[bounds[2]][bounds[0]] == celltype) //bottom row separately
    np++;
  for (int y = bounds[0]+1; y < bounds[1]+1; y++)
    if (tau[bounds[2]][y] == celltype || tau[bounds[2]][y-1] == celltype) //Only consider one celltype
      //add corner point only if a pixel is present at this location or below it.
      {
        np++;
      }
  if (tau[bounds[2]][bounds[1]] == celltype) //add top-left most corner points if there is a pixel there
    np++;


  //Add all 'inner' corner points
  for (int x = bounds[2]+1; x < bounds[3]+1; x++)
  {
    if (tau[x][bounds[0]] == celltype || tau[x-1][bounds[0]] == celltype)
      np++; //special case for bottom row is required
    for (int y = bounds[0]+1; y < bounds[1]+1; y++) //loop over all other rows
    {
      if (tau[x][y] == celltype || tau[x-1][y] == celltype || tau[x][y-1] == celltype ||tau[x-1][y-1] == celltype) //Only consider one celltype
      //and add a corner point on the bottom left of the current pixel if one of the adjacent pixels is present
      {
        np++;
      }
    }
    if (tau[x][bounds[1]] == celltype || tau[x-1][bounds[1]] == celltype)
      np++;
    //add the top-most corner point only if there is a pixel on the top row (or to the left of this pixel)
  }

  //Consider the right-most column separately, only add a corner point if there is a pixel in this column (or to the bottom of it)
  if (tau[bounds[3]][bounds[0]] == celltype) //bottom row separately
    np++;
  for (int y = bounds[0]+1; y < bounds[1]+1; y++)
    if (tau[bounds[3]][y] == celltype || tau[bounds[3]][y-1] == celltype) //Only consider one celltype
      //add corner point only if a pixel is present at this location or below it.
      {
        np++;
      }
  if (tau[bounds[3]][bounds[1]] == celltype)
    np++;

  // Step 2b. Create array which will contain all corner points

  
  Point *p = new Point[np];

  // Step 2c. Fill the array with all lattice points ordered, x-first.
  int pc = 0;

  //First consider the left-most column, a corner point if there is a pixel (or to the bottom of it)
  if (tau[bounds[2]][bounds[0]] == celltype) //bottom row separately
    p[pc++] = Point(bounds[2]-0.5, bounds[0]-0.5); 
  for (int y = bounds[0]+1; y < bounds[1]+1; y++)
    if (tau[bounds[2]][y] == celltype || tau[bounds[2]][y-1] == celltype) //Only consider one celltype
      //add corner point only if a pixel is present at this location or below it.
      {
        p[pc++] = Point(bounds[2]-0.5, y-0.5);
      }
  if (tau[bounds[2]][bounds[1]] == celltype) //add top-left most corner points if there is a pixel there
    p[pc++] = Point(bounds[2]-0.5, bounds[1]+0.5); 


  //Add all 'inner' corner points
  for (int x = bounds[2]+1; x < bounds[3]+1; x++)
  {
    if (tau[x][bounds[0]] == celltype || tau[x-1][bounds[0]] == celltype)
      p[pc++] = Point(x-0.5, bounds[0]-0.5); //special case for bottom row is required 
    for (int y = bounds[0]+1; y < bounds[1]+1; y++) //loop over all other rows
    {
      if (tau[x][y] == celltype || tau[x-1][y] == celltype || tau[x][y-1] == celltype ||tau[x-1][y-1] == celltype) //Only consider one celltype
      //and add a corner point on the bottom left of the current pixel if one of the adjacent pixels is present
      {
        p[pc++] = Point(x-0.5, y-0.5); 
      }
    }
    if (tau[x][bounds[1]] == celltype || tau[x-1][bounds[1]] == celltype)
      p[pc++] = Point(x-0.5, bounds[1]+0.5);
    //add the top-most corner point only if there is a pixel on the top row (or to the left of this pixel)
  }

  //Consider the right-most column separately, only add a corner point if there is a pixel in this column (or to the bottom of it)
  if (tau[bounds[3]][bounds[0]] == celltype) //bottom row separately
    p[pc++] = Point(bounds[3]+0.5, bounds[0]-0.5);
  for (int y = bounds[0]+1; y < bounds[1]+1; y++)
    if (tau[bounds[3]][y] == celltype || tau[bounds[3]][y-1] == celltype) //Only consider one celltype
      //add corner point only if a pixel is present at this location or below it.
      {
        p[pc++] = Point(bounds[3]+0.5, y-0.5);
      }
  if (tau[bounds[3]][bounds[1]] == celltype)
    p[pc++] = Point(bounds[3]+0.5, bounds[1]+0.5);

  // Step 3: call 2D Hull code
  Point *hull = new Point[np];
  int nph = chainHull_2D(p, np, hull);

  // Step 4: calculate area of convex hull

  double hull_area = 0.;
  for (int i = 0; i < nph - 1; i++)
  { 
    hull_area += hull[i].x * hull[i + 1].y - hull[i + 1].x * hull[i].y;
  }
  hull_area /= 2.;

  delete[] p;
  delete[] hull;

  // return compactness
  //cout << "cell_area = " << cell_area << endl;
  //cout << "hull_area = " << hull_area << endl;
  return cell_area / hull_area;
}


void CellularPotts::SetBoundingBox(void)
{
  int min_x = sizex + 2, max_x = 0;
  int min_y = sizey + 2, max_y = 0;
  for (int x = 1; x <= sizex - 2; x++)
  {
    for (int y = 1; y <= sizey - 2; y++)
    {
      if (sigma[x][y])
      {
        if (x < min_x)
        {
          min_x = x;
        }
        if (x > max_x)
        {
          max_x = x;
        }
        if (y < min_y)
        {
          min_y = y;
        }
        if (y > max_y)
        {
          max_y = y;
        }
      }
    }
  }
}

int CellularPotts::BoundaryLength(int start_x, int start_y, int end_x, int end_y){
  cout << "Start boundary length" << endl;
  int loc_x = start_x;
  int loc_y = start_y;
  int MaxBoundaryLength = 0;
  int BoundaryLength = 0;
  char orientation = 'W';
  while (!(loc_x == end_x && loc_y == end_y)){
    
    /*
    cout << "BL = " << BoundaryLength << endl;
    cout << "Max_BL " << MaxBoundaryLength << endl;
    cout << "Orientation = " << orientation << endl;  
    cout << "(" << loc_x << ", " << loc_y << ")" << endl;
    cout << "Orientation = " << orientation << endl;
    cout << "Cell type at (" << loc_x << ", " << loc_y << ") = " << (*cell)[sigma[loc_x][loc_y]].getTau() <<  " and sigma = " << sigma[loc_x][loc_y] << endl;
    cout << "Cell type at (" << loc_x-1 << ", " << loc_y << ") = " << (*cell)[sigma[loc_x-1][loc_y]].getTau() <<  " and sigma = " << sigma[loc_x-1][loc_y] << endl;
    cout << "Cell type at (" << loc_x << ", " << loc_y-1 << ") = " << (*cell)[sigma[loc_x][loc_y-1]].getTau()  <<  " and sigma = " << sigma[loc_x][loc_y-1] << endl;
    cout << "Cell type at (" << loc_x-1 << ", " << loc_y-1 << ") = " << (*cell)[sigma[loc_x-1][loc_y-1]].getTau()  <<  " and sigma = " << sigma[loc_x-1][loc_y-1] << endl;
    */
    
    //North
    if (orientation == 'N'){
      if (sigma[loc_x-1][loc_y] != -1 && (*cell)[sigma[loc_x-1][loc_y]].getTau() == 1){ //if cell to the left belongs to celltype 1
        if (sigma[loc_x-1][loc_y+1] != -1 && (*cell)[sigma[loc_x-1][loc_y+1]].getTau() == 1){ //if cell left forward also belongs to celltype 1
          orientation = 'E'; //Rotate clockwise
          if (sigma[loc_x][loc_y+1] != -1 && (*cell)[sigma[loc_x][loc_y+1]].getTau() == 2){
            BoundaryLength++; //If the pixel we were facing belonged to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
          loc_x = loc_x-1;
          loc_y = loc_y+1;
        }
        else{// move to the left instead
          loc_x = loc_x-1;
          loc_y = loc_y;
          if (sigma[loc_x][loc_y+1] != -1 && (*cell)[sigma[loc_x][loc_y+1]].getTau() == 2){
            BoundaryLength++; //If the new pixel we are facing belongs to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
        }

      }
      else{ // If pixel to the left does not belong to celltype 1
        orientation = 'W'; //Rotate counterclockwise
        if (sigma[loc_x-1][loc_y] != -1 && (*cell)[sigma[loc_x-1][loc_y]].getTau() == 2){ //If new cell we face belongs to celltype 2, increase boundary length
          BoundaryLength++;
        }
        else{ //Else, reset boundary length, keep track of the longest one
          if (BoundaryLength > MaxBoundaryLength)
            MaxBoundaryLength = BoundaryLength;
          BoundaryLength = 0;
        }
      } 
    }

    //East
    else if (orientation == 'E'){
      if (sigma[loc_x][loc_y+1] != -1 && (*cell)[sigma[loc_x][loc_y+1]].getTau() == 1){ //if cell to the left belongs to celltype 1
        if (sigma[loc_x+1][loc_y+1] != -1 && (*cell)[sigma[loc_x+1][loc_y+1]].getTau() == 1){ //if cell left forward also belongs to celltype 1
          orientation = 'S'; //Rotate clockwise
          if (sigma[loc_x+1][loc_y] != -1 && (*cell)[sigma[loc_x+1][loc_y]].getTau() == 2){
            BoundaryLength++; //If the pixel we were facing belonged to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
          loc_x = loc_x+1;
          loc_y = loc_y+1;
        }
        else{// move to the left instead
          loc_x = loc_x;
          loc_y = loc_y+1;
          if (sigma[loc_x+1][loc_y] != -1 && (*cell)[sigma[loc_x+1][loc_y]].getTau() == 2){
            BoundaryLength++; //If the new pixel we are facing belongs to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
        }
      }
      else{ // If pixel to the left does not belong to celltype 1
        orientation = 'N'; //Rotate counterclockwise
        if (sigma[loc_x][loc_y+1] != -1 && (*cell)[sigma[loc_x][loc_y+1]].getTau() == 2){ //If new cell we face belongs to celltype 2, increase boundary length
          BoundaryLength++;
        }
        else{ //Else, reset boundary length, keep track of the longest one
          if (BoundaryLength > MaxBoundaryLength)
            MaxBoundaryLength = BoundaryLength;
          BoundaryLength = 0;
        }
      } 
    }


    //South
    else if (orientation == 'S'){
      if (sigma[loc_x+1][loc_y] != -1 && (*cell)[sigma[loc_x+1][loc_y]].getTau() == 1){ //if cell to the left belongs to celltype 1
        if (sigma[loc_x+1][loc_y-1] != -1 && (*cell)[sigma[loc_x+1][loc_y-1]].getTau() == 1){ //if cell left forward also belongs to celltype 1
          orientation = 'W'; //Rotate clockwise
          if (sigma[loc_x][loc_y-1] != -1 && (*cell)[sigma[loc_x][loc_y-1]].getTau() == 2){
            BoundaryLength++; //If the pixel we were facing belonged to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
          loc_x = loc_x+1;
          loc_y = loc_y-1;
        }
        else{// move to the left instead
          loc_x = loc_x+1;
          loc_y = loc_y;
          if (sigma[loc_x][loc_y-1] != -1 && (*cell)[sigma[loc_x][loc_y-1]].getTau() == 2){
            BoundaryLength++; //If the new pixel we are facing belongs to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
        }
      }
      else{ // If pixel to the left does not belong to celltype 1
        orientation = 'E'; //Rotate counterclockwise
        if (sigma[loc_x+1][loc_y] != -1 && (*cell)[sigma[loc_x+1][loc_y]].getTau() == 2){ //If new cell we face belongs to celltype 2, increase boundary length
          BoundaryLength++;
        }
        else{ //Else, reset boundary length, keep track of the longest one
          if (BoundaryLength > MaxBoundaryLength)
            MaxBoundaryLength = BoundaryLength;
          BoundaryLength = 0;
        } 
      } 
    }


    //West
    else if (orientation == 'W'){
      if (sigma[loc_x][loc_y-1] != -1 && (*cell)[sigma[loc_x][loc_y-1]].getTau() == 1){ //if cell to the left belongs to celltype 1
        if (sigma[loc_x-1][loc_y-1] != -1 && (*cell)[sigma[loc_x-1][loc_y-1]].getTau() == 1){ //if cell left forward also belongs to celltype 1
          orientation = 'N'; //Rotate clockwise
          if (sigma[loc_x-1][loc_y] != -1 && (*cell)[sigma[loc_x-1][loc_y]].getTau() == 2){
            BoundaryLength++; //If the pixel we were facing belonged to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
          loc_x = loc_x-1;
          loc_y = loc_y-1;
        }
        else{// move to the left instead
          loc_x = loc_x;
          loc_y = loc_y-1;
          if (sigma[loc_x-1][loc_y] != -1 && (*cell)[sigma[loc_x-1][loc_y]].getTau() == 2){
            BoundaryLength++; //If the new pixel we are facing belongs to celltype 2, increase boundary length
          }
          else{ //Else, reset boundary length, keep track of the longest one
            if (BoundaryLength > MaxBoundaryLength)
              MaxBoundaryLength = BoundaryLength;
            BoundaryLength = 0;
          }
        }
      }
      else{ // If pixel to the left does not belong to celltype 1
        orientation = 'S'; //Rotate counterclockwise
        if (sigma[loc_x][loc_y-1] != -1 && (*cell)[sigma[loc_x][loc_y-1]].getTau() == 2){ //If new cell we face belongs to celltype 2, increase boundary length
          BoundaryLength++;
        }
        else{ //Else, reset boundary length, keep track of the longest one
          if (BoundaryLength > MaxBoundaryLength)
            MaxBoundaryLength = BoundaryLength;
          BoundaryLength = 0;
        } 
      } 
    }
  }
  return MaxBoundaryLength;
}



// useful to demonstrate large q-Potts
void CellularPotts::RandomSigma(int n_cells)
{
  for (int x = 0; x < sizex; x++)
  {
    for (int y = 0; y < sizey; y++)
    {
      sigma[x][y] = (int)(n_cells * RANDOM());
      tau[x][y] = (*cell)[sigma[x][y]].getTau();
    }
  }
}

bool CellularPotts::plotPos(int x, int y, Graphics *graphics)
{
  int self = sigma[x][y];
  if (par.micropatternmask == string("None"))
  {
    if (self == 0)
      return true;
    graphics->Rectangle((*cell)[self].Colour(), x, y);
    return false;
  }

  else
  {
    /*    if (numberofedges[x][y])
          graphics->Rectangle((*cell)[self].Colour()+numberofedges[x][y], x, y);*/
    // uncomment to view edge boundaries clearly
    if (mask[x][y] == false)
    {
      graphics->Rectangle(11, x, y);
      return false;
    }
    else if (self == 0)
      return true;
    else
    {
      float opacity = 0.1; //0 is fully tranparent, 1 fully opaque
      graphics->Rectangle((*cell)[self].Colour(), x, y, opacity);
      return false;
    }
  }
}

void CellularPotts::linePlotPos(int x, int y, Graphics *graphics)
{
  int self = sigma[x][y];
  int a = self, b = self, c = self, d = self;
  if (x != 0)
    a = sigma[x - 1][y];
  if (y != 0)
    b = sigma[x][y - 1];
  if (x != par.sizex - 1)
    c = sigma[x + 1][y];
  if (y != par.sizey - 1)
    d = sigma[x][y + 1];
  if (self != a)
    graphics->Line(x, y, x, y + 1, 1);
  if (self != b)
    graphics->Line(x, y, x + 1, y, 1);
  if (self != c)
    graphics->Line(x + 1, y, x + 1, y + 1, 1);
  if (self != d)
    graphics->Line(x, y + 1, x + 1, y + 1, 1);
}

void CellularPotts::anneal(int steps)
{
  for (int i = 0; i < steps; i++)
    AmoebaeMove(0, true);
}

int **CellularPotts::get_annealed_sigma(int steps)
{
  int **tmp_a = sigma;
  int **tmp_b;
  AllocateSigma(par.sizex, par.sizey);
  std::copy(*tmp_a, (*tmp_a) + (par.sizex * par.sizey), *sigma);
  anneal(steps);
  tmp_b = sigma;
  sigma = tmp_a;
  return tmp_b;
}

void CellularPotts::StoreMask()
{
  AllocateMask(sizex, sizey);
  // Create an input filestream
  std::ifstream myFile(par.micropatternmask);

  // Make sure the file is open
  if (!myFile.is_open())
    throw std::runtime_error("Could not open mask file");

  std::string line;
  int xcoord, ycoord, val;

  // Read data, line by line
  while (std::getline(myFile, line))
  {
    std::stringstream ss(line);
    int colnr = 0;
    // extract coordinates
    while (ss >> val)
    {
      // Add the current integer to the 'colIdx' column's values vector

      if (!colnr)
        xcoord = val;
      if (colnr)
        ycoord = val;
      colnr++;

      // If the next token is a comma, ignore it and move on
      if (ss.peek() == ',')
        ss.ignore();
    }
    mask[xcoord][ycoord] = 1;
  }

  if (par.second_layer){
    // Create an input filestream
    std::ifstream myFile(par.micropatternlayer2);

    // Make sure the file is open
    if (!myFile.is_open())
      throw std::runtime_error("Could not open mask file");

    // Read data, line by line
    while (std::getline(myFile, line))
    {
      std::stringstream ss(line);
      int colnr = 0;
      // extract coordinates
      while (ss >> val)
      {
        // Add the current integer to the 'colIdx' column's values vector

        if (!colnr)
          xcoord = val;
        if (colnr)
          ycoord = val;
        colnr++;

        // If the next token is a comma, ignore it and move on
        if (ss.peek() == ',')
          ss.ignore();
      }
      mask[xcoord][ycoord] = 2;
    }


  }

}
