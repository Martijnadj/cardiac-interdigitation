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

#ifndef _PDE_HH_
#define _PDE_HH_
#include <stdio.h>
#include <float.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <iostream>
#include <cusparse.h>

#include <MultiCellDS.hpp>
#include <MultiCellDS-pimpl.hpp>
#include <MultiCellDS-simpl.hpp>

#include "cl_manager.hpp"
#include "pdetype.h" 
#include "graph.hpp"

class CellularPotts;
class Dish;
class PDE {

 friend class Info;

 public:

  int sizex;
  int sizey;
  int layers;
  int btype;
  PDEFIELD_TYPE dt;
  PDEFIELD_TYPE min_stepsize;
  PDEFIELD_TYPE dx2;
  bool usePDEorAltPDE;

  /*! \brief Constructor for PDE object containing arbitrary number of planes.
  \param layers: Number of PDE planes
  \param sizex: horizontal size of PDE planes
  \param sizey: vertical size of PDE planes
  */
  PDE(const int layers, const int sizex, 
      const int sizey);
      
  // destructor must also be virtual
  virtual ~PDE();

  /*! \brief Plots one layer of the PDE plane to a Graphics window.
  \param g: Graphics window.
  \param layer: The PDE plane to be plotted. Default layer 0.
  */
  void Plot(Graphics *g, const int layer=0);
  /*! \brief Plots one layer of the PDE to a Graphics window, but not over the cells.
    \param g: Graphics window.
    \param cpm: CellularPotts object containing the cells.
    \param layer: The PDE plane to be plotted. Default layer 0.
  */
  void Plot(Graphics *g, CellularPotts *cpm, const int layer=0);
  
  /*! \brief Plots the PDE field using contour lines.
    
  \param g: Graphics window.
  \param layer: The PDE plane to be plotted. Default layer 0.
  \param colour: Color to use for the contour lines, as defined in the "default.ctb" color map file, which should be in the same directory as the executable. Default color 1 (black in the default color map).
  */
  void ContourPlot(Graphics *g, int layer=0, int colour=1);
  
  //! \brief Returns the horizontal size of the PDE planes.
  inline int SizeX() const {
    return sizex;
  }

  //! \brief Returns the vertical size of the PDE planes.
  inline int SizeY() const {
    return sizey;
  }

  //! \brief Returns the number of PDE layers in the PDE object
  inline int Layers() const {
    return layers;
  }
    
  //! \brief Set the \param name of the species in layer \param l
  void SetSpeciesName(int l, const char *name);
    
  /*! \brief Returns the value of grid point x,y of PDE plane "layer".
    
  Warning, no range checking done.
  
  \param layer: the PDE plane to probe.
  \param x, y: grid point to probe.
  */

  inline PDEFIELD_TYPE PDEVARS(const int layer, const int x, const int y) const {
    return PDEvars[layer*sizex*sizey+x*sizey+y];
  }

  /*! \brief Sets grid point x,y of PDE plane "layer" to value "value".
  \param layer: PDE plane.
  \param x, y: grid point
  \param value: new contents
  */

  
  /*! \brief Adds a number to a PDE grid point.
  \param layer: PDE plane.
  \param x, y: grid point
  \param value: value to add
  */


  /*! \brief Gets the maximum value of PDE layer l.
  \param l: layer
  \return Maximum value in layer l.
  */
  inline PDEFIELD_TYPE Max(int l) {
    PDEFIELD_TYPE max=PDEvars[l*sizex*sizey];
    int loop=sizex*sizey;
    for (int i=1;i<loop;i++)
      if (PDEvars[l*sizex*sizey+i]>max) {
	max=PDEvars[l*l*sizex*sizey+i];
      }
    return max;
  }
  /*! \brief Returns the minimum value of PDE layer l.
  \param l: layer
  \return Minimum value in layer l.
  */
  inline PDEFIELD_TYPE Min(int l) {
    PDEFIELD_TYPE min=PDEvars[l*sizex*sizey];
    int loop=sizex*sizey;
    for (int i=1;i<loop;i++)
      if (PDEvars[l*l*sizex*sizey+i]<min) {
	min=PDEvars[l*l*sizex*sizey+i];
      }
    return min;
  }

  
  /*! \brief Carry out $n$ diffusion steps for all PDE planes.
  We use a forward Euler method here. Can be replaced for better algorithm.
  Function for the Act model. The whole field is initialized, usually with 0
  */
  void InitializeAgeLayer(int l,double value,CellularPotts *cpm);
  void InitializePDEs(CellularPotts * cpm);
  void InitializeCuda(CellularPotts * cpm);
  void InitializePDEvars();

 /* Function for the Act model. All the lattice sites within cells are "aged"
	*  by decreasing their values, usually with 1.
	*/
  void AgeLayer(int l,double value,CellularPotts *cpm, Dish *dish);

  /* Function for the Act model. Plots the values of the activity into the cells.
  */
  void PlotInCells(Graphics *g,CellularPotts *cpm, const int l=0);
  // lymphocyte matrix interaction function

  void MILayerCA(int l,double value,CellularPotts *cpm, Dish *dish);
  /*! \brief Carry out $n$ diffusion steps for all PDE planes.

  We use a forward Euler method here. Can be replaced for better algorithm.

  \param repeat: Number of steps.

  Time step dt, space step dx, diffusion coefficient diff_coeff and
  boundary conditions (bool periodic_boundary) are set as global
  parameters in a parameter file using class Parameter.

  */
  void Diffuse(int repeat);

  /*! \brief Implementation of no-flux boundaries.
    
  Called internally (optionally) by Diffuse(). */
  void NoFluxBoundaries(void);
  
  /*! \brief Implementation of absorbing boundaries.
    
  Called internally (optionally) by Diffuse(). */
  void AbsorbingBoundaries(void);

  /*! \brief Implementation of periodic boundaries.
  Called internally (optionally) by Diffuse(). */
  void PeriodicBoundaries(void);

  /*! \brief Reaction and interaction of CPM plane with PDE planes.
  \param cpm: CellularPotts plane the PDE plane interacts with
  You should implement this member function as part of your main
  simulation code. See for an example vessel.cpp.
  */
  void Secrete(CellularPotts *cpm);

  //Secrete and diffuse functions accelerated using OpenCL
  void ODEstepCL(CellularPotts *cpm, int repeat);
  void cuPDEsteps(CellularPotts *cpm, int repeats);

  /*! \brief Returns cumulative "simulated" time,
    i.e. number of time steps * dt. */
  inline double TheTime(void) const {
    return thetime;
  }
  
  /*! \brief Returns summed amount of chemical in PDE plane "layer".
  \param layer: The PDE plane of which to sum the chemicals. layer=-1 (default) returns the summed amount of chemical in all planes.
  */
  double GetChemAmount(const int layer=-1);

  /*!   Calculates the first and second order gradients, i.e. gradx,
    grady, gradxx, gradxy and gradyy and puts them in the next
    three chemical fields. Not currently used and might need some
    redoing. Make sure you have allocated sufficient fields (this
    method generates five planes).

    \param layer: PDE plane of which to calculate the gradients
    (default 0) \param first_grad_layer: first plane of five in which
    to write the results (default 1).
  */
  void GradC(int layer=0, int first_grad_layer=1); 

  /*!   Plots a field of the first order gradients, i.e. gradx and
    grady; assumes you have called GradC before.
    Not currently used and might need some
    redoing. 
    \param g: Graphics window
    \param stride: Number of grid points between vectors (drawn as lines, currently.
    \param linelength: Length of vector lines, in pixels.
    \param first_grad_layer: first plane of two which contain the
    calculated gradients (default 1).
       
  */
  void PlotVectorField(Graphics &g, int stride, int linelength, int first_grad_layer=1);
  void InitLinearYGradient(int spec, double conc_top, double conc_bottom);
   
  bool plotPos(int x, int y, Graphics * graphics, int layer);

  void reset_plot(){ highest = Max(0); lowest = Min(0);}


  double highest;
  double lowest;

  protected:

  PDEFIELD_TYPE *PDEvars;
  PDEFIELD_TYPE *alt_PDEvars;
  PDEFIELD_TYPE *last_stepsize;
  PDEFIELD_TYPE *d_PDEvars;
  PDEFIELD_TYPE *d_alt_PDEvars;
  PDEFIELD_TYPE **couplingcoefficient;
  PDEFIELD_TYPE *d_couplingcoefficient;
  int **celltype;
  int *d_celltype;

  
  // Used as temporary memory in the diffusion step
  // (addresses will be swapped for every time step, so
  // never directly use them!!! Access is guaranteed to be correct
  // through user interface)

  

  //2d arrays containing the upperdiagonals, lower diagonals and diagonals and vectors B for all rows and columns to solve AX=B equations
  PDEFIELD_TYPE *lowerH, *upperH, *diagH, *BH, *lowerV, *upperV, *diagV, *BV, *XH, *next_stepsize; 

 
 
  // Protected member functions
  /*! \brief Used in Plot. Takes a color and turns it into a grey value.
  \param val: Value from PDE plane.
  Implement this function in you main simulation code. See e.g. vessel.cpp.
  */
  virtual int MapColour(double val);

  //virtual int MapColour3(double val, int l);

  //! empty constructor (necessary for derivation)
  PDE(void);
  
  /*! \brief Allocates a PDE plane (internal use). 
  For internal use, can be reimplemented in derived class to change
  method of memory allocation.
  */   
  //virtual PDEFIELD_TYPE ***AllocatePDEvars(const int layers, const int sx, const int sy);
  //virtual void AllocateTridiagonalvars(int sx, int sy);
  
  //void Tri_diag_inv(PDEFIELD_TYPE *du, PDEFIELD_TYPE *d, PDEFIELD_TYPE *dl, PDEFIELD_TYPE *B, PDEFIELD_TYPE *X, int size); 
  void Diffusionstep(PDEFIELD_TYPE ***PDEvars, CellularPotts *cpm);

  

 
private:
  PDEFIELD_TYPE z[10];
  
  static const int nx[9], ny[9];
  int PDEsteps;
  float thetime;

  inline double Z(double k, int steps);

  
  std::vector<std::string> species_names;
 
  void SetupOpenCL(); 

  
  //Store buffersize for horizontal / vertical ADI sweep


  size_t pbuffersizeH; 
  void *pbufferH;
  //Needed for cuSparse horizontal and vertical sweeps of ADI
  cusparseStatus_t statusH; 
  cusparseHandle_t handleH;
  size_t pbuffersizeV; 
  void *pbufferV;
  //Needed for cuSparse horizontal and vertical sweeps of ADI
  cusparseStatus_t statusV; 
  cusparseHandle_t handleV;



  //OpenCL variables
  bool openclsetup = false;
  cl::Program program;
  cl::Kernel kernel_ODEstep;
  bool first_round = true;
};

#endif
