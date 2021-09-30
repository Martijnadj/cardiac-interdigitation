PDEFIELD_TYPE* ComputeDerivs(PDEFIELD_TYPE* y){
  PDEFIELD_TYPE C_m = 12;
  PDEFIELD_TYPE g_Na_bar = 400; 
  PDEFIELD_TYPE E_An = -60;
  PDEFIELD_TYPE g_An = 0;
    
  
  PDEFIELD_TYPE derivs[4];
  PDEFIELD_TYPE* q = derivs;

  PDEFIELD_TYPE g_K1 = 1.2*exp((-y[0] - 90)/50) + 0.015*exp((y[0]+90)/60);      //eq 5
  PDEFIELD_TYPE g_K2 = 1.2*pow(y[1],4);                                         //eq 6
  PDEFIELD_TYPE alpha_n;
  if (y[0] == -50){
    alpha_n = 0.0001;} 
  else{  
    alpha_n = (0.0001*(-y[0]-50))/(exp((-y[0]-50)/10)-1);}         //eq 8
  PDEFIELD_TYPE beta_n = 0.002*exp((-y[0]-90)/80);                              //eq 9
  PDEFIELD_TYPE I_K = (g_K1 + g_K2)*(y[0]+ 100);                                //eq 10
      
  PDEFIELD_TYPE g_Na = pow(y[2], 3)*y[3]*g_Na_bar;                              //eq 12
  PDEFIELD_TYPE alpha_h = 0.17*exp((-y[0]-90)/20);                              //eq 16
  PDEFIELD_TYPE beta_h = 1/(exp((-y[0]-42)/10) +1);
  PDEFIELD_TYPE alpha_m;
  if (y[0] == -48){
    alpha_m = 0.1;} 
  else{                                                                        //eq 17
    alpha_m = 0.1*(-y[0]-48)/(exp((-y[0]-48)/15) -1);}            //eq 18
  PDEFIELD_TYPE beta_m = 0.12*(y[0]+8)/(exp((y[0]+8)/5)-1);                     //eq 19
  PDEFIELD_TYPE I_Na = (400*pow(y[2],3)*y[3] + 0.14)*(y[0] - 40);               //eq 20
        
  PDEFIELD_TYPE I_An = g_An*(y[0]-E_An);                                        //eq 3    
  derivs[0] =  -(I_Na + I_K + I_An)/C_m;                            //eq 4
  derivs[1] = alpha_n*(1-y[1]) - beta_n*y[1];                       //eq 7
  derivs[2] = alpha_m*(1-y[2]) - beta_m*y[2];                       //eq 13
  derivs[3] = alpha_h*(1-y[3]) - beta_h*y[3];                       //eq 14
  return q;
}



void kernel SecreteAndDiffuse(
global const int* sigmacells,  
global const PDEFIELD_TYPE* sigmaA,  
global PDEFIELD_TYPE* sigmaB, 
int xsize, 
int ysize,
int layers,
PDEFIELD_TYPE decay_rate, 
PDEFIELD_TYPE dt, 
PDEFIELD_TYPE dx2, 
global const PDEFIELD_TYPE* diff_coeff,
PDEFIELD_TYPE secr_rate,
int btype) {



  //ID is used for position in array
  int id = get_global_id(0); 
  //printf("id = %i \n", id);
  // test b 
  //Calculate position in aray
  int layersize = xsize * ysize;
  int zpos = id / layersize;
  int xpos = (id - (zpos * layersize)) / ysize;
  int ypos = id - xpos * ysize - zpos * layersize; 
  //printf("id: %i x: %i y: %i\n", id, xpos, ypos);  
  if (layers == 1){
    //Boundaries
    PDEFIELD_TYPE sum =0;
    if (xpos == 0 || ypos == 0 || xpos == xsize-1 || ypos == ysize-1){
      switch(btype){
      case 1:
      //Noflux gradient
      if (ypos == ysize-1) sigmaB[id] = sigmaA[id-1]; 
      if (ypos == 0) sigmaB[id] = sigmaA[id+1];
      if (xpos == xsize-1) sigmaB[id] = sigmaA[id-ysize];
      if (xpos == 0) sigmaB[id] = sigmaA[id+ysize];
      break;
      //Periodic
      case 2:
      if (xpos == xsize-1) sigmaB[id] = sigmaA[id-layersize+ysize];
      if (xpos == 0) sigmaB[id] = sigmaA[id+layersize-ysize];
      if (ypos == ysize-1) sigmaB[id] = sigmaA[id-ysize+1];
      if (ysize == 0) sigmaB[id] = sigmaB[id+ysize-1]; 
      break;
      //Absorbing
      case 3:
      sigmaB[id] = 0;
      break;
      } 
    }
    else {
      //Retrieve current value in array
      PDEFIELD_TYPE value = sigmaA[id];
      //Secretion
      if (btype != 1){
        if (zpos == 0){
          if (sigmacells[id] > 0){
            value = value + secr_rate * dt;
          }
          else {
            value = value - decay_rate*dt*value;
          }
        }
      }
      //Diffusion
      sum += sigmaA[id-1];
      sum += sigmaA[id+1];
      sum += sigmaA[id-ysize];
      sum += sigmaA[id+ysize];
      sum-=4*value;
      sigmaB[id] = value+sum*dt*diff_coeff[zpos]/dx2;
      //sigmaB[id] =  value;
    }
  }

  else{
    
    
    
    //Boundaries
    PDEFIELD_TYPE sum =0;
    if (xpos == 0 || ypos == 0 || xpos == xsize-1 || ypos == ysize-1){
      switch(btype){
      case 1:
      //Noflux gradient
      if (ypos == ysize-1) sigmaB[id] = sigmaA[id-1]; 
      if (ypos == 0) sigmaB[id] = sigmaA[id+1];
      if (xpos == xsize-1) sigmaB[id] = sigmaA[id-ysize];
      if (xpos == 0) sigmaB[id] = sigmaA[id+ysize];
      break;
      //Periodic
      case 2:
      if (xpos == xsize-1) sigmaB[id] = sigmaA[id-layersize+ysize];
      if (xpos == 0) sigmaB[id] = sigmaA[id+layersize-ysize];
      if (ypos == ysize-1) sigmaB[id] = sigmaA[id-ysize+1];
      if (ysize == 0) sigmaB[id] = sigmaB[id+ysize-1]; 
      break;
      //Absorbing
      case 3:
      sigmaB[id] = 0;
      break;
      } 
    }
    
    else {
      PDEFIELD_TYPE currentvalues[4];
      PDEFIELD_TYPE* p= currentvalues;

      currentvalues[0] = sigmaA[id];
      currentvalues[1] = sigmaA[id + layersize];
      currentvalues[2] = sigmaA[id + 2*layersize];
      currentvalues[3] = sigmaA[id + 3*layersize];
      PDEFIELD_TYPE* derivs = ComputeDerivs(p);
      

      


      //Diffusion
      sum += sigmaA[id-1];
      sum += sigmaA[id+1];
      sum += sigmaA[id-ysize];
      sum += sigmaA[id+ysize];
      sum-=4*currentvalues[0];
      //sum = 0;
      sigmaB[id] = currentvalues[0]+derivs[0]*dt + sum*dt*diff_coeff[0]/dx2;
      
      //if (id < 1000000) printf("id = %i, sigmaB[id] = %.5f --- ", id, sigmaB[id]);

      
      sigmaB[id + 1*layersize] = currentvalues[1] + derivs[1]*dt;
      sigmaB[id + 2*layersize] = currentvalues[2] + derivs[2]*dt;
      sigmaB[id + 3*layersize] = currentvalues[3] + derivs[3]*dt;
    }
    if (id == 13014) printf("V = %.5f ; n = %.2f ; m = %.2f ; h = %.2f \n", sigmaB[id], sigmaB[id + layersize], sigmaB[id + 2*layersize], sigmaB[id + 3*layersize]);
    
  }
  //13485 and 13014
  //printf("id2 = %i \n", id);
}

