int GetLocation(int x, int y, int ysize){
  return x*ysize+y;
}

PDEFIELD_TYPE* ComputeDerivs(PDEFIELD_TYPE* y, int id){
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
    alpha_n = 0.0001*10;} 
  else{  
    alpha_n = (0.0001*(-y[0]-50))/(exp((-y[0]-50)/10)-1);}         //eq 8
  PDEFIELD_TYPE beta_n = 0.002*exp((-y[0]-90)/80);                              //eq 9
  PDEFIELD_TYPE I_K = (g_K1 + g_K2)*(y[0]+ 100);                                //eq 10
      
  PDEFIELD_TYPE g_Na = pow(y[2], 3)*y[3]*g_Na_bar;                              //eq 12
  PDEFIELD_TYPE alpha_h = 0.17*exp((-y[0]-90)/20);                              //eq 16
  PDEFIELD_TYPE beta_h = 1/(exp((-y[0]-42)/10) +1);                             //eq 17
  PDEFIELD_TYPE alpha_m;
  if (y[0] == -48){
    alpha_m = 0.1*15;} 
  else{                                                                        
    alpha_m = 0.1*(-y[0]-48)/(exp((-y[0]-48)/15) -1);}                          //eq 18  

  PDEFIELD_TYPE beta_m;
  if (y[0] == -8){
    beta_m = 0.12*5;} 
  else{                                                                        
    beta_m = 0.12*(y[0]+8)/(exp((y[0]+8)/5)-1);}                          //eq 19  

  PDEFIELD_TYPE I_Na = (400*pow(y[2],3)*y[3] + 0.14)*(y[0] - 40);               //eq 20
        
  PDEFIELD_TYPE I_An = g_An*(y[0]-E_An);                                        //eq 3    
  derivs[0] =  -(I_Na + I_K + I_An)/C_m;                            //eq 4
  derivs[1] = alpha_n*(1-y[1]) - beta_n*y[1];                       //eq 7
  derivs[2] = alpha_m*(1-y[2]) - beta_m*y[2];                       //eq 13
  derivs[3] = alpha_h*(1-y[3]) - beta_h*y[3];                       //eq 14

  //if(id == 190720) printf("I_Na = %.2f, I_K = %.2f, I_An = %.2f \n", I_Na, I_K, I_An);

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
int btype,
global const int* numberofedges,
int PDEsteps) {

  float thetime = PDEsteps * dt;

  //ID is used for position in array
  int id = get_global_id(0); 
  if(fmod(thetime, 250) < 1 && id == 1000);
  //printf("id = %i \n", id);
  // test b 
  //Calculate position in aray
  int layersize = xsize * ysize;
  int zpos = id / layersize;
  int xpos = (id - (zpos * layersize)) / ysize;
  int ypos = id - xpos * ysize - zpos * layersize; 
  //printf("id: %i x: %i y: %i\n", id, xpos, ypos);  
    
  //if (id == 13014) printf("Old:  V = %.5f ; n = %.2f ; m = %.2f ; h = %.2f \n", sigmaA[id], sigmaA[id + layersize], sigmaA[id + 2*layersize], sigmaA[id + 3*layersize]);   
  
  //Boundaries
  PDEFIELD_TYPE sum = 0;
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

    currentvalues[0] = sigmaA[id];
    currentvalues[1] = sigmaA[id + layersize];
    currentvalues[2] = sigmaA[id + 2*layersize];
    currentvalues[3] = sigmaA[id + 3*layersize];

    sum += sigmaA[id-1];
    sum += sigmaA[id+1];
    sum += sigmaA[id-ysize];
    sum += sigmaA[id+ysize];
    sum-=4*currentvalues[0];

    PDEFIELD_TYPE* derivs;
    if (numberofedges[id] > 0){
      derivs = ComputeDerivs(currentvalues, id);

      //Diffusion

      sigmaB[id] = currentvalues[0]+derivs[0]*dt + sum*dt*diff_coeff[0]/dx2;  
      sigmaB[id + 1*layersize] = currentvalues[1] + derivs[1]*dt;
      sigmaB[id + 2*layersize] = currentvalues[2] + derivs[2]*dt;
      sigmaB[id + 3*layersize] = currentvalues[3] + derivs[3]*dt;
    }
    else{

      sigmaB[id] = currentvalues[0]+ sum*dt*diff_coeff[0]/dx2;
      sigmaB[id + 1*layersize] = currentvalues[1];
      sigmaB[id + 2*layersize] = currentvalues[2];
      sigmaB[id + 3*layersize] = currentvalues[3];

    }
    int location;
    if(fmod(thetime, 750) < -50 && xpos > 120 && xpos < 126 && ypos > 175 && ypos < 326){
      sigmaB[id] = 0;
    }
  }
}

