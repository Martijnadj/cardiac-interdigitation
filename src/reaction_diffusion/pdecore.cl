void ComputeDerivs(PDEFIELD_TYPE t, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, int id){
  PDEFIELD_TYPE C_m = 12;
  PDEFIELD_TYPE g_Na_bar = 400; 
  PDEFIELD_TYPE E_An = -60;
  PDEFIELD_TYPE g_An = 0;
    


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
  dydt[0] =  -(I_Na + I_K + I_An)/C_m;                            //eq 4
  dydt[1] = alpha_n*(1-y[1]) - beta_n*y[1];                       //eq 7
  dydt[2] = alpha_m*(1-y[2]) - beta_m*y[2];                       //eq 13
  dydt[3] = alpha_h*(1-y[3]) - beta_h*y[3];                       //eq 14

  //if(id == 190720) printf("I_Na = %.2f, I_K = %.2f, I_An = %.2f \n", I_Na, I_K, I_An);
}


int GetLocation(int x, int y, int ysize){
  return x*ysize+y;
}

void RungeKutta(PDEFIELD_TYPE t, PDEFIELD_TYPE* dyoutdt, PDEFIELD_TYPE* y, PDEFIELD_TYPE dt, int id){
  

  int i;

  PDEFIELD_TYPE a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  PDEFIELD_TYPE dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  //PDEFIELD_TYPE *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
  PDEFIELD_TYPE dydt[4];
  PDEFIELD_TYPE ak2[4];
  PDEFIELD_TYPE ak3[4];
  PDEFIELD_TYPE ak4[4];
  PDEFIELD_TYPE ak5[4];
  PDEFIELD_TYPE ak6[4];
  PDEFIELD_TYPE ytemp[4];
  //PDEFIELD_TYPE dyoutdt[4];
  ComputeDerivs(t,y,dydt,id);
  for (i=0;i<4;i++) //First step.
    ytemp[i]=y[i]+b21*dt*dydt[i];
  ComputeDerivs(t+a2*dt,ytemp,ak2,id);// Second step.
  for (i=0;i<4;i++)
    ytemp[i]=y[i]+dt*(b31*dydt[i]+b32*ak2[i]);
  ComputeDerivs(t+a3*dt,ytemp,ak3,id); //Third step.
  for (i=0;i<4;i++)
    ytemp[i]=y[i]+dt*(b41*dydt[i]+b42*ak2[i]+b43*ak3[i]);
  ComputeDerivs(t+a4*dt,ytemp,ak4,id); //Fourth step.
  for (i=0;i<4;i++)
    ytemp[i]=y[i]+dt*(b51*dydt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  ComputeDerivs(t+a5*dt,ytemp,ak5,id); //Fifth step.
  for (i=0;i<4;i++)
    ytemp[i]=y[i]+dt*(b61*dydt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  ComputeDerivs(t+a6*dt,ytemp,ak6,id); //Sixth step.
  for (i=0;i<4;i++) //Accumulate increments with proper weights.
    dyoutdt[i]=dt*(c1*dydt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);

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
  PDEFIELD_TYPE thetime = PDEsteps * dt;

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

    PDEFIELD_TYPE derivs[4];
    if (numberofedges[id] > 0){
      ComputeDerivs(thetime, currentvalues, derivs, id);
      RungeKutta(thetime, derivs, currentvalues, dt, id);


      //Diffusion

      sigmaB[id] = currentvalues[0]+derivs[0] + sum*dt*diff_coeff[0]/dx2;  
      sigmaB[id + 1*layersize] = currentvalues[1] + derivs[1];
      sigmaB[id + 2*layersize] = currentvalues[2] + derivs[2];
      sigmaB[id + 3*layersize] = currentvalues[3] + derivs[3];
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

