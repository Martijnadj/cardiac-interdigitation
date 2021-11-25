void ComputeDerivs(PDEFIELD_TYPE t, PDEFIELD_TYPE* y, PDEFIELD_TYPE* dydt, int id){
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
  PDEFIELD_TYPE i_stim_Amplitude 		= 5.5e-10;//7.5e-10;   // ampere (in stim_mode)
  PDEFIELD_TYPE i_stim_End 				= 1000.0;   // second (in stim_mode)
  PDEFIELD_TYPE i_stim_PulseDuration	= 0.005;   // second (in stim_mode)
  PDEFIELD_TYPE i_stim_Start 			= 0.0;   // second (in stim_mode)
  PDEFIELD_TYPE i_stim_frequency        = 60.0;   // per_second (in stim_mode)
  //PDEFIELD_TYPE stim_flag 				= stimFlag;   // dimensionless (in stim_mode)
  PDEFIELD_TYPE i_stim_Period 			= 60.0/i_stim_frequency;

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
  /*
  //// Output variables
  IK1     = i_K1;
  Ito     = i_to;
  IKr     = i_Kr;
  IKs     = i_Ks;
  ICaL    = i_CaL;
  INaK    = i_NaK;
  INa     = i_Na;
  INaCa   = i_NaCa;
  IpCa    = i_PCa;
  If      = i_f;
  IbNa    = i_b_Na;
  IbCa    = i_b_Ca;
  Irel    = i_rel;
  Iup     = i_up;
  Ileak   = i_leak;
  Istim   = i_stim;
  INaL    = i_NaL; */

  /*dati = [INa, If, ICaL, Ito, IKs, IKr, IK1, INaCa, INaK, IpCa, IbNa, IbCa, Irel, Iup, Ileak, Istim, E_K, E_Na, INaL];



   //Equations for Noble1962 \
  PDEFIELD_TYPE C_m = 12;
  PDEFIELD_TYPE g_Na_bar = 400; 
  PDEFIELD_TYPE E_An = -60;
  PDEFIELD_TYPE g_An = 0;
    


  PDEFIELD_TYPE g_K1 = 1.2*exp((-y[0] - 90)/50) + 0.015*exp((y[0]+90)/60);      //eq 5
  PDEFIELD_TYPE g_K2 = 1.2*pow((float)y[1],4f);                                         //eq 6
  PDEFIELD_TYPE alpha_n;
  if (y[0] == -50){
    alpha_n = 0.0001*10;} 
  else{  
    alpha_n = (0.0001*(-y[0]-50))/(exp((-y[0]-50)/10)-1);}         //eq 8
  PDEFIELD_TYPE beta_n = 0.002*exp((-y[0]-90)/80);                              //eq 9
  PDEFIELD_TYPE I_K = (g_K1 + g_K2)*(y[0]+ 100);                                //eq 10
      
  PDEFIELD_TYPE g_Na = pow((float)y[2], 3.0f)*y[3]*g_Na_bar;                              //eq 12
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

  PDEFIELD_TYPE I_Na = (400*pow((float)y[2],3.0f)*y[3] + 0.14)*(y[0] - 40);               //eq 20
        
  PDEFIELD_TYPE I_An = g_An*(y[0]-E_An);                                        //eq 3    
  dydt[0] =  -(I_Na + I_K + I_An)/C_m;                            //eq 4
  dydt[1] = alpha_n*(1-y[1]) - beta_n*y[1];                       //eq 7
  dydt[2] = alpha_m*(1-y[2]) - beta_m*y[2];                       //eq 13
  dydt[3] = alpha_h*(1-y[3]) - beta_h*y[3];                       //eq 14 
  */

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
  /* //Noble1962
  PDEFIELD_TYPE dydt[4];
  PDEFIELD_TYPE ak2[4];
  PDEFIELD_TYPE ak3[4];
  PDEFIELD_TYPE ak4[4];
  PDEFIELD_TYPE ak5[4];
  PDEFIELD_TYPE ak6[4];
  PDEFIELD_TYPE ytemp[4];
  
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
  */

  //Paci2020
  PDEFIELD_TYPE dydt[23];
  PDEFIELD_TYPE ak2[23];
  PDEFIELD_TYPE ak3[23];
  PDEFIELD_TYPE ak4[23];
  PDEFIELD_TYPE ak5[23];
  PDEFIELD_TYPE ak6[23];
  PDEFIELD_TYPE ytemp[23];  
  ComputeDerivs(t,y,dydt,id);
  for (i=0;i<23;i++) //First step.
    ytemp[i]=y[i]+b21*dt*dydt[i];
  ComputeDerivs(t+a2*dt,ytemp,ak2,id);// Second step.

  for (i=0;i<23;i++)
    ytemp[i]=y[i]+dt*(b31*dydt[i]+b32*ak2[i]);
  ComputeDerivs(t+a3*dt,ytemp,ak3,id); //Third step.
  for (i=0;i<23;i++)
    ytemp[i]=y[i]+dt*(b41*dydt[i]+b42*ak2[i]+b43*ak3[i]);
  ComputeDerivs(t+a4*dt,ytemp,ak4,id); //Fourth step.
  for (i=0;i<23;i++)
    ytemp[i]=y[i]+dt*(b51*dydt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  ComputeDerivs(t+a5*dt,ytemp,ak5,id); //Fifth step.
  for (i=0;i<23;i++)
    ytemp[i]=y[i]+dt*(b61*dydt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  ComputeDerivs(t+a6*dt,ytemp,ak6,id); //Sixth step.
  for (i=0;i<23;i++) //Accumulate increments with proper weights.
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
global const PDEFIELD_TYPE* couplingcoefficient,
int PDEsteps) {
  PDEFIELD_TYPE thetime = PDEsteps * dt;

  //ID is used for position in array
  int id = get_global_id(0); 
  //Calculate position in aray
  int layersize = xsize * ysize;
  int zpos = id / layersize;
  int xpos = (id - (zpos * layersize)) / ysize;
  int ypos = id - xpos * ysize - zpos * layersize; 
  //printf("id: %i x: %i y: %i\n", id, xpos, ypos);  
    
  
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
    //Paci2020
    int nvar = 23;

    int i;
    PDEFIELD_TYPE currentvalues[23];

    for (i = 0; i<nvar; i++)
      currentvalues[i] = sigmaA[id + i*layersize];

    sum += couplingcoefficient[id-1]*sigmaA[id-1];
    sum += couplingcoefficient[id+1]*sigmaA[id+1];
    sum += couplingcoefficient[id-ysize] * sigmaA[id-ysize];
    sum += couplingcoefficient[id+ysize] * sigmaA[id+ysize];
    sum-=(couplingcoefficient[id-1] + couplingcoefficient[id+1] + couplingcoefficient[id-ysize] + couplingcoefficient[id+ysize])*sigmaA[id];
    PDEFIELD_TYPE derivs[23];
    if (sigmacells[id] > 0){
      ComputeDerivs(thetime, currentvalues, derivs, id);
      RungeKutta(thetime, derivs, currentvalues, dt, id);
  
        //Diffusion

      sigmaB[id] = currentvalues[0]+derivs[0] + sum*dt/dx2;
      for (i = 1; i<nvar; i++)
        sigmaB[id + i*layersize] = currentvalues[i] + derivs[i];
      
    }
    else{

      sigmaB[id] = currentvalues[0]+ sum*dt/dx2;

      for (i = 1; i<nvar; i++)
        sigmaB[id + i*layersize] = currentvalues[i];

    }
    //sigmaB[id] = sigmaB[id]+0.1;



    if(fmod(thetime, 750) < -50 && xpos > 120 && xpos < 126 && ypos > 175 && ypos < 326){
      sigmaB[id] = 0;
    }
    /*Noble1962 
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
  */


    if (id == 189734)
      printf("Time = %.5f, V = %.5f \n", thetime, sigmaB[id]);
    if (!((sigmaB[id] > -1) && (sigmaB[id] < 1)) && thetime == 0.03206){
      printf("id = %i \n", id);
    }  

  }
}

