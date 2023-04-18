#include <iostream>
#include <fstream>
#include <math.h>
#include <bits/stdc++.h> 
using namespace std;

#define PDEFIELD_FLOAT
//#define PDEFIELD_DOUBLE
 #ifdef PDEFIELD_FLOAT
#define PDEFIELD_TYPE float
#endif
#ifdef PDEFIELD_DOUBLE
#define PDEFIELD_TYPE double
#endif

void derivsMaleckar(PDEFIELD_TYPE VOI, PDEFIELD_TYPE* STATES,PDEFIELD_TYPE* RATES, PDEFIELD_TYPE pacing_interval, int id, PDEFIELD_TYPE activation_strength, PDEFIELD_TYPE t_A){


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


void ComputeQthr(PDEFIELD_TYPE* y, PDEFIELD_TYPE t_A, PDEFIELD_TYPE ddt, PDEFIELD_TYPE &Q_thr, int length){
  //cout << "t_A = " << t_A << endl;
  PDEFIELD_TYPE activation_strength = 50;
  PDEFIELD_TYPE stepsize_activation = activation_strength/2;
  PDEFIELD_TYPE elapsed_time = 0;
  PDEFIELD_TYPE dydt[length];
  PDEFIELD_TYPE y_copy[length];


  while(stepsize_activation > 1e-7){
    //cout << stepsize_activation << endl;
    for (int i = 0; i < length; i ++){
      y_copy[i] = y[i];
    }
    elapsed_time = 0;
    while (y_copy[0] < 0 && elapsed_time < 1){
      //derivsMaleckar_host(elapsed_time,y_copy,dydt,0, 0, -activation_strength, t_A); //linear
      derivsMaleckar(elapsed_time,y_copy,dydt,0, 0, -activation_strength, t_A);  //parabolic
      elapsed_time += ddt;
      for (int i = 0; i < length; i++) {
        y_copy[i]=y_copy[i]+ddt*dydt[i];  
      }
      //cout  << "y_copy[0] = " << y_copy[0] << " after " << elapsed_time << " has elapsed." << endl;

    }
    if (y_copy[0] < 0){
      activation_strength = activation_strength + stepsize_activation;
      //cout << "For activation_strength = " << activation_strength << ". Failure after " << elapsed_time << " seconds, y_copy[0] = " << y_copy[0] << endl;
    }
    else{
      activation_strength = activation_strength - stepsize_activation;
      //cout << "For activation_strength = " << activation_strength << ". Success after " << elapsed_time << " seconds, y_copy[0] = " << y_copy[0] << endl;
    }
    stepsize_activation = stepsize_activation / 2;
    //cout << "activation_strength of " << activation_strength << " gives y_copy[0] = " << y_copy[0] << endl;
  }



  Q_thr = activation_strength *1000*50e-9*t_A;
}
 


int main(int argc, char* argv[]) {
	ofstream SF_file;
    SF_file.open("SF_file.txt", std::ios_base::app);
    int array_length = 39;
    int counter = 0;
    char file_loc[500];
    strcat(file_loc,"Q_tot_data_"); 
    strcat(file_loc,argv[1]);
    strcat(file_loc, "_");
    strcat(file_loc,argv[2]);
    strcat(file_loc,".txt");
    cout << file_loc << endl;
	
	ifstream line_count(file_loc);
    if (!line_count.is_open()){
        cout << "No SF available" << endl;
        SF_file << argv[1] << ", " << argv[2] << ", 0" <<  endl;
        throw runtime_error("Could not open mask file");
    }

    

    

    string line;
    //Count the length of the file
    while (getline(line_count, line)){
        counter++;
    }
    int Q_tot_length = counter-array_length-1;
    PDEFIELD_TYPE Q_tot[Q_tot_length];

    PDEFIELD_TYPE val;
    PDEFIELD_TYPE ddt;
    PDEFIELD_TYPE PDEvars[array_length];

   

    
    //Read in the file
    counter = 0;
    ifstream Q_tot_file(file_loc);
    if (!Q_tot_file.is_open())
        throw runtime_error("Could not open mask file");
    while (getline(Q_tot_file, line)){
        stringstream ss(line);
        while (ss >> val){
            if (counter == 0)
                ddt = val;
            else if (counter > 0 && counter <= array_length)
                PDEvars[counter-1] = val;
            else
                Q_tot[counter-array_length] = val;
            counter++;
            //cout << "counter = " << counter << endl;
        }
    }

    bool Q_tot_max_achieved = false;

    int Q_tot_max_loc = -1;
    double Q_tot_max_val = -1;

    for (int i = 0; i < Q_tot_length; i++){
      if (Q_tot[i] > Q_tot_max_val){
        Q_tot_max_val = Q_tot[i];
        Q_tot_max_loc = i;
      }
    }

    if (Q_tot_max_achieved != Q_tot_length-1 && Q_tot_length  > Q_tot_max_loc*5)
      Q_tot_length = Q_tot_max_loc*5;


    int loga = (int) (log((double)Q_tot_length)/log(2.0));
    Q_tot_length = pow(2,loga);

    int stepper = Q_tot_length/2;
    int current_loc = Q_tot_length/2;
    PDEFIELD_TYPE Q_thr;
    double SF_lower;
    double SF_upper;
    double SF;
    while (stepper != 0){
    //for (int i = 100; i < 5000; i += 100){
      stepper /= 2;
      //cout << "stepper = " << stepper << "\n";
      //cout << "current_loc = " << current_loc << "\n";
      
      ComputeQthr(PDEvars, current_loc*ddt, ddt, Q_thr, array_length);
      SF_lower = Q_tot[current_loc]/Q_thr;
      ComputeQthr(PDEvars, (current_loc+1)*ddt, ddt, Q_thr,array_length);
      SF_upper = Q_tot[current_loc+1]/Q_thr; 
      if(SF_lower < SF_upper)
        current_loc += stepper;
      else 
        current_loc -= stepper;
      if (stepper == 0){
        if (SF_lower < SF_upper)
          SF = SF_upper;
        else
          SF = SF_lower;
      }
      
      //ComputeQthr(PDEvars, i*ddt, ddt, Q_thr, array_length);
      //SF = Q_tot[i]/Q_thr;
      //cout << "i = " << i << " and SF = " << SF << "\n";
    }

    SF_file << argv[1] << ", " << argv[2] << ", " <<  SF << endl;


    exit(0);
}
