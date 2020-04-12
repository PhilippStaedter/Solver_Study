#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_leloup1_Fig2A(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*PT_complex_degradation_k_dC;
    dwdx[1] = -1.0*PT_complex_formation_k4;
    dwdx[2] = 1.0*PT_complex_nucleation_k1;
    dwdx[3] = -1.0*pow(Cn, Mp_production_n)*pow(Mp_production_K_IP, Mp_production_n)*Mp_production_n*Mp_production_v_sP/(Cn*pow(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n), 2));
    dwdx[4] = -1.0*pow(Cn, Mt_production_n)*pow(Mt_production_K_IT, Mt_production_n)*Mt_production_V_sT*Mt_production_n/(Cn*pow(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n), 2));
    dwdx[5] = -1.0*PT_complex_nucleation_k2;
    dwdx[6] = 1.0*PTnucl_complex_degradation_k_dN;
    dwdx[7] = -1.0*Mp*Mp_degradation_V_mP/pow(Mp + Mp_degradation_K_mP, 2) + 1.0*Mp_degradation_V_mP/(Mp + Mp_degradation_K_mP) + 1.0*Mp_degradation_k_d;
    dwdx[8] = 1.0*P0_production_k_sP;
    dwdx[9] = -1.0*Mt*V_mT/pow(Mt + Mt_degradation_K_mT, 2) + 1.0*Mt_degradation_k_d + 1.0*V_mT/(Mt + Mt_degradation_K_mT);
    dwdx[10] = 1.0*T0_production_k_sT;
    dwdx[11] = 1.0*P0_degradation_k_d;
    dwdx[12] = -1.0*P0*P0_to_P1_V_1P/pow(P0 + P0_to_P1_K1_P, 2) + 1.0*P0_to_P1_V_1P/(P0 + P0_to_P1_K1_P);
    dwdx[13] = 1.0*P1_degradation_k_d;
    dwdx[14] = -1.0*P1*P1_to_P0_V_2P/pow(P1 + P1_to_P0_K_2P, 2) + 1.0*P1_to_P0_V_2P/(P1 + P1_to_P0_K_2P);
    dwdx[15] = -1.0*P1*P1_to_P2_V_3P/pow(P1 + P1_to_P2_K_3P, 2) + 1.0*P1_to_P2_V_3P/(P1 + P1_to_P2_K_3P);
    dwdx[16] = -1.0*P2*P2_degradation_V_dP/pow(P2 + P2_degradation_K_dP, 2) + 1.0*P2_degradation_V_dP/(P2 + P2_degradation_K_dP) + 1.0*P2_degradation_k_d;
    dwdx[17] = -1.0*P2*P2_to_P1_V_4P/pow(P2 + P2_to_P1_K_4P, 2) + 1.0*P2_to_P1_V_4P/(P2 + P2_to_P1_K_4P);
    dwdx[18] = 1.0*PT_complex_formation_k3*T2;
    dwdx[19] = 1.0*T0_degradation_k_d;
    dwdx[20] = -1.0*T0*T0_to_T1_V_1T/pow(T0 + T0_to_T1_K_1T, 2) + 1.0*T0_to_T1_V_1T/(T0 + T0_to_T1_K_1T);
    dwdx[21] = 1.0*T1_degradation_k_d;
    dwdx[22] = -1.0*T1*T1_to_T0_V_2T/pow(T1 + T1_to_T0_K_2T, 2) + 1.0*T1_to_T0_V_2T/(T1 + T1_to_T0_K_2T);
    dwdx[23] = -1.0*T1*T1_to_T2_V_3T/pow(T1 + T1_to_T2_K_3T, 2) + 1.0*T1_to_T2_V_3T/(T1 + T1_to_T2_K_3T);
    dwdx[24] = 1.0*P2*PT_complex_formation_k3;
    dwdx[25] = -1.0*T2*V_dT/pow(T2 + T2_degradation_K_dT, 2) + 1.0*T2_degradation_k_d + 1.0*V_dT/(T2 + T2_degradation_K_dT);
    dwdx[26] = -1.0*T2*T2_to_T1_V_4T/pow(T2 + T2_to_T1_K_4T, 2) + 1.0*T2_to_T1_V_4T/(T2 + T2_to_T1_K_4T);
}