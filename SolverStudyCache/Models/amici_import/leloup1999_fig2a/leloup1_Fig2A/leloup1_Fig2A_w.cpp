#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_leloup1_Fig2A(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*Mp*Mp_degradation_V_mP/(Mp + Mp_degradation_K_mP) + 1.0*Mp*Mp_degradation_k_d;
    w[1] = 1.0*pow(Mp_production_K_IP, Mp_production_n)*Mp_production_v_sP/(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n));
    w[2] = 1.0*Mt*Mt_degradation_k_d + 1.0*Mt*V_mT/(Mt + Mt_degradation_K_mT);
    w[3] = 1.0*pow(Mt_production_K_IT, Mt_production_n)*Mt_production_V_sT/(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n));
    w[4] = 1.0*P0*P0_degradation_k_d;
    w[5] = 1.0*Mp*P0_production_k_sP;
    w[6] = 1.0*P0*P0_to_P1_V_1P/(P0 + P0_to_P1_K1_P);
    w[7] = 1.0*P1*P1_degradation_k_d;
    w[8] = 1.0*P1*P1_to_P0_V_2P/(P1 + P1_to_P0_K_2P);
    w[9] = 1.0*P1*P1_to_P2_V_3P/(P1 + P1_to_P2_K_3P);
    w[10] = 1.0*P2*P2_degradation_V_dP/(P2 + P2_degradation_K_dP) + 1.0*P2*P2_degradation_k_d;
    w[11] = 1.0*P2*P2_to_P1_V_4P/(P2 + P2_to_P1_K_4P);
    w[12] = 1.0*CC*PT_complex_degradation_k_dC;
    w[13] = -1.0*CC*PT_complex_formation_k4 + 1.0*P2*PT_complex_formation_k3*T2;
    w[14] = 1.0*CC*PT_complex_nucleation_k1 - 1.0*Cn*PT_complex_nucleation_k2;
    w[15] = 1.0*Cn*PTnucl_complex_degradation_k_dN;
    w[16] = 1.0*T0*T0_degradation_k_d;
    w[17] = 1.0*Mt*T0_production_k_sT;
    w[18] = 1.0*T0*T0_to_T1_V_1T/(T0 + T0_to_T1_K_1T);
    w[19] = 1.0*T1*T1_degradation_k_d;
    w[20] = 1.0*T1*T1_to_T0_V_2T/(T1 + T1_to_T0_K_2T);
    w[21] = 1.0*T1*T1_to_T2_V_3T/(T1 + T1_to_T2_K_3T);
    w[22] = 1.0*T2*T2_degradation_k_d + 1.0*T2*V_dT/(T2 + T2_degradation_K_dT);
    w[23] = 1.0*T2*T2_to_T1_V_4T/(T2 + T2_to_T1_K_4T);
}