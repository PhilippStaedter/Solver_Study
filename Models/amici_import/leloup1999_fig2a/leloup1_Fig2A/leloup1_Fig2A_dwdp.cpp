#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_leloup1_Fig2A(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[22] = 1.0*T2/(T2 + T2_degradation_K_dT);
            break;
        case 1:
            dwdp[2] = 1.0*Mt/(Mt + Mt_degradation_K_mT);
            break;
        case 2:
            dwdp[0] = -1.0*Mp*Mp_degradation_V_mP/pow(Mp + Mp_degradation_K_mP, 2);
            break;
        case 3:
            dwdp[0] = 1.0*Mp/(Mp + Mp_degradation_K_mP);
            break;
        case 4:
            dwdp[0] = 1.0*Mp;
            break;
        case 5:
            dwdp[1] = 1.0*pow(Mp_production_K_IP, Mp_production_n)*Mp_production_v_sP*log(Mp_production_K_IP)/(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n)) + 1.0*pow(Mp_production_K_IP, Mp_production_n)*Mp_production_v_sP*(-pow(Cn, Mp_production_n)*log(Cn) - pow(Mp_production_K_IP, Mp_production_n)*log(Mp_production_K_IP))/pow(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n), 2);
            break;
        case 6:
            dwdp[1] = -1.0*pow(Mp_production_K_IP, 2*Mp_production_n)*Mp_production_n*Mp_production_v_sP/(Mp_production_K_IP*pow(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n), 2)) + 1.0*pow(Mp_production_K_IP, Mp_production_n)*Mp_production_n*Mp_production_v_sP/(Mp_production_K_IP*(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n)));
            break;
        case 7:
            dwdp[1] = 1.0*pow(Mp_production_K_IP, Mp_production_n)/(pow(Cn, Mp_production_n) + pow(Mp_production_K_IP, Mp_production_n));
            break;
        case 8:
            dwdp[2] = -1.0*Mt*V_mT/pow(Mt + Mt_degradation_K_mT, 2);
            break;
        case 9:
            dwdp[2] = 1.0*Mt;
            break;
        case 10:
            dwdp[3] = 1.0*pow(Mt_production_K_IT, Mt_production_n)*Mt_production_V_sT*log(Mt_production_K_IT)/(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n)) + 1.0*pow(Mt_production_K_IT, Mt_production_n)*Mt_production_V_sT*(-pow(Cn, Mt_production_n)*log(Cn) - pow(Mt_production_K_IT, Mt_production_n)*log(Mt_production_K_IT))/pow(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n), 2);
            break;
        case 11:
            dwdp[3] = -1.0*pow(Mt_production_K_IT, 2*Mt_production_n)*Mt_production_V_sT*Mt_production_n/(Mt_production_K_IT*pow(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n), 2)) + 1.0*pow(Mt_production_K_IT, Mt_production_n)*Mt_production_V_sT*Mt_production_n/(Mt_production_K_IT*(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n)));
            break;
        case 12:
            dwdp[3] = 1.0*pow(Mt_production_K_IT, Mt_production_n)/(pow(Cn, Mt_production_n) + pow(Mt_production_K_IT, Mt_production_n));
            break;
        case 13:
            dwdp[4] = 1.0*P0;
            break;
        case 14:
            dwdp[5] = 1.0*Mp;
            break;
        case 15:
            dwdp[6] = 1.0*P0/(P0 + P0_to_P1_K1_P);
            break;
        case 16:
            dwdp[6] = -1.0*P0*P0_to_P1_V_1P/pow(P0 + P0_to_P1_K1_P, 2);
            break;
        case 17:
            dwdp[7] = 1.0*P1;
            break;
        case 18:
            dwdp[8] = 1.0*P1/(P1 + P1_to_P0_K_2P);
            break;
        case 19:
            dwdp[8] = -1.0*P1*P1_to_P0_V_2P/pow(P1 + P1_to_P0_K_2P, 2);
            break;
        case 20:
            dwdp[9] = 1.0*P1/(P1 + P1_to_P2_K_3P);
            break;
        case 21:
            dwdp[9] = -1.0*P1*P1_to_P2_V_3P/pow(P1 + P1_to_P2_K_3P, 2);
            break;
        case 22:
            dwdp[10] = -1.0*P2*P2_degradation_V_dP/pow(P2 + P2_degradation_K_dP, 2);
            break;
        case 23:
            dwdp[10] = 1.0*P2/(P2 + P2_degradation_K_dP);
            break;
        case 24:
            dwdp[10] = 1.0*P2;
            break;
        case 25:
            dwdp[11] = 1.0*P2/(P2 + P2_to_P1_K_4P);
            break;
        case 26:
            dwdp[11] = -1.0*P2*P2_to_P1_V_4P/pow(P2 + P2_to_P1_K_4P, 2);
            break;
        case 27:
            dwdp[12] = 1.0*CC;
            break;
        case 28:
            dwdp[13] = -1.0*CC;
            break;
        case 29:
            dwdp[13] = 1.0*P2*T2;
            break;
        case 30:
            dwdp[14] = -1.0*Cn;
            break;
        case 31:
            dwdp[14] = 1.0*CC;
            break;
        case 32:
            dwdp[15] = 1.0*Cn;
            break;
        case 33:
            dwdp[16] = 1.0*T0;
            break;
        case 34:
            dwdp[17] = 1.0*Mt;
            break;
        case 35:
            dwdp[18] = 1.0*T0/(T0 + T0_to_T1_K_1T);
            break;
        case 36:
            dwdp[18] = -1.0*T0*T0_to_T1_V_1T/pow(T0 + T0_to_T1_K_1T, 2);
            break;
        case 37:
            dwdp[19] = 1.0*T1;
            break;
        case 38:
            dwdp[20] = 1.0*T1/(T1 + T1_to_T0_K_2T);
            break;
        case 39:
            dwdp[20] = -1.0*T1*T1_to_T0_V_2T/pow(T1 + T1_to_T0_K_2T, 2);
            break;
        case 40:
            dwdp[21] = 1.0*T1/(T1 + T1_to_T2_K_3T);
            break;
        case 41:
            dwdp[21] = -1.0*T1*T1_to_T2_V_3T/pow(T1 + T1_to_T2_K_3T, 2);
            break;
        case 42:
            dwdp[22] = -1.0*T2*V_dT/pow(T2 + T2_degradation_K_dT, 2);
            break;
        case 43:
            dwdp[22] = 1.0*T2;
            break;
        case 44:
            dwdp[23] = 1.0*T2/(T2 + T2_to_T1_K_4T);
            break;
        case 45:
            dwdp[23] = -1.0*T2*T2_to_T1_V_4T/pow(T2 + T2_to_T1_K_4T, 2);
            break;
    }
}