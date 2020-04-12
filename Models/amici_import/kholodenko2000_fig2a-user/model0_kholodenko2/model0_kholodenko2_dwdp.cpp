#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kholodenko2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*J0_V1*MKKK/(pow(J0_K1 + MKKK, 2)*(pow(MAPK_PP/J0_Ki, J0_n) + 1));
            break;
        case 1:
            dwdp[0] = -1.0*J0_V1*MKKK*pow(MAPK_PP/J0_Ki, J0_n)*log(MAPK_PP/J0_Ki)/((J0_K1 + MKKK)*pow(pow(MAPK_PP/J0_Ki, J0_n) + 1, 2));
            break;
        case 2:
            dwdp[0] = 1.0*J0_V1*J0_n*MKKK*pow(MAPK_PP/J0_Ki, J0_n)/(J0_Ki*(J0_K1 + MKKK)*pow(pow(MAPK_PP/J0_Ki, J0_n) + 1, 2));
            break;
        case 3:
            dwdp[0] = 1.0*MKKK/((J0_K1 + MKKK)*(pow(MAPK_PP/J0_Ki, J0_n) + 1));
            break;
        case 4:
            dwdp[1] = -1.0*J1_V2*MKKK_P/pow(J1_KK2 + MKKK_P, 2);
            break;
        case 5:
            dwdp[1] = 1.0*MKKK_P/(J1_KK2 + MKKK_P);
            break;
        case 6:
            dwdp[2] = -1.0*J2_k3*MKK*MKKK_P/pow(J2_KK3 + MKK, 2);
            break;
        case 7:
            dwdp[2] = 1.0*MKK*MKKK_P/(J2_KK3 + MKK);
            break;
        case 8:
            dwdp[3] = -1.0*J3_k4*MKKK_P*MKK_P/pow(J3_KK4 + MKK_P, 2);
            break;
        case 9:
            dwdp[3] = 1.0*MKKK_P*MKK_P/(J3_KK4 + MKK_P);
            break;
        case 10:
            dwdp[4] = -1.0*J4_V5*MKK_PP/pow(J4_KK5 + MKK_PP, 2);
            break;
        case 11:
            dwdp[4] = 1.0*MKK_PP/(J4_KK5 + MKK_PP);
            break;
        case 12:
            dwdp[5] = -1.0*J5_V6*MKK_P/pow(J5_KK6 + MKK_P, 2);
            break;
        case 13:
            dwdp[5] = 1.0*MKK_P/(J5_KK6 + MKK_P);
            break;
        case 14:
            dwdp[6] = -1.0*J6_k7*MAPK*MKK_PP/pow(J6_KK7 + MAPK, 2);
            break;
        case 15:
            dwdp[6] = 1.0*MAPK*MKK_PP/(J6_KK7 + MAPK);
            break;
        case 16:
            dwdp[7] = -1.0*J7_k8*MAPK_P*MKK_PP/pow(J7_KK8 + MAPK_P, 2);
            break;
        case 17:
            dwdp[7] = 1.0*MAPK_P*MKK_PP/(J7_KK8 + MAPK_P);
            break;
        case 18:
            dwdp[8] = -1.0*J8_V9*MAPK_PP/pow(J8_KK9 + MAPK_PP, 2);
            break;
        case 19:
            dwdp[8] = 1.0*MAPK_PP/(J8_KK9 + MAPK_PP);
            break;
        case 20:
            dwdp[9] = -1.0*J9_V10*MAPK_P/pow(J9_KK10 + MAPK_P, 2);
            break;
        case 21:
            dwdp[9] = 1.0*MAPK_P/(J9_KK10 + MAPK_P);
            break;
    }
}