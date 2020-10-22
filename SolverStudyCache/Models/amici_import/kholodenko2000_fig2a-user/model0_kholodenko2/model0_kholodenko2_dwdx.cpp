#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_kholodenko2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -1.0*J6_k7*MAPK*MKK_PP/pow(J6_KK7 + MAPK, 2) + 1.0*J6_k7*MKK_PP/(J6_KK7 + MAPK);
    dwdx[1] = -1.0*J7_k8*MAPK_P*MKK_PP/pow(J7_KK8 + MAPK_P, 2) + 1.0*J7_k8*MKK_PP/(J7_KK8 + MAPK_P);
    dwdx[2] = -1.0*J9_V10*MAPK_P/pow(J9_KK10 + MAPK_P, 2) + 1.0*J9_V10/(J9_KK10 + MAPK_P);
    dwdx[3] = -1.0*J0_V1*J0_n*MKKK*pow(MAPK_PP/J0_Ki, J0_n)/(MAPK_PP*(J0_K1 + MKKK)*pow(pow(MAPK_PP/J0_Ki, J0_n) + 1, 2));
    dwdx[4] = -1.0*J8_V9*MAPK_PP/pow(J8_KK9 + MAPK_PP, 2) + 1.0*J8_V9/(J8_KK9 + MAPK_PP);
    dwdx[5] = -1.0*J2_k3*MKK*MKKK_P/pow(J2_KK3 + MKK, 2) + 1.0*J2_k3*MKKK_P/(J2_KK3 + MKK);
    dwdx[6] = -1.0*J0_V1*MKKK/(pow(J0_K1 + MKKK, 2)*(pow(MAPK_PP/J0_Ki, J0_n) + 1)) + 1.0*J0_V1/((J0_K1 + MKKK)*(pow(MAPK_PP/J0_Ki, J0_n) + 1));
    dwdx[7] = -1.0*J1_V2*MKKK_P/pow(J1_KK2 + MKKK_P, 2) + 1.0*J1_V2/(J1_KK2 + MKKK_P);
    dwdx[8] = 1.0*J2_k3*MKK/(J2_KK3 + MKK);
    dwdx[9] = 1.0*J3_k4*MKK_P/(J3_KK4 + MKK_P);
    dwdx[10] = -1.0*J3_k4*MKKK_P*MKK_P/pow(J3_KK4 + MKK_P, 2) + 1.0*J3_k4*MKKK_P/(J3_KK4 + MKK_P);
    dwdx[11] = -1.0*J5_V6*MKK_P/pow(J5_KK6 + MKK_P, 2) + 1.0*J5_V6/(J5_KK6 + MKK_P);
    dwdx[12] = -1.0*J4_V5*MKK_PP/pow(J4_KK5 + MKK_PP, 2) + 1.0*J4_V5/(J4_KK5 + MKK_PP);
    dwdx[13] = 1.0*J6_k7*MAPK/(J6_KK7 + MAPK);
    dwdx[14] = 1.0*J7_k8*MAPK_P/(J7_KK8 + MAPK_P);
}