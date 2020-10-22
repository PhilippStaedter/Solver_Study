#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_kholodenko2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*J0_V1*MKKK/((J0_K1 + MKKK)*(pow(MAPK_PP/J0_Ki, J0_n) + 1));
    w[1] = 1.0*J1_V2*MKKK_P/(J1_KK2 + MKKK_P);
    w[2] = 1.0*J2_k3*MKK*MKKK_P/(J2_KK3 + MKK);
    w[3] = 1.0*J3_k4*MKKK_P*MKK_P/(J3_KK4 + MKK_P);
    w[4] = 1.0*J4_V5*MKK_PP/(J4_KK5 + MKK_PP);
    w[5] = 1.0*J5_V6*MKK_P/(J5_KK6 + MKK_P);
    w[6] = 1.0*J6_k7*MAPK*MKK_PP/(J6_KK7 + MAPK);
    w[7] = 1.0*J7_k8*MAPK_P*MKK_PP/(J7_KK8 + MAPK_P);
    w[8] = 1.0*J8_V9*MAPK_PP/(J8_KK9 + MAPK_PP);
    w[9] = 1.0*J9_V10*MAPK_P/(J9_KK10 + MAPK_P);
}