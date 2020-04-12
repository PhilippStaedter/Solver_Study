#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_borghans2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*pow(Z, 2)*v2_Vm2/(pow(Z, 2) + pow(v2_K2, 2));
    w[1] = 1.0*pow(A, 4)*pow(Y, 2)*pow(Z, 4)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4) + pow(v3_Kz, 4)));
    w[2] = 1.0*Y*v4_Kf;
    w[3] = 1.0*Z*v5_K;
    w[4] = 1.0*beta*v6_Vp;
    w[5] = 1.0*pow(A, 2)*pow(Z, v7_n)*v7_Vd/((pow(A, 2) + pow(v7_Kp, 2))*(pow(Z, v7_n) + pow(v7_Kd, v7_n)));
    w[6] = 1.0*A*v8_epsilon;
    w[7] = 1.0*beta*vin_v1 + 1.0*vin_v0;
}