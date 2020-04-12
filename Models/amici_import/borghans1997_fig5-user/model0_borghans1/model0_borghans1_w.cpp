#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_borghans1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*pow(Z, 2)*v2_Vm2/(pow(Z, 2) + pow(v2_K2, 2));
    w[1] = 1.0*Rho*pow(Y, 2)*pow(Z, 4)*a*beta*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
    w[2] = 1.0*Y*v4_Kf;
    w[3] = 1.0*Z*v5_K;
    w[4] = 1.0*Rho*pow(Z, 4)*v6_Kd;
    w[5] = 1.0*v7_Kr*(-Rho + 1);
    w[6] = 1.0*beta*vin_v1 + 1.0*vin_v0;
}