#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_borghans3(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*pow(Z, 2)*v2i_Vm2i/(pow(Z, 2) + pow(v2i_K2i, 2));
    w[1] = 1.0*pow(Z, 2)*v2s_Vm2s/(pow(Z, 2) + pow(v2s_K2s, 2));
    w[2] = 1.0*pow(Y, 2)*pow(Z, 2)*v3i_Vm3i/((pow(Y, 2) + pow(v3i_K3y, 2))*(pow(Z, 2) + pow(v3i_K3z, 2)));
    w[3] = 1.0*pow(X, 2)*beta*v3s_Vm3s/(pow(X, 2) + pow(v3s_K3s, 2));
    w[4] = 1.0*Y*v4_Kf;
    w[5] = 1.0*Z*v5_K;
    w[6] = 1.0*X*v6_Kf;
    w[7] = 1.0*beta*vin_v1 + 1.0*vin_v0;
}