#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_borghans1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*pow(Y, 2)*pow(Z, 4)*a*beta*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
    dwdx[1] = 1.0*pow(Z, 4)*v6_Kd;
    dwdx[2] = -1.0*v7_Kr;
    dwdx[3] = -2.0*Rho*pow(Y, 3)*pow(Z, 4)*a*beta*v3_Vm3/(d*pow(pow(Y, 2) + pow(v3_Ky, 2), 2)*(pow(Z, 4)*a/d + 1)) + 2.0*Rho*Y*pow(Z, 4)*a*beta*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
    dwdx[4] = 1.0*v4_Kf;
    dwdx[5] = -2.0*pow(Z, 3)*v2_Vm2/pow(pow(Z, 2) + pow(v2_K2, 2), 2) + 2.0*Z*v2_Vm2/(pow(Z, 2) + pow(v2_K2, 2));
    dwdx[6] = -4.0*Rho*pow(Y, 2)*pow(Z, 7)*pow(a, 2)*beta*v3_Vm3/(pow(d, 2)*(pow(Y, 2) + pow(v3_Ky, 2))*pow(pow(Z, 4)*a/d + 1, 2)) + 4.0*Rho*pow(Y, 2)*pow(Z, 3)*a*beta*v3_Vm3/(d*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4)*a/d + 1));
    dwdx[7] = 1.0*v5_K;
    dwdx[8] = 4.0*Rho*pow(Z, 3)*v6_Kd;
}