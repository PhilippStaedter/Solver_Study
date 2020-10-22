#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_borghans3(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -2.0*pow(X, 3)*beta*v3s_Vm3s/pow(pow(X, 2) + pow(v3s_K3s, 2), 2) + 2.0*X*beta*v3s_Vm3s/(pow(X, 2) + pow(v3s_K3s, 2));
    dwdx[1] = 1.0*v6_Kf;
    dwdx[2] = -2.0*pow(Y, 3)*pow(Z, 2)*v3i_Vm3i/(pow(pow(Y, 2) + pow(v3i_K3y, 2), 2)*(pow(Z, 2) + pow(v3i_K3z, 2))) + 2.0*Y*pow(Z, 2)*v3i_Vm3i/((pow(Y, 2) + pow(v3i_K3y, 2))*(pow(Z, 2) + pow(v3i_K3z, 2)));
    dwdx[3] = 1.0*v4_Kf;
    dwdx[4] = -2.0*pow(Z, 3)*v2i_Vm2i/pow(pow(Z, 2) + pow(v2i_K2i, 2), 2) + 2.0*Z*v2i_Vm2i/(pow(Z, 2) + pow(v2i_K2i, 2));
    dwdx[5] = -2.0*pow(Z, 3)*v2s_Vm2s/pow(pow(Z, 2) + pow(v2s_K2s, 2), 2) + 2.0*Z*v2s_Vm2s/(pow(Z, 2) + pow(v2s_K2s, 2));
    dwdx[6] = -2.0*pow(Y, 2)*pow(Z, 3)*v3i_Vm3i/((pow(Y, 2) + pow(v3i_K3y, 2))*pow(pow(Z, 2) + pow(v3i_K3z, 2), 2)) + 2.0*pow(Y, 2)*Z*v3i_Vm3i/((pow(Y, 2) + pow(v3i_K3y, 2))*(pow(Z, 2) + pow(v3i_K3z, 2)));
    dwdx[7] = 1.0*v5_K;
}