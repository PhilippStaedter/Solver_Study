#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_borghans2(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -4.0*pow(A, 7)*pow(Y, 2)*pow(Z, 4)*v3_Vm3/(pow(pow(A, 4) + pow(v3_Ka, 4), 2)*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4) + pow(v3_Kz, 4))) + 4.0*pow(A, 3)*pow(Y, 2)*pow(Z, 4)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4) + pow(v3_Kz, 4)));
    dwdx[1] = -2.0*pow(A, 3)*pow(Z, v7_n)*v7_Vd/(pow(pow(A, 2) + pow(v7_Kp, 2), 2)*(pow(Z, v7_n) + pow(v7_Kd, v7_n))) + 2.0*A*pow(Z, v7_n)*v7_Vd/((pow(A, 2) + pow(v7_Kp, 2))*(pow(Z, v7_n) + pow(v7_Kd, v7_n)));
    dwdx[2] = 1.0*v8_epsilon;
    dwdx[3] = -2.0*pow(A, 4)*pow(Y, 3)*pow(Z, 4)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*pow(pow(Y, 2) + pow(v3_Ky, 2), 2)*(pow(Z, 4) + pow(v3_Kz, 4))) + 2.0*pow(A, 4)*Y*pow(Z, 4)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4) + pow(v3_Kz, 4)));
    dwdx[4] = 1.0*v4_Kf;
    dwdx[5] = -2.0*pow(Z, 3)*v2_Vm2/pow(pow(Z, 2) + pow(v2_K2, 2), 2) + 2.0*Z*v2_Vm2/(pow(Z, 2) + pow(v2_K2, 2));
    dwdx[6] = -4.0*pow(A, 4)*pow(Y, 2)*pow(Z, 7)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*(pow(Y, 2) + pow(v3_Ky, 2))*pow(pow(Z, 4) + pow(v3_Kz, 4), 2)) + 4.0*pow(A, 4)*pow(Y, 2)*pow(Z, 3)*v3_Vm3/((pow(A, 4) + pow(v3_Ka, 4))*(pow(Y, 2) + pow(v3_Ky, 2))*(pow(Z, 4) + pow(v3_Kz, 4)));
    dwdx[7] = 1.0*v5_K;
    dwdx[8] = -1.0*pow(A, 2)*pow(Z, 2*v7_n)*v7_Vd*v7_n/(Z*(pow(A, 2) + pow(v7_Kp, 2))*pow(pow(Z, v7_n) + pow(v7_Kd, v7_n), 2)) + 1.0*pow(A, 2)*pow(Z, v7_n)*v7_Vd*v7_n/(Z*(pow(A, 2) + pow(v7_Kp, 2))*(pow(Z, v7_n) + pow(v7_Kd, v7_n)));
}