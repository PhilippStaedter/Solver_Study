#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_zhao1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = d;
    dwdx[1] = Sp*(-I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) + (-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    dwdx[2] = S*(-I*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) + (-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    dwdx[3] = mu;
    dwdx[4] = Sp*(-I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) + (-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    dwdx[5] = S*(-I*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) + (-alphai + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    dwdx[6] = mu;
    dwdx[7] = d;
    dwdx[8] = Sp*(-I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2));
    dwdx[9] = I*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + S*(-I*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2));
    dwdx[10] = mu;
    dwdx[11] = I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + Sp*(-I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2));
    dwdx[12] = mu;
    dwdx[13] = S*(-I*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2) - Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/pow(I + Ip + S + Sp, 2));
}