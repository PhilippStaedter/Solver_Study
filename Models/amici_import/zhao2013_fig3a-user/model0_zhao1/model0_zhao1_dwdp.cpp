#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_zhao1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = amici_k;
            dwdp[4] = -amici_k + 1;
            break;
        case 1:
            dwdp[2] = -Ip*Sp*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp);
            dwdp[5] = -Ip*S*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp);
            break;
        case 2:
            dwdp[2] = Sp*(-I*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) - Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
            break;
        case 3:
            dwdp[2] = Sp*(I*n*(-alphas + 1)*pow(-ba + 1, n)/((-ba + 1)*(I + Ip + S + Sp)) + Ip*n*(-alphai + 1)*(-alphas + 1)*pow(-ba + 1, n)/((-ba + 1)*(I + Ip + S + Sp)));
            dwdp[5] = S*(I*n*pow(-ba + 1, n)/((-ba + 1)*(I + Ip + S + Sp)) + Ip*n*(-alphai + 1)*pow(-ba + 1, n)/((-ba + 1)*(I + Ip + S + Sp)));
            break;
        case 4:
            dwdp[1] = I;
            dwdp[8] = Ip;
            break;
        case 5:
            dwdp[0] = Lambda;
            dwdp[4] = -Lambda;
            break;
        case 6:
            dwdp[3] = Sp;
            dwdp[6] = S;
            dwdp[7] = Ip;
            dwdp[9] = I;
            break;
        case 7:
            dwdp[2] = Sp*(-I*(-alphas + 1)*pow(-ba + 1, n)*log(-ba + 1)/(I + Ip + S + Sp) - Ip*(-alphai + 1)*(-alphas + 1)*pow(-ba + 1, n)*log(-ba + 1)/(I + Ip + S + Sp));
            dwdp[5] = S*(-I*pow(-ba + 1, n)*log(-ba + 1)/(I + Ip + S + Sp) - Ip*(-alphai + 1)*pow(-ba + 1, n)*log(-ba + 1)/(I + Ip + S + Sp));
            break;
    }
}