#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_fraser2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -pow(a1, 2)*eta1/pow(a1 + at, 2);
            dwdp[3] = -a3*x3*(amici_p - mu)/pow(a3 + at, 2);
            dwdp[9] = -a1*x1*(amici_p - mu)/pow(a1 + at, 2);
            dwdp[10] = -pow(a2, 2)*eta2/pow(a2 + at, 2);
            dwdp[13] = -a2*x2*(amici_p - mu)/pow(a2 + at, 2);
            dwdp[14] = -pow(a3, 2)*eta3/pow(a3 + at, 2);
            break;
        case 1:
            dwdp[0] = pow(a1, 2)/(a1 + at);
            break;
        case 2:
            dwdp[10] = pow(a2, 2)/(a2 + at);
            break;
        case 3:
            dwdp[14] = pow(a3, 2)/(a3 + at);
            break;
        case 4:
            dwdp[4] = -mu*pow(x1, 2)/pow(gamma, 2);
            dwdp[5] = -mu*pow(x2, 2)/pow(gamma, 2);
            dwdp[6] = -mu*pow(x3, 2)/pow(gamma, 2);
            break;
        case 5:
            dwdp[2] = x3;
            dwdp[3] = -a3*x3/(a3 + at);
            dwdp[4] = pow(x1, 2)/gamma;
            dwdp[5] = pow(x2, 2)/gamma;
            dwdp[6] = pow(x3, 2)/gamma;
            dwdp[8] = x1;
            dwdp[9] = -a1*x1/(a1 + at);
            dwdp[12] = x2;
            dwdp[13] = -a2*x2/(a2 + at);
            break;
        case 6:
            dwdp[3] = a3*x3/(a3 + at);
            dwdp[9] = a1*x1/(a1 + at);
            dwdp[13] = a2*x2/(a2 + at);
            break;
        case 7:
            dwdp[1] = a3*x3;
            dwdp[7] = a1*x1;
            dwdp[11] = a2*x2;
            break;
    }
}