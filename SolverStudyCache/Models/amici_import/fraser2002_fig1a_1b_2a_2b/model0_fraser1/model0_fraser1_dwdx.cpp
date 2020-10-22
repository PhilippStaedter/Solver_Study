#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_fraser1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = eta1;
    dwdx[1] = sigma*x1;
    dwdx[2] = -a1*x1/pow(a1 + at, 2) + x1/(a1 + at);
    dwdx[3] = eta2;
    dwdx[4] = sigma*x2;
    dwdx[5] = -a2*x2*(amici_p - mu)/pow(a2 + at, 2) + x2*(amici_p - mu)/(a2 + at);
    dwdx[6] = sigma*x3;
    dwdx[7] = -a3*x3*(amici_p - mu)/pow(a3 + at, 2) + x3*(amici_p - mu)/(a3 + at);
    dwdx[8] = eta3;
    dwdx[9] = 2*mu*x1/gamma;
    dwdx[10] = a1*sigma;
    dwdx[11] = mu;
    dwdx[12] = a1/(a1 + at);
    dwdx[13] = 2*mu*x2/gamma;
    dwdx[14] = a2*sigma;
    dwdx[15] = mu;
    dwdx[16] = a2*(amici_p - mu)/(a2 + at);
    dwdx[17] = a3*sigma;
    dwdx[18] = mu;
    dwdx[19] = a3*(amici_p - mu)/(a3 + at);
    dwdx[20] = 2*mu*x3/gamma;
}