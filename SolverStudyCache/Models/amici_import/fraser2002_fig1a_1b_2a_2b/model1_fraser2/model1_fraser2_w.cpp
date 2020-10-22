#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model1_fraser2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = pow(a1, 2)*eta1/(a1 + at);
    w[1] = a3*sigma*x3;
    w[2] = mu*x3;
    w[3] = a3*x3*(amici_p - mu)/(a3 + at);
    w[4] = mu*pow(x1, 2)/gamma;
    w[5] = mu*pow(x2, 2)/gamma;
    w[6] = mu*pow(x3, 2)/gamma;
    w[7] = a1*sigma*x1;
    w[8] = mu*x1;
    w[9] = a1*x1*(amici_p - mu)/(a1 + at);
    w[10] = pow(a2, 2)*eta2/(a2 + at);
    w[11] = a2*sigma*x2;
    w[12] = mu*x2;
    w[13] = a2*x2*(amici_p - mu)/(a2 + at);
    w[14] = pow(a3, 2)*eta3/(a3 + at);
}