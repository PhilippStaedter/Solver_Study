#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_zhao1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = Lambda*amici_k;
    w[1] = I*d;
    w[2] = Sp*(I*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + Ip*(-alphai + 1)*(-alphas + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    w[3] = Sp*mu;
    w[4] = Lambda*(-amici_k + 1);
    w[5] = S*(I*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp) + Ip*(-alphai + 1)*(-pow(-ba + 1, n) + 1)/(I + Ip + S + Sp));
    w[6] = S*mu;
    w[7] = Ip*mu;
    w[8] = Ip*d;
    w[9] = I*mu;
}