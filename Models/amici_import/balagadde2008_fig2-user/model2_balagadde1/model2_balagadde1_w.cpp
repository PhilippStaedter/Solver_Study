#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model2_balagadde1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*C1*kc1*(1 - (C1 + C2)/Cm);
    w[1] = 1.0*C1*(D + K1*(pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.5)/(pow(A2, 2) + K1));
    w[2] = 1.0*C2*kc2*(1 - (C1 + C2)/Cm);
    w[3] = 1.0*C2*(pow(A1, 2)*d2/(pow(A1, 2) + K2) + D);
    w[4] = 1.0*C1*kA1;
    w[5] = 1.0*A1*(D + dAA1);
    w[6] = 1.0*C2*(0.029999999999999999*pow(IPTG, 2)/(pow(IPTG, 2) + 25) + 0.02);
    w[7] = 1.0*A2*(D + dAA2);
}