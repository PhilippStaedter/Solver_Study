#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_ma(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = CAR1*k1;
    w[1] = REGA*incAMP*k10;
    w[2] = ACA*k11;
    w[3] = excAMP*k12;
    w[4] = excAMP*k13;
    w[5] = CAR1*k14;
    w[6] = ACA*PKA*k2;
    w[7] = incAMP*k3;
    w[8] = PKA*k4;
    w[9] = CAR1*k5;
    w[10] = ERK2*PKA*k6;
    w[11] = k7;
    w[12] = ERK2*REGA*k8;
    w[13] = ACA*k9;
}