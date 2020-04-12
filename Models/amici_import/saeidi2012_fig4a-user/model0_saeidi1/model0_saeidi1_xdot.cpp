#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_saeidi1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r0;
    xdot[1] = -1.0*flux_r5;
    xdot[2] = -1.0*flux_r1 + 1.0*flux_r5;
    xdot[3] = 1.0*flux_r2 - 1.0*flux_r3 - 1.0*flux_r4;
    xdot[4] = 1.0*flux_r0 - 1.0*flux_r2;
    xdot[5] = 1.0*flux_r3;
    xdot[6] = -1.0*flux_r4 - 1.0*flux_r6;
    xdot[7] = 1.0*flux_r4 - 1.0*flux_r5;
    xdot[8] = 1.0*flux_r1;
    xdot[9] = 1.0*flux_r6;
}