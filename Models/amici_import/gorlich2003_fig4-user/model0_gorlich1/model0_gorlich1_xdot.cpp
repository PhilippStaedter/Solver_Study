#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_gorlich1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[2] = -83333333333.333328*flux_r4 + 83333333333.333328*flux_r8;
    xdot[3] = 83333333333.333328*flux_r1 - 83333333333.333328*flux_r2;
    xdot[4] = -83333333333.333328*flux_r1 + 83333333333.333328*flux_r4;
    xdot[5] = 83333333333.333328*flux_r2 - 83333333333.333328*flux_r8;
    xdot[6] = 55555555555.555557*flux_r5 - 55555555555.555557*flux_r7;
    xdot[8] = 55555555555.555557*flux_r3 + 55555555555.555557*flux_r5 + 55555555555.555557*flux_r6;
    xdot[9] = -83333333333.333328*flux_r3 - 83333333333.333328*flux_r4;
    xdot[10] = -55555555555.555557*flux_r5 + 55555555555.555557*flux_r7;
    xdot[11] = 55555555555.555557*flux_r0 - 55555555555.555557*flux_r6 - 55555555555.555557*flux_r7;
    xdot[12] = -83333333333.333328*flux_r0 + 83333333333.333328*flux_r8;
}