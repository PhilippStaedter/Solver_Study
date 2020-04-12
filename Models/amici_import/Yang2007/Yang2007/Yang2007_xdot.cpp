#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Yang2007(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 1.0*flux_r0 - 1.0*flux_r1 - 1.0*flux_r3 - 1.0*flux_r31 - 1.0*flux_r5 - 1.0*flux_r9;
    xdot[1] = -1.0*flux_r10 - 1.0*flux_r11 + 1.0*flux_r9;
    xdot[2] = 1.0*flux_r10 - 1.0*flux_r27;
    xdot[3] = 1.0*flux_r11 - 1.0*flux_r12 - 1.0*flux_r23;
    xdot[4] = 1.0*flux_r12 - 1.0*flux_r13 - 1.0*flux_r22;
    xdot[5] = 1.0*flux_r13;
    xdot[7] = 1.0*flux_r14 - 1.0*flux_r26;
    xdot[8] = -1.0*flux_r15;
    xdot[11] = 1.0*flux_r1 - 1.0*flux_r2 - 1.0*flux_r25;
    xdot[12] = -1.0*flux_r16 - 1.0*flux_r20;
    xdot[13] = 1.0*flux_r17 - 1.0*flux_r28 - 1.0*flux_r29 - 1.0*flux_r30;
    xdot[14] = -1.0*flux_r18;
    xdot[18] = 1.0*flux_r2 - 1.0*flux_r24;
    xdot[19] = 1.0*flux_r3 - 1.0*flux_r4;
    xdot[20] = 1.0*flux_r4;
    xdot[21] = 1.0*flux_r5 - 1.0*flux_r6 - 1.0*flux_r7;
    xdot[22] = 1.0*flux_r6;
    xdot[23] = -1.0*flux_r21 + 1.0*flux_r7;
    xdot[24] = -1.0*flux_r19 + 1.0*flux_r8;
}