#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model1_mcauley1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[1] = -1.0*flux_r24 + 1.0*flux_r34 + 1.0*flux_r5;
    xdot[2] = 1.0*flux_r16 + 1.0*flux_r17 - 1.0*flux_r20 + 1.0*flux_r21 - 1.0*flux_r23 - 1.0*flux_r24 + 1.0*flux_r25;
    xdot[4] = 1.0*flux_r3 - 1.0*flux_r4;
    xdot[8] = -1.0*flux_r10 + 1.0*flux_r26 + 1.0*flux_r6 - 1.0*flux_r9;
    xdot[9] = 1.0*flux_r7 - 1.0*flux_r8;
    xdot[11] = 1.0*flux_r0 + 1.0*flux_r1 + 1.0*flux_r11 - 1.0*flux_r32 - 1.0*flux_r33;
    xdot[12] = 1.0*flux_r8;
    xdot[13] = 1.0*flux_r10 - 1.0*flux_r12 - 1.0*flux_r13;
    xdot[15] = 1.0*flux_r13 - 1.0*flux_r14 - 1.0*flux_r15 - 1.0*flux_r16 - 1.0*flux_r17 + 1.0*flux_r27;
    xdot[17] = 1.0*flux_r18 - 1.0*flux_r19;
    xdot[19] = 1.0*flux_r19;
    xdot[20] = 1.0*flux_r20 - 1.0*flux_r21;
    xdot[21] = 1.0*flux_r23;
    xdot[23] = 1.0*flux_r24 - 1.0*flux_r26 - 1.0*flux_r27 - 1.0*flux_r28;
    xdot[28] = -1.0*flux_r22 + 1.0*flux_r29 + 1.0*flux_r31;
    xdot[29] = 1.0*flux_r22 - 1.0*flux_r29 - 1.0*flux_r30;
    xdot[30] = 1.0*flux_r30;
    xdot[31] = -1.0*flux_r1 + 1.0*flux_r12 + 1.0*flux_r14 + 1.0*flux_r15 + 1.0*flux_r2 + 1.0*flux_r28 - 1.0*flux_r3 - 1.0*flux_r31 + 1.0*flux_r32 + 1.0*flux_r4 - 1.0*flux_r6 + 1.0*flux_r9;
    xdot[32] = 1.0*flux_r33;
}