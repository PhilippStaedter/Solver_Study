#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_xie1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r0 + 1.0*flux_r17 - 1.0*flux_r19 + 1.0*flux_r27 - 1.0*flux_r4;
    xdot[1] = 1.0*flux_r0 - 1.0*flux_r27 - 1.0*flux_r7;
    xdot[2] = 1.0*flux_r1 - 1.0*flux_r17 - 1.0*flux_r18 + 1.0*flux_r19;
    xdot[4] = 1.0*flux_r13 - 1.0*flux_r14;
    xdot[5] = -1.0*flux_r20 + 1.0*flux_r21 - 1.0*flux_r22 + 1.0*flux_r3;
    xdot[6] = -1.0*flux_r0 - 1.0*flux_r15 + 1.0*flux_r20 - 1.0*flux_r21 + 1.0*flux_r27;
    xdot[7] = -1.0*flux_r20 + 1.0*flux_r21 + 1.0*flux_r25 - 1.0*flux_r26;
    xdot[8] = -1.0*flux_r10 + 1.0*flux_r9;
    xdot[9] = 1.0*flux_r16 - 1.0*flux_r40;
    xdot[11] = 1.0*flux_r11 - 1.0*flux_r12;
    xdot[13] = -1.0*flux_r2 + 1.0*flux_r5;
    xdot[15] = -1.0*flux_r32 + 1.0*flux_r33;
    xdot[16] = 1.0*flux_r28 - 1.0*flux_r29;
    xdot[17] = 1.0*flux_r38 - 1.0*flux_r39;
    xdot[18] = -1.0*flux_r30 + 1.0*flux_r31;
    xdot[19] = 1.0*flux_r36 - 1.0*flux_r37;
    xdot[20] = 1.0*flux_r34 - 1.0*flux_r35;
    xdot[21] = 1.0*flux_r23 - 1.0*flux_r24;
    xdot[23] = 1.0*flux_r6 - 1.0*flux_r8;
}