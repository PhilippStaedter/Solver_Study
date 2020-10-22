#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Leber2015(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 1.0*flux_r10 - 1.0*flux_r21 - 1.0*flux_r7 - 1.0*flux_r8;
    xdot[1] = -1.0*flux_r19;
    xdot[2] = 1.0*flux_r19 - 1.0*flux_r28;
    xdot[3] = 1.0*flux_r21 - 1.0*flux_r22;
    xdot[4] = 1.0*flux_r22 - 1.0*flux_r23 - 1.0*flux_r26;
    xdot[5] = -1.0*flux_r27;
    xdot[6] = -1.0*flux_r4 + 1.0*flux_r9;
    xdot[7] = -0.25*flux_r15 + 0.25*flux_r20 - 0.25*flux_r5;
    xdot[8] = 0.25*flux_r16 - 0.25*flux_r20 + 0.25*flux_r29 + 0.25*flux_r5;
    xdot[10] = 0.25*flux_r15 - 0.25*flux_r16 - 0.25*flux_r29;
    xdot[11] = 0.25*flux_r17 - 0.25*flux_r18;
    xdot[12] = -14.285714285714285*flux_r6 + 14.285714285714285*flux_r7;
    xdot[15] = -14.285714285714285*flux_r13 + 14.285714285714285*flux_r14 - 14.285714285714285*flux_r2;
    xdot[16] = 14.285714285714285*flux_r12 - 14.285714285714285*flux_r3;
    xdot[17] = -14.285714285714285*flux_r0 + 14.285714285714285*flux_r11 + 14.285714285714285*flux_r13;
    xdot[18] = -1.0*flux_r1 - 1.0*flux_r24 - 1.0*flux_r25 + 1.0*flux_r6;
    xdot[19] = -1.0*flux_r11 + 1.0*flux_r26;
    xdot[21] = -1.0*flux_r14 + 1.0*flux_r24;
    xdot[22] = -1.0*flux_r12 + 1.0*flux_r25;
}