#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_bachmann(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 2.5*flux_r16 - 2.5*flux_r17 + 2.5*flux_r18;
    xdot[1] = 2.5*flux_r14 - 2.5*flux_r15;
    xdot[2] = 3.6363636363636362*flux_r8 - 3.6363636363636362*flux_r9;
    xdot[3] = -3.6363636363636362*flux_r10 + 3.6363636363636362*flux_r9;
    xdot[4] = 3.6363636363636362*flux_r10 - 3.6363636363636362*flux_r12;
    xdot[5] = 3.6363636363636362*flux_r12 - 3.6363636363636362*flux_r13;
    xdot[6] = 3.6363636363636362*flux_r13 - 3.6363636363636362*flux_r14;
    xdot[8] = -2.5*flux_r0 + 2.5*flux_r11 + 2.5*flux_r33 + 2.5*flux_r34 + 2.5*flux_r35;
    xdot[9] = -2.5*flux_r1;
    xdot[10] = 2.5*flux_r0 - 2.5*flux_r11 - 2.5*flux_r22 - 2.5*flux_r30;
    xdot[11] = -2.5*flux_r2 + 2.5*flux_r3;
    xdot[12] = 2.5*flux_r2 - 2.5*flux_r3;
    xdot[13] = 2.5*flux_r27 - 2.5*flux_r28 + 2.5*flux_r29;
    xdot[14] = 2.5*flux_r25 - 2.5*flux_r26;
    xdot[15] = 3.6363636363636362*flux_r19 - 3.6363636363636362*flux_r20;
    xdot[16] = 3.6363636363636362*flux_r20 - 3.6363636363636362*flux_r21;
    xdot[17] = 3.6363636363636362*flux_r21 - 3.6363636363636362*flux_r23;
    xdot[18] = 3.6363636363636362*flux_r23 - 3.6363636363636362*flux_r24;
    xdot[19] = 3.6363636363636362*flux_r24 - 3.6363636363636362*flux_r25;
    xdot[20] = -2.5*flux_r4 - 2.5*flux_r5 + 2.5*flux_r7;
    xdot[21] = 3.6363636363636362*flux_r6 - 3.6363636363636362*flux_r7;
    xdot[22] = 2.5*flux_r31 + 2.5*flux_r32 - 2.5*flux_r35;
    xdot[23] = 2.5*flux_r22 - 2.5*flux_r31 - 2.5*flux_r33;
    xdot[24] = 2.5*flux_r30 - 2.5*flux_r32 - 2.5*flux_r34;
    xdot[25] = 2.5*flux_r4 + 2.5*flux_r5 - 2.5*flux_r6;
}