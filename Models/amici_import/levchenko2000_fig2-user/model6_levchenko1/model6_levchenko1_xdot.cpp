#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model6_levchenko1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r10 + 1.0*flux_r12 + 1.0*flux_r16;
    xdot[1] = 1.0*flux_r10 - 1.0*flux_r12 - 1.0*flux_r13;
    xdot[2] = -1.0*flux_r14 + 1.0*flux_r15 + 1.0*flux_r16 - 1.0*flux_r20 + 1.0*flux_r21 + 1.0*flux_r23;
    xdot[3] = 1.0*flux_r13 - 1.0*flux_r14 + 1.0*flux_r15 - 1.0*flux_r17 + 1.0*flux_r18 + 1.0*flux_r23;
    xdot[4] = 1.0*flux_r14 - 1.0*flux_r15 - 1.0*flux_r16;
    xdot[5] = 1.0*flux_r17 - 1.0*flux_r18 - 1.0*flux_r19;
    xdot[6] = 1.0*flux_r19 - 1.0*flux_r20 + 1.0*flux_r21;
    xdot[7] = 1.0*flux_r20 - 1.0*flux_r21 - 1.0*flux_r23;
    xdot[8] = -1.0*flux_r27 + 1.0*flux_r28 + 1.0*flux_r3;
    xdot[9] = -1.0*flux_r1 + 1.0*flux_r2 + 1.0*flux_r3 - 1.0*flux_r7 + 1.0*flux_r8 + 1.0*flux_r9;
    xdot[10] = 1.0*flux_r27 - 1.0*flux_r28 - 1.0*flux_r29;
    xdot[11] = -1.0*flux_r1 + 1.0*flux_r2 + 1.0*flux_r29 - 1.0*flux_r4 + 1.0*flux_r5 + 1.0*flux_r9;
    xdot[12] = 1.0*flux_r1 - 1.0*flux_r2 - 1.0*flux_r3;
    xdot[13] = 1.0*flux_r4 - 1.0*flux_r5 - 1.0*flux_r6;
    xdot[14] = -1.0*flux_r10 + 1.0*flux_r12 + 1.0*flux_r13 - 1.0*flux_r17 + 1.0*flux_r18 + 1.0*flux_r19 + 1.0*flux_r6 - 1.0*flux_r7 + 1.0*flux_r8;
    xdot[15] = 1.0*flux_r7 - 1.0*flux_r8 - 1.0*flux_r9;
    xdot[16] = -1.0*flux_r0 + 1.0*flux_r11 + 1.0*flux_r26;
    xdot[17] = -1.0*flux_r0 + 1.0*flux_r11 + 1.0*flux_r22;
    xdot[18] = -1.0*flux_r24 + 1.0*flux_r25 + 1.0*flux_r26;
    xdot[19] = 1.0*flux_r0 - 1.0*flux_r11 - 1.0*flux_r22;
    xdot[20] = 1.0*flux_r22 - 1.0*flux_r24 + 1.0*flux_r25 - 1.0*flux_r27 + 1.0*flux_r28 + 1.0*flux_r29 - 1.0*flux_r4 + 1.0*flux_r5 + 1.0*flux_r6;
    xdot[21] = 1.0*flux_r24 - 1.0*flux_r25 - 1.0*flux_r26;
}