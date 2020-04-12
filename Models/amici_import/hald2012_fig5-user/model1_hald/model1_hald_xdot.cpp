#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model1_hald(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 1.0*flux_r1 + 1.0*flux_r12 - 59.0*flux_r20 - 1.0*flux_r4;
    xdot[1] = -1.0*flux_r0 - 1.0*flux_r11 + 1.0*flux_r20;
    xdot[2] = 1.0*flux_r10 + 1.0*flux_r13 - 1.0*flux_r15 + 1.0*flux_r18 + 2.0*flux_r2 - 1.0*flux_r26 + 1.0*flux_r5;
    xdot[3] = -1.0*flux_r2;
    xdot[4] = -1.0*flux_r10 - 1.0*flux_r13 + 1.0*flux_r15 - 1.0*flux_r18 - 1.0*flux_r2 + 1.0*flux_r26 - 1.0*flux_r5;
    xdot[5] = -1.0*flux_r19 - 1.0*flux_r25 + 1.0*flux_r3 - 1.0*flux_r6;
    xdot[6] = 1.0*flux_r6;
    xdot[7] = -1.0*flux_r26 + 1.0*flux_r7;
    xdot[8] = -1.0*flux_r1 - 59.0*flux_r21;
    xdot[9] = 1.0*flux_r21;
    xdot[10] = -1.0*flux_r13 + 1.0*flux_r14;
    xdot[11] = 1.0*flux_r13 - 1.0*flux_r3;
    xdot[12] = 1.0*flux_r10 - 1.0*flux_r14 - 1.0*flux_r18;
    xdot[13] = 1.0*flux_r19 + 1.0*flux_r3 - 1.0*flux_r7;
    xdot[14] = -1.0*flux_r10 + 59.0*flux_r8;
    xdot[15] = -1.0*flux_r8;
    xdot[16] = -59.0*flux_r22 + 1.0*flux_r25;
    xdot[17] = 1.0*flux_r22;
    xdot[18] = -1.0*flux_r16 - 59.0*flux_r23 - 1.0*flux_r6;
    xdot[19] = -1.0*flux_r11 + 1.0*flux_r23 - 1.0*flux_r9;
    xdot[20] = -1.0*flux_r1 + 1.0*flux_r25 - 1.0*flux_r4 - 1.0*flux_r7;
    xdot[21] = 1.0*flux_r1 - 1.0*flux_r25 + 1.0*flux_r4 + 1.0*flux_r7;
    xdot[22] = -59.0*flux_r24 + 1.0*flux_r4;
    xdot[23] = 1.0*flux_r24;
    xdot[24] = -1.0*flux_r15 + 1.0*flux_r26;
    xdot[25] = -1.0*flux_r12 + 1.0*flux_r15 - 1.0*flux_r16 - 1.0*flux_r17;
    xdot[26] = 1.0*flux_r16;
    xdot[29] = 1.0*flux_r18;
    xdot[30] = 1.0*flux_r11;
}