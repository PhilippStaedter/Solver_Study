#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Bungay2006a(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r0;
    xdot[1] = 1.0*flux_r0 - 1.0*flux_r24;
    xdot[2] = -1.0*flux_r1 - 1.0*flux_r41;
    xdot[3] = 1.0*flux_r1 - 1.0*flux_r25 - 1.0*flux_r36 + 1.0*flux_r37;
    xdot[4] = -1.0*flux_r2;
    xdot[5] = 1.0*flux_r2 - 1.0*flux_r20 - 1.0*flux_r22 - 1.0*flux_r36;
    xdot[6] = -1.0*flux_r3;
    xdot[7] = -1.0*flux_r19 + 1.0*flux_r21 + 1.0*flux_r23 + 1.0*flux_r3 - 1.0*flux_r30 + 1.0*flux_r37;
    xdot[8] = -1.0*flux_r4;
    xdot[9] = -1.0*flux_r13 - 1.0*flux_r28 + 1.0*flux_r4;
    xdot[10] = -1.0*flux_r5;
    xdot[11] = -1.0*flux_r12 + 1.0*flux_r29 + 1.0*flux_r5;
    xdot[12] = -1.0*flux_r6;
    xdot[13] = -1.0*flux_r14 + 1.0*flux_r6;
    xdot[14] = -1.0*flux_r32 - 1.0*flux_r34 - 1.0*flux_r44 - 1.0*flux_r7;
    xdot[15] = 1.0*flux_r16 - 1.0*flux_r17 + 1.0*flux_r18 - 1.0*flux_r19 - 1.0*flux_r20 + 1.0*flux_r21 - 1.0*flux_r28 + 1.0*flux_r29 + 1.0*flux_r7;
    xdot[16] = -1.0*flux_r8;
    xdot[17] = 1.0*flux_r40 - 1.0*flux_r42 + 1.0*flux_r8;
    xdot[18] = -1.0*flux_r9;
    xdot[19] = -1.0*flux_r42 + 1.0*flux_r9;
    xdot[20] = -1.0*flux_r10;
    xdot[21] = 1.0*flux_r10 + 1.0*flux_r31;
    xdot[22] = -1.0*flux_r11;
    xdot[23] = 1.0*flux_r11 - 1.0*flux_r39;
    xdot[24] = -1.0*flux_r12 - 1.0*flux_r13;
    xdot[25] = 1.0*flux_r12 - 1.0*flux_r14 + 1.0*flux_r16 + 1.0*flux_r18 - 1.0*flux_r33;
    xdot[26] = 1.0*flux_r13 - 1.0*flux_r17;
    xdot[27] = 1.0*flux_r14 - 1.0*flux_r15;
    xdot[28] = 1.0*flux_r15 - 1.0*flux_r16;
    xdot[29] = 1.0*flux_r17 - 1.0*flux_r18;
    xdot[30] = 1.0*flux_r19 - 1.0*flux_r24 - 1.0*flux_r25 + 1.0*flux_r27;
    xdot[31] = 1.0*flux_r20 - 1.0*flux_r21;
    xdot[32] = -1.0*flux_r22 + 1.0*flux_r23 + 1.0*flux_r27 - 1.0*flux_r35 - 1.0*flux_r38 - 1.0*flux_r43;
    xdot[33] = 1.0*flux_r22 - 1.0*flux_r23;
    xdot[34] = 1.0*flux_r24 - 1.0*flux_r26;
    xdot[35] = 1.0*flux_r25 + 1.0*flux_r26 - 1.0*flux_r27;
    xdot[36] = -1.0*flux_r30 + 1.0*flux_r31 + 1.0*flux_r42;
    xdot[37] = -1.0*flux_r32;
    xdot[38] = -1.0*flux_r34 - 1.0*flux_r35 - 1.0*flux_r41;
    xdot[39] = 1.0*flux_r35;
    xdot[40] = 1.0*flux_r32 - 1.0*flux_r33;
    xdot[41] = 1.0*flux_r33;
    xdot[42] = 1.0*flux_r30 - 1.0*flux_r31;
    xdot[43] = 1.0*flux_r34;
    xdot[44] = 1.0*flux_r28 - 1.0*flux_r29;
    xdot[45] = 1.0*flux_r36 - 1.0*flux_r37;
    xdot[46] = -1.0*flux_r38;
    xdot[47] = 1.0*flux_r38 - 1.0*flux_r39 + 1.0*flux_r40;
    xdot[48] = 1.0*flux_r39 - 1.0*flux_r40;
    xdot[49] = 1.0*flux_r41;
    xdot[50] = -100.0*flux_r0 - 100.0*flux_r1 - 100.0*flux_r10 - 100.0*flux_r11 - 100.0*flux_r2 + 100.0*flux_r27 - 100.0*flux_r3 - 100.0*flux_r4 - 100.0*flux_r5 - 100.0*flux_r6 - 100.0*flux_r7 - 100.0*flux_r8 - 100.0*flux_r9;
    xdot[51] = -1.0*flux_r43 - 1.0*flux_r44;
    xdot[52] = 1.0*flux_r43;
    xdot[53] = 1.0*flux_r44;
}