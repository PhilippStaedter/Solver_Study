#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_model0_ihekwaba1(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r1 - 1.0*flux_r10 + 1.0*flux_r11 + 1.0*flux_r12 - 1.0*flux_r13 + 1.0*flux_r14 - 1.0*flux_r15 + 1.0*flux_r16 + 1.0*flux_r2 + 1.0*flux_r3 - 1.0*flux_r33 + 1.0*flux_r44 + 1.0*flux_r5 + 1.0*flux_r55 - 1.0*flux_r60 + 1.0*flux_r62 - 1.0*flux_r7 + 1.0*flux_r8 + 1.0*flux_r9;
    xdot[1] = -1.0*flux_r0 + 1.0*flux_r22 + 1.0*flux_r33 - 1.0*flux_r44 - 1.0*flux_r55;
    xdot[2] = 1.0*flux_r0 + 1.0*flux_r10 - 1.0*flux_r11 - 1.0*flux_r12 - 1.0*flux_r22;
    xdot[3] = 1.0*flux_r1 - 1.0*flux_r2 - 1.0*flux_r3 - 1.0*flux_r61 + 1.0*flux_r63;
    xdot[4] = 1.0*flux_r13 - 1.0*flux_r14 + 1.0*flux_r61 - 1.0*flux_r62 - 1.0*flux_r63;
    xdot[5] = -1.0*flux_r4 + 1.0*flux_r6 + 1.0*flux_r7 - 1.0*flux_r8 - 1.0*flux_r9;
    xdot[6] = 1.0*flux_r15 - 1.0*flux_r16 + 1.0*flux_r4 - 1.0*flux_r5 - 1.0*flux_r6;
    xdot[7] = -1.0*flux_r17 + 1.0*flux_r19 + 1.0*flux_r26 - 1.0*flux_r27 + 1.0*flux_r28 - 1.0*flux_r33 + 1.0*flux_r55 - 1.0*flux_r57;
    xdot[8] = -1.0*flux_r10 + 1.0*flux_r12 + 1.0*flux_r17 - 1.0*flux_r19 + 1.0*flux_r23 - 1.0*flux_r24;
    xdot[9] = -1.0*flux_r20 + 1.0*flux_r21 - 1.0*flux_r26 + 1.0*flux_r27;
    xdot[10] = 1.0*flux_r20 - 1.0*flux_r21 - 1.0*flux_r23;
    xdot[11] = 1.0*flux_r18 - 1.0*flux_r25 + 1.0*flux_r53;
    xdot[12] = -1.0*flux_r1 - 1.0*flux_r29 + 1.0*flux_r3 + 1.0*flux_r30 + 1.0*flux_r37 - 1.0*flux_r38 + 1.0*flux_r39 - 1.0*flux_r58;
    xdot[13] = -1.0*flux_r13 + 1.0*flux_r14 + 1.0*flux_r29 - 1.0*flux_r30 + 1.0*flux_r34 - 1.0*flux_r35;
    xdot[14] = -1.0*flux_r31 + 1.0*flux_r32 - 1.0*flux_r37 + 1.0*flux_r38;
    xdot[15] = 1.0*flux_r31 - 1.0*flux_r32 - 1.0*flux_r34;
    xdot[16] = -1.0*flux_r36 + 1.0*flux_r54;
    xdot[17] = -1.0*flux_r40 + 1.0*flux_r41 + 1.0*flux_r48 - 1.0*flux_r49 + 1.0*flux_r50 - 1.0*flux_r59 - 1.0*flux_r7 + 1.0*flux_r9;
    xdot[18] = -1.0*flux_r15 + 1.0*flux_r16 + 1.0*flux_r40 - 1.0*flux_r41 + 1.0*flux_r45 - 1.0*flux_r46;
    xdot[19] = -1.0*flux_r42 + 1.0*flux_r43 - 1.0*flux_r48 + 1.0*flux_r49;
    xdot[20] = 1.0*flux_r42 - 1.0*flux_r43 - 1.0*flux_r45;
    xdot[21] = -1.0*flux_r47 + 1.0*flux_r56;
    xdot[22] = -1.0*flux_r0 + 1.0*flux_r11 - 1.0*flux_r17 + 1.0*flux_r19 + 1.0*flux_r22 + 1.0*flux_r24 - 1.0*flux_r29 + 1.0*flux_r30 + 1.0*flux_r35 - 1.0*flux_r4 - 1.0*flux_r40 + 1.0*flux_r41 + 1.0*flux_r46 + 1.0*flux_r5 + 1.0*flux_r51 - 1.0*flux_r52 + 1.0*flux_r6 - 1.0*flux_r61 + 1.0*flux_r62 + 1.0*flux_r63;
    xdot[23] = -1.0*flux_r20 + 1.0*flux_r21 - 1.0*flux_r31 + 1.0*flux_r32 - 1.0*flux_r42 + 1.0*flux_r43 - 1.0*flux_r51 + 1.0*flux_r52;
}