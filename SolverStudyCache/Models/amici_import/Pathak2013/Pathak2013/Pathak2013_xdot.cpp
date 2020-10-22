#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Pathak2013(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r0;
    xdot[1] = -1.0*flux_r1 - 1.0*flux_r2 - 1.0*flux_r3;
    xdot[2] = -1.0*flux_r7 - 1.0*flux_r8;
    xdot[3] = -1.0*flux_r9;
    xdot[4] = -1.0*flux_r10;
    xdot[5] = -1.0*flux_r4 - 1.0*flux_r5 - 1.0*flux_r6;
    xdot[6] = 1.0*flux_r0 + 1.0*flux_r1 + 1.0*flux_r10 - 1.0*flux_r12 - 1.0*flux_r17 + 1.0*flux_r8 + 1.0*flux_r9;
    xdot[7] = -1.0*flux_r13 + 1.0*flux_r2;
    xdot[8] = -1.0*flux_r14 + 1.0*flux_r3 + 1.0*flux_r7;
    xdot[9] = -1.0*flux_r15 + 1.0*flux_r6;
    xdot[10] = -1.0*flux_r19 + 1.0*flux_r5;
    xdot[11] = -1.0*flux_r20 + 1.0*flux_r4;
    xdot[12] = -1.0*flux_r11 + 1.0*flux_r12 + 1.0*flux_r13 + 1.0*flux_r14 + 1.0*flux_r15;
    xdot[13] = 1.0*flux_r11 - 1.0*flux_r16 - 1.0*flux_r18 - 1.0*flux_r22;
    xdot[14] = 1.0*flux_r16 + 1.0*flux_r17 - 1.0*flux_r33 - 1.0*flux_r34;
    xdot[15] = 1.0*flux_r18 + 1.0*flux_r19 + 1.0*flux_r20 - 1.0*flux_r35;
    xdot[16] = -1.0*flux_r21 + 1.0*flux_r22;
    xdot[17] = 1.0*flux_r21 - 1.0*flux_r23 - 1.0*flux_r24 - 1.0*flux_r25 - 1.0*flux_r26 - 1.0*flux_r27 - 1.0*flux_r28 - 1.0*flux_r29 - 1.0*flux_r30 - 1.0*flux_r32;
    xdot[18] = -1.0*flux_r31 + 1.0*flux_r32;
    xdot[19] = 1.0*flux_r31 - 1.0*flux_r36 - 1.0*flux_r37 - 1.0*flux_r38 - 1.0*flux_r39 - 1.0*flux_r65 - 1.0*flux_r66 - 1.0*flux_r67 - 1.0*flux_r68 - 1.0*flux_r69 - 1.0*flux_r70;
    xdot[20] = 1.0*flux_r23 + 1.0*flux_r33;
    xdot[21] = 1.0*flux_r24 + 1.0*flux_r34 - 1.0*flux_r40 - 1.0*flux_r41 - 1.0*flux_r42;
    xdot[22] = 1.0*flux_r25;
    xdot[23] = 1.0*flux_r26;
    xdot[24] = 1.0*flux_r27;
    xdot[25] = 1.0*flux_r28;
    xdot[26] = 1.0*flux_r29;
    xdot[27] = 1.0*flux_r30 + 1.0*flux_r35 - 1.0*flux_r43;
    xdot[28] = 1.0*flux_r36 - 1.0*flux_r54 - 1.0*flux_r55;
    xdot[29] = 1.0*flux_r37 + 1.0*flux_r40 + 1.0*flux_r43 - 1.0*flux_r56 - 1.0*flux_r57 - 1.0*flux_r58 - 1.0*flux_r78;
    xdot[30] = 1.0*flux_r38 + 1.0*flux_r41 - 1.0*flux_r59 - 1.0*flux_r60 - 1.0*flux_r61;
    xdot[31] = 1.0*flux_r39 + 1.0*flux_r42 - 1.0*flux_r62 - 1.0*flux_r63 - 1.0*flux_r64;
    xdot[32] = -1.0*flux_r44 + 1.0*flux_r55 + 1.0*flux_r59;
    xdot[33] = 1.0*flux_r44 - 1.0*flux_r84;
    xdot[34] = -1.0*flux_r45 + 1.0*flux_r56 + 1.0*flux_r64;
    xdot[35] = 1.0*flux_r45 - 1.0*flux_r83;
    xdot[36] = -1.0*flux_r46 + 1.0*flux_r54;
    xdot[37] = 1.0*flux_r46 - 1.0*flux_r82;
    xdot[38] = -1.0*flux_r47 + 1.0*flux_r61;
    xdot[39] = 1.0*flux_r47 - 1.0*flux_r71;
    xdot[40] = -1.0*flux_r48 + 1.0*flux_r57;
    xdot[41] = 1.0*flux_r48 - 1.0*flux_r80;
    xdot[42] = -1.0*flux_r49 + 1.0*flux_r78;
    xdot[43] = 1.0*flux_r49 - 1.0*flux_r81;
    xdot[44] = -1.0*flux_r50 + 1.0*flux_r60 + 1.0*flux_r63;
    xdot[45] = 1.0*flux_r50 - 1.0*flux_r85;
    xdot[46] = -1.0*flux_r51 + 1.0*flux_r58 + 1.0*flux_r62;
    xdot[47] = 1.0*flux_r51 - 1.0*flux_r77;
    xdot[48] = -1.0*flux_r52 + 1.0*flux_r66;
    xdot[49] = 1.0*flux_r52 - 1.0*flux_r75;
    xdot[50] = -1.0*flux_r53 + 1.0*flux_r67;
    xdot[51] = 1.0*flux_r53 - 1.0*flux_r74;
    xdot[52] = 1.0*flux_r68 - 1.0*flux_r72;
    xdot[53] = 1.0*flux_r69 - 1.0*flux_r73;
    xdot[54] = 1.0*flux_r70 - 1.0*flux_r79;
    xdot[55] = 1.0*flux_r65 - 1.0*flux_r76;
    xdot[56] = 1.0*flux_r71 + 1.0*flux_r72 + 1.0*flux_r73 + 1.0*flux_r74 + 1.0*flux_r75 + 1.0*flux_r76 + 1.0*flux_r77 + 1.0*flux_r79 + 1.0*flux_r80 + 1.0*flux_r81 + 1.0*flux_r82 + 1.0*flux_r83 + 1.0*flux_r84 + 1.0*flux_r85;
}