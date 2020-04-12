#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Proctor2013a(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = 1.0*flux_r47 + 1.0*flux_r48 - 1.0*flux_r50;
    xdot[1] = 1.0*flux_r64 - 1.0*flux_r65 - 1.0*flux_r66 + 1.0*flux_r67 + 1.0*flux_r71;
    xdot[2] = 1.0*flux_r62 - 1.0*flux_r63 + 1.0*flux_r88;
    xdot[3] = -1.0*flux_r12 + 1.0*flux_r13 + 1.0*flux_r20 - 1.0*flux_r21;
    xdot[4] = 1.0*flux_r16 + 1.0*flux_r17 + 1.0*flux_r18 - 1.0*flux_r19;
    xdot[5] = 1.0*flux_r54 + 1.0*flux_r55 - 1.0*flux_r57;
    xdot[6] = -1.0*flux_r5 + 1.0*flux_r6 + 1.0*flux_r8;
    xdot[7] = -1.0*flux_r59 + 1.0*flux_r61 + 1.0*flux_r7 - 1.0*flux_r8;
    xdot[8] = 1.0*flux_r59 - 1.0*flux_r61;
    xdot[9] = -1.0*flux_r78 + 1.0*flux_r79 + 1.0*flux_r80;
    xdot[10] = 1.0*flux_r78 - 1.0*flux_r79 - 1.0*flux_r80;
    xdot[11] = 1.0*flux_r10 + 1.0*flux_r11 - 1.0*flux_r9;
    xdot[12] = -1.0*flux_r10 - 1.0*flux_r11 + 1.0*flux_r9;
    xdot[13] = 1.0*flux_r109 - 1.0*flux_r110;
    xdot[14] = 1.0*flux_r68 + 1.0*flux_r69 - 1.0*flux_r70;
    xdot[15] = 1.0*flux_r25 + 1.0*flux_r26 - 1.0*flux_r28;
    xdot[16] = 1.0*flux_r32 + 1.0*flux_r33 - 1.0*flux_r35;
    xdot[17] = 1.0*flux_r38 + 1.0*flux_r39 - 1.0*flux_r41;
    xdot[18] = -1.0*flux_r22 + 1.0*flux_r23 + 1.0*flux_r24;
    xdot[19] = 1.0*flux_r22 - 1.0*flux_r23 - 1.0*flux_r24;
    xdot[20] = 1.0*flux_r52 + 1.0*flux_r53 - 1.0*flux_r56 - 1.0*flux_r58 - 1.0*flux_r59 + 1.0*flux_r60 + 1.0*flux_r61;
    xdot[21] = 1.0*flux_r27 - 1.0*flux_r29 - 1.0*flux_r30;
    xdot[22] = 1.0*flux_r34 - 1.0*flux_r36;
    xdot[23] = 1.0*flux_r40 - 1.0*flux_r42;
    xdot[24] = 1.0*flux_r89 - 1.0*flux_r90;
    xdot[25] = 1.0*flux_r93 - 1.0*flux_r94 - 1.0*flux_r95 + 1.0*flux_r96;
    xdot[26] = 1.0*flux_r91 - 1.0*flux_r92;
    xdot[27] = -1.0*flux_r81 + 1.0*flux_r82 + 1.0*flux_r83 + 1.0*flux_r87;
    xdot[28] = 1.0*flux_r81 - 1.0*flux_r82 - 1.0*flux_r83 - 1.0*flux_r86;
    xdot[29] = 1.0*flux_r115 + 1.0*flux_r116 + 1.0*flux_r117 - 1.0*flux_r45;
    xdot[30] = 1.0*flux_r118 + 1.0*flux_r119 + 1.0*flux_r120 - 1.0*flux_r122;
    xdot[31] = -1.0*flux_r58 + 1.0*flux_r60 - 1.0*flux_r7 + 1.0*flux_r8;
    xdot[32] = 1.0*flux_r58 - 1.0*flux_r60;
    xdot[33] = -1.0*flux_r104 + 1.0*flux_r105 - 1.0*flux_r124 + 1.0*flux_r125 + 1.0*flux_r49 - 1.0*flux_r51;
    xdot[34] = 1.0*flux_r104 - 1.0*flux_r105;
    xdot[35] = 1.0*flux_r124 - 1.0*flux_r125;
    xdot[37] = -1.0*flux_r108;
    xdot[38] = 1.0*flux_r108;
    xdot[39] = 1.0*flux_r106 + 1.0*flux_r107;
    xdot[40] = -1.0*flux_r106 - 1.0*flux_r107 + 1.0*flux_r108;
    xdot[41] = -1.0*flux_r0 + 1.0*flux_r1 - 1.0*flux_r2 + 1.0*flux_r3 - 1.0*flux_r4;
    xdot[42] = -1.0*flux_r126 + 1.0*flux_r127 + 1.0*flux_r29 + 1.0*flux_r30 - 1.0*flux_r31 - 1.0*flux_r98 + 1.0*flux_r99;
    xdot[43] = 1.0*flux_r98 - 1.0*flux_r99;
    xdot[44] = 1.0*flux_r126 - 1.0*flux_r127;
    xdot[45] = -1.0*flux_r100 + 1.0*flux_r101 - 1.0*flux_r128 + 1.0*flux_r129 + 1.0*flux_r36 - 1.0*flux_r37;
    xdot[46] = 1.0*flux_r100 - 1.0*flux_r101;
    xdot[47] = 1.0*flux_r128 - 1.0*flux_r129;
    xdot[48] = -1.0*flux_r102 + 1.0*flux_r103 - 1.0*flux_r130 + 1.0*flux_r131 + 1.0*flux_r42 - 1.0*flux_r43;
    xdot[49] = 1.0*flux_r102 - 1.0*flux_r103;
    xdot[50] = 1.0*flux_r130 - 1.0*flux_r131;
    xdot[51] = -1.0*flux_r74 + 1.0*flux_r75 - 1.0*flux_r76 + 1.0*flux_r77 - 1.0*flux_r97;
    xdot[52] = -1.0*flux_r100 + 1.0*flux_r101 - 1.0*flux_r102 + 1.0*flux_r103 - 1.0*flux_r104 + 1.0*flux_r105 + 1.0*flux_r44 - 1.0*flux_r46 - 1.0*flux_r98 + 1.0*flux_r99;
    xdot[53] = 1.0*flux_r121 - 1.0*flux_r123 - 1.0*flux_r124 + 1.0*flux_r125 - 1.0*flux_r126 + 1.0*flux_r127 - 1.0*flux_r128 + 1.0*flux_r129 - 1.0*flux_r130 + 1.0*flux_r131;
    xdot[54] = 1.0*flux_r0 - 1.0*flux_r1 - 1.0*flux_r5 + 1.0*flux_r6 + 1.0*flux_r7;
    xdot[55] = 1.0*flux_r2 - 1.0*flux_r3;
    xdot[56] = 1.0*flux_r5 - 1.0*flux_r6 - 1.0*flux_r7;
    xdot[57] = -1.0*flux_r0 + 1.0*flux_r1;
    xdot[58] = -1.0*flux_r2 + 1.0*flux_r3;
    xdot[59] = 1.0*flux_r74 - 1.0*flux_r75;
    xdot[60] = 1.0*flux_r76 - 1.0*flux_r77;
    xdot[61] = 1.0*flux_r95 - 1.0*flux_r96;
    xdot[62] = -1.0*flux_r74 + 1.0*flux_r75 - 1.0*flux_r95 + 1.0*flux_r96;
    xdot[63] = -1.0*flux_r76 + 1.0*flux_r77;
    xdot[64] = 1.0*flux_r72 - 1.0*flux_r73;
    xdot[65] = 1.0*flux_r66 - 1.0*flux_r67 - 1.0*flux_r71 - 1.0*flux_r72 + 1.0*flux_r73;
    xdot[66] = 1.0*flux_r12 - 1.0*flux_r13 - 2.0*flux_r14 + 2.0*flux_r15 - 1.0*flux_r72 + 1.0*flux_r73;
    xdot[67] = 1.0*flux_r14 - 1.0*flux_r15;
    xdot[68] = 1.0*flux_r111 - 1.0*flux_r112 - 1.0*flux_r113 + 1.0*flux_r114;
    xdot[69] = 1.0*flux_r113 - 1.0*flux_r114;
    xdot[70] = 1.0*flux_r84 + 1.0*flux_r85 - 1.0*flux_r87;
    xdot[71] = -1.0*flux_r84 - 1.0*flux_r85 + 1.0*flux_r86;
    xdot[72] = -1.0*flux_r113 + 1.0*flux_r114;
}