#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void xdot_Ouzounoglou2014(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot[0] = -1.0*flux_r17 - 1.0*flux_r24 - 1.0*flux_r29 + 1.0*flux_r36 + 1.0*flux_r37 - 1.0*flux_r4 + 1.0*flux_r69 - 1.0*flux_r7 + 1.0*flux_r74 + 1.0*flux_r75 + 1.0*flux_r76 + 1.0*flux_r77 + 1.0*flux_r78 + 1.0*flux_r79 + 1.0*flux_r80 - 1.0*flux_r81 - 1.0*flux_r82 - 1.0*flux_r83 - 1.0*flux_r84;
    xdot[1] = 1.0*flux_r36 - 1.0*flux_r68;
    xdot[2] = 1.0*flux_r37 - 1.0*flux_r67;
    xdot[3] = -1.0*flux_r36 + 1.0*flux_r4 - 1.0*flux_r92;
    xdot[4] = -1.0*flux_r37 + 1.0*flux_r7 + 1.0*flux_r92 - 1.0*flux_r93;
    xdot[5] = 1.0*flux_r67;
    xdot[6] = 1.0*flux_r68;
    xdot[7] = 1.0*flux_r32 - 1.0*flux_r70;
    xdot[8] = 1.0*flux_r70 - 1.0*flux_r71;
    xdot[9] = 1.0*flux_r71 - 1.0*flux_r72;
    xdot[10] = 1.0*flux_r33 - 1.0*flux_r34;
    xdot[11] = -1.0*flux_r33 + 1.0*flux_r73;
    xdot[12] = 1.0*flux_r72 - 1.0*flux_r73;
    xdot[13] = 1.0*flux_r34 - 1.0*flux_r35;
    xdot[14] = 1.0*flux_r35;
    xdot[15] = -1.0*flux_r74 + 1.0*flux_r81 + 1.0*flux_r93 - 1.0*flux_r94;
    xdot[16] = -1.0*flux_r75 + 1.0*flux_r82 + 1.0*flux_r94 - 1.0*flux_r95;
    xdot[17] = 1.0*flux_r17 - 1.0*flux_r76 + 1.0*flux_r95 - 1.0*flux_r96;
    xdot[18] = -1.0*flux_r77 + 1.0*flux_r83 + 1.0*flux_r96 - 1.0*flux_r97;
    xdot[19] = 1.0*flux_r24 - 1.0*flux_r78 + 1.0*flux_r97 - 1.0*flux_r98;
    xdot[20] = -1.0*flux_r79 - 1.0*flux_r80 + 1.0*flux_r98 - 1.0*flux_r99;
    xdot[21] = 1.0*flux_r29 - 1.0*flux_r69 + 1.0*flux_r99;
    xdot[22] = -1.0*flux_r32 + 1.0*flux_r84;
    xdot[23] = 1.0*flux_r133 + 1.0*flux_r41;
    xdot[24] = 1.0*flux_r40;
    xdot[25] = 1.0*flux_r38;
    xdot[26] = 1.0*flux_r43;
    xdot[27] = 1.0*flux_r44;
    xdot[28] = 1.0*flux_r42;
    xdot[29] = 1.0*flux_r39;
    xdot[30] = 1.0*flux_r135 + 1.0*flux_r87;
    xdot[31] = 1.0*flux_r86;
    xdot[32] = 1.0*flux_r85;
    xdot[33] = 1.0*flux_r88;
    xdot[34] = 1.0*flux_r91;
    xdot[35] = 1.0*flux_r11 - 1.0*flux_r38;
    xdot[36] = 1.0*flux_r14 - 1.0*flux_r43;
    xdot[37] = 1.0*flux_r18 - 1.0*flux_r44;
    xdot[38] = -1.0*flux_r40 + 1.0*flux_r8;
    xdot[39] = -1.0*flux_r41 + 1.0*flux_r5;
    xdot[40] = 1.0*flux_r21 - 1.0*flux_r42;
    xdot[41] = 1.0*flux_r25 - 1.0*flux_r39;
    xdot[42] = 1.0*flux_r89 + 1.0*flux_r90;
    xdot[43] = 1.0*flux_r63 - 1.0*flux_r91;
    xdot[44] = 1.0*flux_r60 - 1.0*flux_r90;
    xdot[45] = 1.0*flux_r57 - 1.0*flux_r89;
    xdot[46] = 1.0*flux_r54 - 1.0*flux_r88;
    xdot[47] = 1.0*flux_r51 - 1.0*flux_r85;
    xdot[48] = 1.0*flux_r48 - 1.0*flux_r86;
    xdot[49] = 1.0*flux_r132 - 1.0*flux_r133;
    xdot[50] = 1.0*flux_r134 - 1.0*flux_r135;
    xdot[51] = 1.0*flux_r46 - 1.0*flux_r87;
    xdot[52] = 1.0*flux_r118 - 1.0*flux_r120 + 1.0*flux_r52 - 1.0*flux_r54 - 1.0*flux_r55 - 1.0*flux_r56;
    xdot[53] = -1.0*flux_r119 + 1.0*flux_r120 + 1.0*flux_r49 - 1.0*flux_r51 - 1.0*flux_r52 - 1.0*flux_r53;
    xdot[55] = 1.0*flux_r119 - 1.0*flux_r121 + 1.0*flux_r47 - 1.0*flux_r48 - 1.0*flux_r49 - 1.0*flux_r50;
    xdot[56] = 1.0*flux_r121 - 1.0*flux_r122 + 1.0*flux_r31 - 1.0*flux_r46 - 1.0*flux_r47;
    xdot[57] = 1.0*flux_r115 + 1.0*flux_r116 + 1.0*flux_r117 + 1.0*flux_r118 + 1.0*flux_r119 + 1.0*flux_r120 + 1.0*flux_r121 + 2.0*flux_r122 - 1.0*flux_r134 + 1.0*flux_r3 - 2.0*flux_r31 - 1.0*flux_r32 - 1.0*flux_r33 - 1.0*flux_r34 - 1.0*flux_r35 - 1.0*flux_r47 - 1.0*flux_r49 - 1.0*flux_r52 - 1.0*flux_r55 - 1.0*flux_r58 - 1.0*flux_r61 - 1.0*flux_r64 - 1.0*flux_r70 - 1.0*flux_r71 - 1.0*flux_r72 - 1.0*flux_r73 - 1.0*flux_r84;
    xdot[58] = 1.0*flux_r0 - 1.0*flux_r12 + 1.0*flux_r123 + 1.0*flux_r124 + 1.0*flux_r125 + 1.0*flux_r126 + 1.0*flux_r127 + 1.0*flux_r128 + 1.0*flux_r129 + 2.0*flux_r130 - 1.0*flux_r131 - 1.0*flux_r132 - 1.0*flux_r15 - 1.0*flux_r19 - 2.0*flux_r2 - 1.0*flux_r22 - 1.0*flux_r26 - 1.0*flux_r3 - 1.0*flux_r4 - 1.0*flux_r45 - 1.0*flux_r6 - 1.0*flux_r9 - 1.0*flux_r92 - 1.0*flux_r93 - 1.0*flux_r94 - 1.0*flux_r95 - 1.0*flux_r96 - 1.0*flux_r97 - 1.0*flux_r98 - 1.0*flux_r99;
    xdot[59] = 1.0*flux_r129 - 1.0*flux_r130 + 1.0*flux_r2 - 1.0*flux_r5 - 1.0*flux_r6 - 1.0*flux_r7;
    xdot[60] = -1.0*flux_r10 + 1.0*flux_r128 - 1.0*flux_r129 + 1.0*flux_r6 + 1.0*flux_r74 - 1.0*flux_r8 - 1.0*flux_r81 - 1.0*flux_r9;
    xdot[61] = 1.0*flux_r117 - 1.0*flux_r118 + 1.0*flux_r55 - 1.0*flux_r57 - 1.0*flux_r58 - 1.0*flux_r59;
    xdot[62] = 1.0*flux_r1 - 1.0*flux_r3 - 1.0*flux_r30;
    xdot[63] = 1.0*flux_r12 + 1.0*flux_r126 - 1.0*flux_r127 - 1.0*flux_r14 - 1.0*flux_r15 - 1.0*flux_r16 - 1.0*flux_r17 + 1.0*flux_r76;
    xdot[64] = -1.0*flux_r11 - 1.0*flux_r12 + 1.0*flux_r127 - 1.0*flux_r128 - 1.0*flux_r13 + 1.0*flux_r75 - 1.0*flux_r82 + 1.0*flux_r9;
    xdot[65] = 1.0*flux_r116 - 1.0*flux_r117 + 1.0*flux_r58 - 1.0*flux_r60 - 1.0*flux_r61 - 1.0*flux_r62;
    xdot[66] = 1.0*flux_r115 - 1.0*flux_r116 + 1.0*flux_r61 - 1.0*flux_r63 - 1.0*flux_r64 - 1.0*flux_r65;
    xdot[67] = -1.0*flux_r115 + 1.0*flux_r64 - 1.0*flux_r66;
    xdot[68] = -1.0*flux_r123 + 1.0*flux_r26 - 1.0*flux_r28 - 1.0*flux_r29 - 1.0*flux_r45 + 1.0*flux_r69;
    xdot[69] = 1.0*flux_r123 - 1.0*flux_r124 + 1.0*flux_r22 - 1.0*flux_r25 - 1.0*flux_r26 - 1.0*flux_r27 + 1.0*flux_r79 - 1.0*flux_r80;
    xdot[70] = 1.0*flux_r124 - 1.0*flux_r125 + 1.0*flux_r19 - 1.0*flux_r21 - 1.0*flux_r22 - 1.0*flux_r23 - 1.0*flux_r24 + 1.0*flux_r78;
    xdot[71] = 1.0*flux_r125 - 1.0*flux_r126 + 1.0*flux_r15 - 1.0*flux_r18 - 1.0*flux_r19 - 1.0*flux_r20 + 1.0*flux_r77 - 1.0*flux_r83;
    xdot[72] = -1.0*flux_r114 + 1.0*flux_r45;
    xdot[73] = -1.0*flux_r10 + 1.0*flux_r100 + 1.0*flux_r101 + 1.0*flux_r102 + 1.0*flux_r103 + 1.0*flux_r104 + 1.0*flux_r105 + 1.0*flux_r106 + 1.0*flux_r107 + 1.0*flux_r108 + 1.0*flux_r109 + 1.0*flux_r110 + 1.0*flux_r111 + 1.0*flux_r112 + 1.0*flux_r113 - 1.0*flux_r114 - 1.0*flux_r13 - 1.0*flux_r16 - 1.0*flux_r20 - 1.0*flux_r23 - 1.0*flux_r27 - 1.0*flux_r28 - 1.0*flux_r50 - 1.0*flux_r53 - 1.0*flux_r56 - 1.0*flux_r59 - 1.0*flux_r62 - 1.0*flux_r65 - 1.0*flux_r66;
    xdot[74] = 1.0*flux_r30;
    xdot[75] = 1.0*flux_r10 - 1.0*flux_r100;
    xdot[76] = -1.0*flux_r101 + 1.0*flux_r13;
    xdot[77] = -1.0*flux_r102 + 1.0*flux_r16;
    xdot[78] = -1.0*flux_r103 + 1.0*flux_r20;
    xdot[79] = -1.0*flux_r104 + 1.0*flux_r23;
    xdot[80] = -1.0*flux_r105 + 1.0*flux_r27;
    xdot[81] = -1.0*flux_r106 + 1.0*flux_r28;
    xdot[82] = -1.0*flux_r107 + 1.0*flux_r50;
    xdot[83] = -1.0*flux_r108 + 1.0*flux_r53;
    xdot[84] = -1.0*flux_r109 + 1.0*flux_r56;
    xdot[85] = -1.0*flux_r110 + 1.0*flux_r59;
    xdot[86] = -1.0*flux_r111 + 1.0*flux_r62;
    xdot[87] = -1.0*flux_r112 + 1.0*flux_r65;
    xdot[88] = -1.0*flux_r113 + 1.0*flux_r66;
    xdot[89] = 1.0*flux_r114;
}