#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Singh2006(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = -1.0*dwdx1 + 1.0*dwdx2;
    JDiag[2] = -1.0*dwdx3 - 2.0*dwdx4;
    JDiag[3] = 1.0*dwdx5 + 1.0*dwdx6;
    JDiag[4] = -1.0*dwdx7 - 1.0*dwdx8;
    JDiag[5] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 + 1.0*dwdx13 + 1.0*dwdx14 + 1.0*dwdx15 - 1.0*dwdx9;
    JDiag[6] = -1.0*dwdx16 + 1.0*dwdx17 - 1.0*dwdx18;
    JDiag[7] = -1.0*dwdx19 - 1.0*dwdx20;
    JDiag[8] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23;
    JDiag[9] = 1.0*dwdx24 - 1.0*dwdx25;
    JDiag[10] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29;
    JDiag[11] = 1.0*dwdx30;
    JDiag[12] = -1.0*dwdx31 - 1.0*dwdx32;
    JDiag[13] = 1.0*dwdx33 - 1.0*dwdx34;
    JDiag[14] = -1.0*dwdx35 + 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38;
    JDiag[15] = 1.0*dwdx39 - 1.0*dwdx40 - 1.0*dwdx41 + 1.0*dwdx42;
    JDiag[16] = -1.0*dwdx43 + 1.0*dwdx44 - 1.0*dwdx45;
    JDiag[17] = 1.0*dwdx46 + 1.0*dwdx47 - 1.0*dwdx48;
    JDiag[18] = -1.0*dwdx49 + 1.0*dwdx50 - 1.0*dwdx51 + 1.0*dwdx52 + 1.0*dwdx53 + 1.0*dwdx54 - 1.0*dwdx55;
    JDiag[19] = 1.0*dwdx56 + 1.0*dwdx57 - 1.0*dwdx58;
    JDiag[20] = -1.0*dwdx59 + 1.0*dwdx60;
    JDiag[21] = -1.0*dwdx61 - 1.0*dwdx62;
    JDiag[22] = 1.0*dwdx63;
    JDiag[23] = -1.0*dwdx64 + 1.0*dwdx65;
    JDiag[24] = -1.0*dwdx66 + 1.0*dwdx67 - 1.0*dwdx68;
    JDiag[25] = 0.5*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71;
    JDiag[26] = -1.0*dwdx72 - 1.0*dwdx73;
    JDiag[27] = 2.0*dwdx75 - 1.0*dwdx76 + 1.0*dwdx77;
    JDiag[28] = -1.0*dwdx78 - 1.0*dwdx79;
    JDiag[29] = 1.0*dwdx80 - 1.0*dwdx81;
    JDiag[30] = -1.0*dwdx82;
    JDiag[31] = -1.0*dwdx83;
    JDiag[32] = -1.0*dwdx84 + 1.0*dwdx85 - 1.0*dwdx86;
    JDiag[33] = -1.0*dwdx87 + 1.0*dwdx88;
    JDiag[34] = 1.0*dwdx89 - 1.0*dwdx90;
    JDiag[35] = 1.0*dwdx91 - 1.0*dwdx92;
    JDiag[36] = -1.0*dwdx93;
    JDiag[37] = 1.0*dwdx94 - 1.0*dwdx95;
    JDiag[38] = -1.0*dwdx96 + 1.0*dwdx97;
    JDiag[39] = -1.0*dwdx100 - 1.0*dwdx101 + 1.0*dwdx98 - 1.0*dwdx99;
    JDiag[40] = -1.0*dwdx102;
    JDiag[41] = -1.0*dwdx103;
    JDiag[42] = 1.0*dwdx104 - 1.0*dwdx105;
    JDiag[43] = -1.0*dwdx106 - 1.0*dwdx107;
    JDiag[44] = -1.0*dwdx108 + 1.0*dwdx109;
    JDiag[45] = -1.0*dwdx110 - 1.0*dwdx111;
    JDiag[46] = -1.0*dwdx112 - 1.0*dwdx113 + 1.0*dwdx114;
    JDiag[47] = -1.0*dwdx115;
    JDiag[48] = 1.0*dwdx116 - 1.0*dwdx117;
    JDiag[49] = 1.0*dwdx118 - 1.0*dwdx119;
    JDiag[50] = -1.0*dwdx120;
    JDiag[51] = -1.0*dwdx121 + 1.0*dwdx122;
    JDiag[52] = 1.0*dwdx123 - 1.0*dwdx124 - 1.0*dwdx125;
    JDiag[53] = 1.0*dwdx126 - 1.0*dwdx127;
    JDiag[54] = 1.0*dwdx128 - 1.0*dwdx129;
    JDiag[55] = 1.0*dwdx130 - 1.0*dwdx131;
    JDiag[56] = -1.0*dwdx132 + 1.0*dwdx133;
    JDiag[58] = -1.0*dwdx135 + 1.0*dwdx136;
    JDiag[59] = -1.0*dwdx137;
    JDiag[60] = -1.0*dwdx138;
    JDiag[61] = -1.0*dwdx140 + 1.0*dwdx141;
    JDiag[63] = -1.0*dwdx142 - 1.0*dwdx143;
    JDiag[64] = 1.0*dwdx144 - 1.0*dwdx145;
    JDiag[65] = -1.0*dwdx146 - 1.0*dwdx147 + 1.0*dwdx148;
    JDiag[66] = -1.0*dwdx149 - 1.0*dwdx150;
    JDiag[67] = -1.0*dwdx151 + 1.0*dwdx152;
}