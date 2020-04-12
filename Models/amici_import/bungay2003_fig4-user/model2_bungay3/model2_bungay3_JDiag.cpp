#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_bungay3(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = 1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4 - 1.0*dwdx5 + 1.0*dwdx6;
    JDiag[3] = -1.0*dwdx7;
    JDiag[4] = 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13;
    JDiag[6] = -1.0*dwdx14;
    JDiag[7] = 1.0*dwdx15 - 1.0*dwdx16;
    JDiag[9] = 1.0*dwdx17 - 1.0*dwdx18;
    JDiag[10] = 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[11] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25;
    JDiag[12] = -1.0*dwdx26;
    JDiag[13] = 1.0*dwdx27 - 1.0*dwdx28 - 1.0*dwdx29;
    JDiag[15] = 1.0*dwdx30 - 1.0*dwdx31;
    JDiag[16] = 1.0*dwdx32 - 1.0*dwdx33;
    JDiag[17] = -1.0*dwdx34 - 1.0*dwdx35;
    JDiag[18] = 1.0*dwdx36 - 1.0*dwdx37;
    JDiag[19] = -100.0*dwdx38 - 100.0*dwdx39 - 100.0*dwdx40 - 100.0*dwdx41 - 100.0*dwdx42 - 100.0*dwdx43 - 100.0*dwdx44 - 100.0*dwdx45 - 100.0*dwdx46 - 100.0*dwdx47 - 100.0*dwdx48 - 100.0*dwdx49 - 100.0*dwdx50 - 100.0*dwdx51 - 100.0*dwdx52 - 100.0*dwdx53 - 100.0*dwdx54;
    JDiag[20] = -1.0*dwdx55;
    JDiag[21] = 1.0*dwdx56 - 1.0*dwdx57;
    JDiag[22] = -1.0*dwdx58;
    JDiag[23] = 1.0*dwdx59 - 1.0*dwdx60;
    JDiag[24] = 1.0*dwdx61;
    JDiag[25] = 1.0*dwdx62 - 1.0*dwdx63;
    JDiag[26] = -1.0*dwdx64;
    JDiag[27] = 1.0*dwdx65 - 1.0*dwdx66;
    JDiag[28] = 1.0*dwdx67 - 1.0*dwdx68;
    JDiag[29] = 1.0*dwdx69 - 1.0*dwdx70;
    JDiag[31] = 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[32] = -1.0*dwdx73;
    JDiag[33] = 1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77;
    JDiag[34] = -1.0*dwdx78 - 1.0*dwdx79;
    JDiag[35] = -1.0*dwdx80;
    JDiag[36] = 1.0*dwdx81 - 1.0*dwdx82;
    JDiag[37] = 1.0*dwdx83 - 1.0*dwdx84;
    JDiag[38] = -1.0*dwdx85;
    JDiag[39] = 1.0*dwdx86 - 1.0*dwdx87 - 1.0*dwdx88 - 1.0*dwdx89;
    JDiag[40] = 1.0*dwdx90 - 1.0*dwdx91;
    JDiag[41] = -1.0*dwdx92;
    JDiag[42] = 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95;
    JDiag[43] = -1.0*dwdx96;
    JDiag[44] = 1.0*dwdx97;
    JDiag[45] = 1.0*dwdx98 - 1.0*dwdx99;
    JDiag[46] = -1.0*dwdx100;
    JDiag[47] = 1.0*dwdx101 - 1.0*dwdx102 - 1.0*dwdx103;
    JDiag[48] = -1.0*dwdx104;
    JDiag[49] = 1.0*dwdx105 - 1.0*dwdx106;
    JDiag[50] = 1.0*dwdx107 - 1.0*dwdx108;
    JDiag[51] = 1.0*dwdx109 - 1.0*dwdx110;
    JDiag[52] = -1.0*dwdx111;
    JDiag[53] = 1.0*dwdx112 - 1.0*dwdx113 - 1.0*dwdx114 - 1.0*dwdx115;
    JDiag[54] = 1.0*dwdx116 - 1.0*dwdx117;
    JDiag[55] = -1.0*dwdx118;
    JDiag[56] = 1.0*dwdx119 - 1.0*dwdx120 - 1.0*dwdx121;
    JDiag[57] = -1.0*dwdx122;
    JDiag[58] = 1.0*dwdx123;
    JDiag[59] = 1.0*dwdx124 - 1.0*dwdx125;
    JDiag[60] = -1.0*dwdx126;
    JDiag[61] = 1.0*dwdx127 - 1.0*dwdx128;
    JDiag[62] = -1.0*dwdx129;
    JDiag[63] = -1.0*dwdx130;
    JDiag[64] = 1.0*dwdx131 - 1.0*dwdx132 - 1.0*dwdx133;
    JDiag[66] = 1.0*dwdx134 - 1.0*dwdx135;
    JDiag[67] = -1.0*dwdx136 - 1.0*dwdx137 + 1.0*dwdx138;
    JDiag[68] = 1.0*dwdx139 - 1.0*dwdx140;
    JDiag[69] = -1.0*dwdx141 - 1.0*dwdx142 - 1.0*dwdx143;
    JDiag[70] = 1.0*dwdx144 - 1.0*dwdx145 - 1.0*dwdx146 - 1.0*dwdx147 - 1.0*dwdx148 - 1.0*dwdx149;
    JDiag[72] = -1.0*dwdx150 - 1.0*dwdx151;
    JDiag[73] = 1.0*dwdx152 - 1.0*dwdx153 - 1.0*dwdx154 - 1.0*dwdx155;
}