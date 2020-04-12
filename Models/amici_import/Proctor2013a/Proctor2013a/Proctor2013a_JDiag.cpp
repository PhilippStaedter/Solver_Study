#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Proctor2013a(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx1;
    JDiag[1] = -1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4;
    JDiag[3] = -1.0*dwdx6 - 1.0*dwdx7;
    JDiag[4] = -1.0*dwdx8;
    JDiag[5] = -1.0*dwdx11;
    JDiag[6] = -1.0*dwdx13;
    JDiag[7] = -1.0*dwdx14 - 1.0*dwdx17;
    JDiag[8] = -1.0*dwdx18;
    JDiag[9] = -1.0*dwdx19;
    JDiag[10] = -1.0*dwdx20 - 1.0*dwdx21;
    JDiag[11] = -1.0*dwdx23;
    JDiag[12] = -1.0*dwdx24 - 1.0*dwdx25;
    JDiag[13] = -1.0*dwdx29;
    JDiag[14] = -1.0*dwdx31;
    JDiag[15] = -1.0*dwdx33;
    JDiag[16] = -1.0*dwdx35;
    JDiag[17] = -1.0*dwdx37;
    JDiag[18] = -1.0*dwdx38;
    JDiag[19] = -1.0*dwdx39 - 1.0*dwdx40;
    JDiag[20] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    JDiag[21] = -1.0*dwdx45 - 1.0*dwdx46;
    JDiag[22] = -1.0*dwdx47;
    JDiag[23] = -1.0*dwdx48;
    JDiag[24] = -1.0*dwdx52;
    JDiag[25] = -1.0*dwdx53 - 1.0*dwdx54;
    JDiag[26] = -1.0*dwdx55;
    JDiag[27] = -1.0*dwdx57;
    JDiag[28] = -1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60;
    JDiag[29] = -1.0*dwdx62;
    JDiag[30] = -1.0*dwdx64;
    JDiag[31] = -1.0*dwdx65 - 1.0*dwdx66;
    JDiag[32] = -1.0*dwdx67;
    JDiag[33] = -1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx71;
    JDiag[34] = -1.0*dwdx72;
    JDiag[35] = -1.0*dwdx73;
    JDiag[37] = -1.0*dwdx74;
    JDiag[40] = -1.0*dwdx75 - 1.0*dwdx76;
    JDiag[41] = -1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79;
    JDiag[42] = -1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx83;
    JDiag[43] = -1.0*dwdx84;
    JDiag[44] = -1.0*dwdx85;
    JDiag[45] = -1.0*dwdx87 - 1.0*dwdx89 - 1.0*dwdx90;
    JDiag[46] = -1.0*dwdx91;
    JDiag[47] = -1.0*dwdx92;
    JDiag[48] = -1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx96;
    JDiag[49] = -1.0*dwdx97;
    JDiag[50] = -1.0*dwdx98;
    JDiag[51] = -1.0*dwdx100 - 1.0*dwdx101 - 1.0*dwdx99;
    JDiag[52] = -1.0*dwdx102 - 1.0*dwdx103 - 1.0*dwdx104 - 1.0*dwdx105 - 1.0*dwdx106;
    JDiag[53] = -1.0*dwdx107 - 1.0*dwdx108 - 1.0*dwdx109 - 1.0*dwdx110 - 1.0*dwdx111;
    JDiag[54] = -1.0*dwdx112 - 1.0*dwdx113;
    JDiag[55] = -1.0*dwdx114;
    JDiag[56] = -1.0*dwdx115 - 1.0*dwdx116;
    JDiag[57] = -1.0*dwdx117;
    JDiag[58] = -1.0*dwdx118;
    JDiag[59] = -1.0*dwdx119;
    JDiag[60] = -1.0*dwdx121;
    JDiag[61] = -1.0*dwdx122;
    JDiag[62] = -1.0*dwdx123 - 1.0*dwdx124;
    JDiag[63] = -1.0*dwdx125;
    JDiag[64] = -1.0*dwdx135;
    JDiag[65] = -1.0*dwdx140 - 1.0*dwdx141 - 1.0*dwdx142;
    JDiag[66] = -1.0*dwdx143 - 2.0*dwdx144 - 1.0*dwdx145;
    JDiag[67] = -1.0*dwdx146;
    JDiag[68] = -1.0*dwdx155 - 1.0*dwdx156;
    JDiag[69] = -1.0*dwdx157;
    JDiag[70] = -1.0*dwdx158;
    JDiag[71] = -1.0*dwdx159 - 1.0*dwdx160;
    JDiag[72] = -1.0*dwdx166;
}