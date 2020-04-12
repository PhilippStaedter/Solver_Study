#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_bungay2(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = -1.0*dwdx2 + 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4;
    JDiag[3] = 1.0*dwdx5 - 1.0*dwdx6;
    JDiag[4] = -1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx10;
    JDiag[6] = 1.0*dwdx11 - 1.0*dwdx12;
    JDiag[8] = 1.0*dwdx13 - 1.0*dwdx14;
    JDiag[9] = 1.0*dwdx15 - 1.0*dwdx16;
    JDiag[10] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[11] = -100.0*dwdx21 - 100.0*dwdx22 - 100.0*dwdx23 - 100.0*dwdx24 - 100.0*dwdx25 - 100.0*dwdx26 - 100.0*dwdx27 - 100.0*dwdx28 - 100.0*dwdx29 - 100.0*dwdx30 - 100.0*dwdx31 - 100.0*dwdx32;
    JDiag[12] = -1.0*dwdx33;
    JDiag[13] = 1.0*dwdx34 - 1.0*dwdx35;
    JDiag[14] = -1.0*dwdx36;
    JDiag[15] = 1.0*dwdx37 - 1.0*dwdx38;
    JDiag[16] = 1.0*dwdx39;
    JDiag[17] = 1.0*dwdx40 - 1.0*dwdx41;
    JDiag[18] = -1.0*dwdx42;
    JDiag[19] = 1.0*dwdx43 - 1.0*dwdx44;
    JDiag[20] = 1.0*dwdx45 - 1.0*dwdx46;
    JDiag[21] = 1.0*dwdx47 - 1.0*dwdx48;
    JDiag[22] = -1.0*dwdx49;
    JDiag[23] = 1.0*dwdx50 - 1.0*dwdx51 - 1.0*dwdx52;
    JDiag[24] = -1.0*dwdx53 - 1.0*dwdx54;
    JDiag[25] = -1.0*dwdx55;
    JDiag[26] = 1.0*dwdx56 - 1.0*dwdx57;
    JDiag[27] = -1.0*dwdx58;
    JDiag[28] = 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[29] = -1.0*dwdx62;
    JDiag[30] = 1.0*dwdx63 - 1.0*dwdx64;
    JDiag[31] = 1.0*dwdx65 - 1.0*dwdx66;
    JDiag[32] = 1.0*dwdx67 - 1.0*dwdx68;
    JDiag[33] = -1.0*dwdx69;
    JDiag[34] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72 - 1.0*dwdx73;
    JDiag[35] = 1.0*dwdx74 - 1.0*dwdx75;
    JDiag[36] = -1.0*dwdx76;
    JDiag[37] = 1.0*dwdx77 - 1.0*dwdx78 - 1.0*dwdx79;
    JDiag[38] = -1.0*dwdx80;
    JDiag[39] = 1.0*dwdx81;
    JDiag[40] = -1.0*dwdx82;
    JDiag[41] = 1.0*dwdx83 - 1.0*dwdx84;
    JDiag[43] = 1.0*dwdx85 - 1.0*dwdx86;
    JDiag[44] = -1.0*dwdx87 - 1.0*dwdx88 + 1.0*dwdx89;
    JDiag[45] = 1.0*dwdx90 - 1.0*dwdx91;
    JDiag[46] = -1.0*dwdx92 - 1.0*dwdx93 - 1.0*dwdx94 - 1.0*dwdx95;
    JDiag[47] = -1.0*dwdx100 + 1.0*dwdx96 - 1.0*dwdx97 - 1.0*dwdx98 - 1.0*dwdx99;
    JDiag[50] = -1.0*dwdx101 - 1.0*dwdx102;
    JDiag[52] = -1.0*dwdx103 - 1.0*dwdx104;
    JDiag[53] = 1.0*dwdx105 - 1.0*dwdx106 - 1.0*dwdx107;
}