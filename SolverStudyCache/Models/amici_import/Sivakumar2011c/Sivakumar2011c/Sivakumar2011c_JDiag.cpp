#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Sivakumar2011c(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 + 1.0*dwdx1;
    JDiag[1] = 1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4;
    JDiag[3] = 1.0*dwdx5;
    JDiag[4] = -1.0*dwdx6;
    JDiag[5] = 1.0*dwdx7;
    JDiag[6] = -1.0*dwdx10;
    JDiag[7] = -1.0*dwdx11 - 1.0*dwdx12 - 3.0*dwdx13;
    JDiag[8] = 1.0*dwdx14 + 1.0*dwdx15 + 3.0*dwdx16;
    JDiag[9] = 1.0*dwdx17;
    JDiag[10] = -1.0*dwdx19;
    JDiag[11] = -1.0*dwdx20;
    JDiag[12] = 1.0*dwdx21;
    JDiag[13] = -1.0*dwdx22;
    JDiag[14] = -1.0*dwdx23;
    JDiag[15] = -1.0*dwdx24;
    JDiag[16] = 1.0*dwdx25;
    JDiag[17] = 1.0*dwdx26;
    JDiag[18] = -1.0*dwdx27;
    JDiag[19] = -1.0*dwdx28;
    JDiag[20] = -1.0*dwdx29;
    JDiag[21] = -1.0*dwdx30;
    JDiag[23] = -1.0*dwdx31;
    JDiag[24] = -1.0*dwdx32;
    JDiag[25] = -1.0*dwdx33 - 1.0*dwdx34;
    JDiag[26] = 1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37;
    JDiag[27] = 1.0*dwdx38 - 1.0*dwdx39;
    JDiag[28] = 1.0*dwdx40 - 1.0*dwdx41;
    JDiag[29] = 1.0*dwdx42 - 1.0*dwdx43;
    JDiag[30] = 1.0*dwdx44 - 1.0*dwdx45;
    JDiag[31] = -1.0*dwdx46 + 1.0*dwdx47;
    JDiag[32] = -1.0*dwdx48 + 1.0*dwdx49;
    JDiag[33] = 1.0*dwdx50 - 1.0*dwdx51;
    JDiag[34] = 1.0*dwdx52 - 1.0*dwdx53;
    JDiag[35] = -1.0*dwdx54 + 1.0*dwdx55;
    JDiag[36] = -1.0*dwdx56 - 1.0*dwdx57;
    JDiag[37] = 1.0*dwdx58 - 1.0*dwdx59;
    JDiag[38] = 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[39] = 1.0*dwdx62 - 1.0*dwdx63;
    JDiag[40] = 1.0*dwdx64 - 1.0*dwdx65;
    JDiag[41] = -1.0*dwdx66 - 1.0*dwdx67 + 1.0*dwdx68;
    JDiag[42] = 1.0*dwdx69 - 1.0*dwdx70;
    JDiag[43] = 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[44] = -1.0*dwdx73 + 1.0*dwdx74;
    JDiag[45] = 1.0*dwdx75 - 1.0*dwdx76;
    JDiag[46] = 1.0*dwdx77 - 1.0*dwdx78;
    JDiag[47] = 1.0*dwdx79 - 1.0*dwdx80;
    JDiag[49] = -1.0*dwdx82 + 1.0*dwdx83;
}