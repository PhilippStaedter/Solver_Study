#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model1_liu2(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    JDiag[1] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[2] = -1.0*dwdx12;
    JDiag[4] = -1.0*dwdx13;
    JDiag[6] = -1.0*dwdx14 - 1.0*dwdx15;
    JDiag[7] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23;
    JDiag[8] = -1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx27 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32;
    JDiag[9] = 1.0*dwdx33;
    JDiag[10] = 1.0*dwdx34;
    JDiag[11] = 1.0*dwdx35;
    JDiag[12] = 1.0*dwdx36;
    JDiag[14] = -1.0*dwdx37 - 1.0*dwdx38;
    JDiag[15] = 1.0*dwdx39 - 1.0*dwdx41 - 1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44 - 1.0*dwdx45;
    JDiag[16] = 1.0*dwdx46;
    JDiag[17] = -1.0*dwdx47 - 1.0*dwdx48;
    JDiag[18] = -1.0*dwdx49;
    JDiag[19] = 1.0*dwdx50 - 1.0*dwdx51;
    JDiag[20] = 1.0*dwdx52;
    JDiag[21] = 1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57;
    JDiag[23] = 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[24] = 1.0*dwdx62;
    JDiag[25] = 1.0*dwdx65;
    JDiag[26] = 1.0*dwdx68;
    JDiag[27] = -1.0*dwdx71;
    JDiag[28] = -1.0*dwdx72 - 1.0*dwdx73;
    JDiag[29] = -1.0*dwdx74 - 1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx77 - 1.0*dwdx78;
    JDiag[30] = -1.0*dwdx79;
    JDiag[31] = 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83;
    JDiag[32] = 1.0*dwdx84;
    JDiag[33] = 1.0*dwdx87 - 1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90;
    JDiag[34] = 1.0*dwdx91 - 1.0*dwdx94;
    JDiag[35] = 1.0*dwdx95 + 1.0*dwdx96;
    JDiag[36] = -1.0*dwdx102 + 1.0*dwdx99;
    JDiag[37] = -1.0*dwdx103;
    JDiag[38] = 1.0*dwdx104;
    JDiag[39] = 1.0*dwdx105 - 1.0*dwdx106;
    JDiag[40] = 1.0*dwdx107;
}