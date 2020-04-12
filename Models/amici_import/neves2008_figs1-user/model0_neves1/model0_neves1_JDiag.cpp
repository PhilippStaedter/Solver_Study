#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_neves1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = 5.0*dwdx0;
    JDiag[2] = -5.0*dwdx2;
    JDiag[4] = -1.0*dwdx4 - 1.0*dwdx5;
    JDiag[5] = 5.0*dwdx6 - 5.0*dwdx7;
    JDiag[6] = -5.0*dwdx8 - 5.0*dwdx9;
    JDiag[7] = -1.0*dwdx10;
    JDiag[8] = -1.0*dwdx12;
    JDiag[9] = 1.0*dwdx14;
    JDiag[10] = -1.0*dwdx16;
    JDiag[11] = 1.0*dwdx17 - 1.0*dwdx18;
    JDiag[12] = -1.0*dwdx19 - 1.0*dwdx20 + 1.0*dwdx21;
    JDiag[13] = -1.0*dwdx22 - 1.0*dwdx23 + 1.0*dwdx24;
    JDiag[14] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27;
    JDiag[15] = -1.0*dwdx28;
    JDiag[16] = -1.0*dwdx30;
    JDiag[17] = -1.0*dwdx31;
    JDiag[18] = -1.0*dwdx33;
    JDiag[19] = -1.0*dwdx35;
    JDiag[21] = 1.0*dwdx37;
    JDiag[24] = -1.0*dwdx46;
    JDiag[26] = -1.0*dwdx48;
    JDiag[27] = -1.0*dwdx50;
    JDiag[28] = 1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53;
    JDiag[29] = -1.0*dwdx54 + 1.0*dwdx55;
    JDiag[30] = 1.0*dwdx56 - 1.0*dwdx57;
    JDiag[31] = -1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63 - 1.0*dwdx64;
    JDiag[32] = 1.0*dwdx65 - 1.0*dwdx66;
    JDiag[33] = 5.0*dwdx67 - 5.0*dwdx68 + 5.0*dwdx69;
    JDiag[34] = -5.0*dwdx70 - 5.0*dwdx71 - 5.0*dwdx72 + 5.0*dwdx73 + 5.0*dwdx74;
    JDiag[36] = -9.0000000000000089*dwdx75 - 9.0000000000000089*dwdx76;
}