#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_kofahl(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx1 - 1.0*dwdx2;
    JDiag[3] = -1.0*dwdx5;
    JDiag[4] = -1.0*dwdx6 - 1.0*dwdx7;
    JDiag[5] = -1.0*dwdx10 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[7] = -1.0*dwdx11 - 1.0*dwdx12;
    JDiag[8] = -1.0*dwdx13 - 1.0*dwdx14;
    JDiag[9] = -1.0*dwdx17;
    JDiag[10] = -1.0*dwdx18 - 1.0*dwdx19;
    JDiag[11] = -1.0*dwdx20;
    JDiag[12] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23;
    JDiag[13] = -1.0*dwdx24;
    JDiag[14] = -1.0*dwdx26;
    JDiag[15] = -1.0*dwdx27;
    JDiag[16] = -1.0*dwdx28;
    JDiag[17] = -1.0*dwdx30 - 1.0*dwdx31;
    JDiag[18] = -1.0*dwdx32;
    JDiag[19] = -1.0*dwdx33 - 1.0*dwdx34;
    JDiag[20] = -1.0*dwdx36;
    JDiag[21] = -1.0*dwdx37;
    JDiag[22] = -1.0*dwdx38;
    JDiag[23] = -1.0*dwdx40 - 1.0*dwdx41;
    JDiag[24] = -1.0*dwdx42 - 1.0*dwdx43;
    JDiag[25] = -1.0*dwdx44 - 1.0*dwdx45;
    JDiag[26] = -1.0*dwdx46 - 1.0*dwdx47;
    JDiag[27] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[28] = -1.0*dwdx51 - 1.0*dwdx52;
    JDiag[29] = -1.0*dwdx53 - 1.0*dwdx54;
    JDiag[30] = -1.0*dwdx55 - 1.0*dwdx56;
    JDiag[31] = -1.0*dwdx57;
    JDiag[32] = -1.0*dwdx58 - 1.0*dwdx59;
    JDiag[33] = -1.0*dwdx60 - 1.0*dwdx61;
    JDiag[34] = -1.0*dwdx62;
    JDiag[35] = -1.0*dwdx63;
}