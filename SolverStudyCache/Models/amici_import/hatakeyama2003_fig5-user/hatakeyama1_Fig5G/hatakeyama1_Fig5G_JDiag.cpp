#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_hatakeyama1_Fig5G(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5;
    JDiag[2] = 1.0*dwdx6 - 1.0*dwdx7;
    JDiag[3] = -1.0*dwdx13;
    JDiag[5] = -1.0*dwdx15;
    JDiag[6] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
    JDiag[7] = -1.0*dwdx22;
    JDiag[8] = -1.0*dwdx23 + 1.0*dwdx24;
    JDiag[9] = -1.0*dwdx25;
    JDiag[10] = -1.0*dwdx26;
    JDiag[11] = 1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 + 1.0*dwdx31;
    JDiag[12] = -1.0*dwdx35;
    JDiag[14] = -1.0*dwdx42;
    JDiag[15] = 1.0*dwdx43 - 1.0*dwdx44;
    JDiag[16] = -1.0*dwdx46 - 1.0*dwdx47;
    JDiag[18] = -1.0*dwdx52;
    JDiag[19] = -1.0*dwdx53;
    JDiag[20] = 1.0*dwdx54 - 2.0*dwdx55;
    JDiag[21] = 1.0*dwdx56 - 1.0*dwdx57;
    JDiag[22] = 1.0*dwdx58 - 1.0*dwdx59 - 1.0*dwdx60 + 1.0*dwdx61 - 1.0*dwdx62 + 1.0*dwdx63 - 1.0*dwdx64;
    JDiag[23] = 1.0*dwdx65 - 1.0*dwdx66;
    JDiag[24] = 1.0*dwdx67 - 1.0*dwdx68;
    JDiag[25] = 1.0*dwdx69 - 1.0*dwdx70;
    JDiag[26] = 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[27] = 1.0*dwdx73 - 1.0*dwdx74;
    JDiag[28] = -1.0*dwdx75;
    JDiag[29] = -1.0*dwdx76;
    JDiag[30] = -1.0*dwdx79;
    JDiag[31] = -1.0*dwdx80;
    JDiag[32] = 1.0*dwdx82 - 1.0*dwdx83;
    JDiag[33] = 1.0*dwdx85 - 1.0*dwdx86;
    JDiag[34] = -1.0*dwdx87;
    JDiag[35] = 1.0*dwdx88;
}