#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Liu2011(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = -1.0*dwdx2;
    JDiag[2] = 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    JDiag[3] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx15 - 1.0*dwdx16;
    JDiag[6] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24;
    JDiag[7] = -1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JDiag[8] = 1.0*dwdx29;
    JDiag[9] = -1.0*dwdx32;
    JDiag[11] = 1.0*dwdx33 - 1.0*dwdx35 - 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38 - 1.0*dwdx39;
    JDiag[12] = -1.0*dwdx40;
    JDiag[14] = -1.0*dwdx41 - 1.0*dwdx42;
    JDiag[15] = 1.0*dwdx43;
    JDiag[16] = -1.0*dwdx44 - 1.0*dwdx45 - 1.0*dwdx46 - 1.0*dwdx47 - 1.0*dwdx48;
    JDiag[17] = 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[18] = -1.0*dwdx51;
    JDiag[19] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JDiag[20] = -1.0*dwdx55 - 1.0*dwdx56;
    JDiag[21] = 1.0*dwdx57;
    JDiag[22] = 1.0*dwdx60 - 1.0*dwdx61 - 1.0*dwdx62 - 1.0*dwdx63;
    JDiag[23] = 1.0*dwdx64 - 1.0*dwdx67;
    JDiag[24] = 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70 - 1.0*dwdx71;
    JDiag[25] = 1.0*dwdx72;
    JDiag[26] = -1.0*dwdx75 - 1.0*dwdx76 - 1.0*dwdx78 - 1.0*dwdx80 - 1.0*dwdx81 - 1.0*dwdx82 - 1.0*dwdx83;
    JDiag[27] = 1.0*dwdx84;
    JDiag[28] = 1.0*dwdx85;
    JDiag[30] = 1.0*dwdx86;
    JDiag[31] = 1.0*dwdx87;
    JDiag[32] = 1.0*dwdx88;
    JDiag[33] = 1.0*dwdx89 - 1.0*dwdx92;
    JDiag[34] = 1.0*dwdx93;
    JDiag[35] = 1.0*dwdx94;
    JDiag[36] = 1.0*dwdx97 + 1.0*dwdx98;
    JDiag[38] = 1.0*dwdx101 - 1.0*dwdx102;
    JDiag[39] = -1.0*dwdx103;
    JDiag[40] = 1.0*dwdx104;
    JDiag[41] = -1.0*dwdx107;
}