#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_hald(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0 - 1.0*dwdx1 - 59.0*dwdx2;
    JDiag[1] = -1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5;
    JDiag[2] = 2.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8;
    JDiag[3] = -1.0*dwdx9;
    JDiag[4] = -1.0*dwdx11 - 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 + 1.0*dwdx16;
    JDiag[5] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[6] = 1.0*dwdx21;
    JDiag[7] = 1.0*dwdx22 - 1.0*dwdx23;
    JDiag[8] = -1.0*dwdx24 - 59.0*dwdx25;
    JDiag[9] = 1.0*dwdx26;
    JDiag[10] = -1.0*dwdx27 + 1.0*dwdx28;
    JDiag[11] = -1.0*dwdx29;
    JDiag[12] = -1.0*dwdx32 - 1.0*dwdx33;
    JDiag[13] = 1.0*dwdx34 - 1.0*dwdx35 + 1.0*dwdx36;
    JDiag[14] = 59.0*dwdx37 - 1.0*dwdx38;
    JDiag[15] = -1.0*dwdx39;
    JDiag[16] = -59.0*dwdx40;
    JDiag[17] = 1.0*dwdx41;
    JDiag[18] = -1.0*dwdx42 - 1.0*dwdx43 - 59.0*dwdx44;
    JDiag[19] = -1.0*dwdx45 - 1.0*dwdx46 + 1.0*dwdx47;
    JDiag[20] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51;
    JDiag[21] = 1.0*dwdx52 + 1.0*dwdx53 - 1.0*dwdx54;
    JDiag[22] = -59.0*dwdx55;
    JDiag[23] = 1.0*dwdx56;
    JDiag[24] = -1.0*dwdx57 + 1.0*dwdx58;
    JDiag[25] = -1.0*dwdx59 - 1.0*dwdx60 - 1.0*dwdx61;
    JDiag[26] = 1.0*dwdx62;
    JDiag[30] = 1.0*dwdx63;
}