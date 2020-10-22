#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model3_hockin1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4;
    JDiag[1] = -1.0*dwdx5 - 1.0*dwdx6;
    JDiag[2] = -1.0*dwdx9;
    JDiag[4] = -1.0*dwdx11;
    JDiag[5] = -1.0*dwdx12 - 1.0*dwdx13;
    JDiag[7] = 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    JDiag[8] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    JDiag[9] = -1.0*dwdx20 - 1.0*dwdx21;
    JDiag[10] = -1.0*dwdx22 - 1.0*dwdx23;
    JDiag[11] = 1.0*dwdx24;
    JDiag[12] = 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    JDiag[14] = 1.0*dwdx32 - 1.0*dwdx33;
    JDiag[15] = 1.0*dwdx34 - 1.0*dwdx35;
    JDiag[16] = -1.0*dwdx36 + 1.0*dwdx37;
    JDiag[17] = 1.0*dwdx38;
    JDiag[18] = -1.0*dwdx39;
    JDiag[19] = -1.0*dwdx40 - 1.0*dwdx41 - 1.0*dwdx42 - 1.0*dwdx43;
    JDiag[20] = -1.0*dwdx44;
    JDiag[21] = -1.0*dwdx45 - 1.0*dwdx46;
    JDiag[22] = 1.0*dwdx47;
    JDiag[23] = 1.0*dwdx48;
    JDiag[24] = -1.0*dwdx49;
    JDiag[25] = -1.0*dwdx50;
    JDiag[26] = -1.0*dwdx51 - 1.0*dwdx52;
    JDiag[27] = -1.0*dwdx53 - 1.0*dwdx54 - 1.0*dwdx55 - 1.0*dwdx57;
    JDiag[29] = 1.0*dwdx59 - 1.0*dwdx60;
    JDiag[30] = 1.0*dwdx61 - 1.0*dwdx62;
    JDiag[31] = 1.0*dwdx64 - 1.0*dwdx65;
    JDiag[32] = -1.0*dwdx66 - 1.0*dwdx67;
}