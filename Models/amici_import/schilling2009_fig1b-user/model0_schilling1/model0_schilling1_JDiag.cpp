#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_schilling1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx1;
    JDiag[2] = -1.0*dwdx2;
    JDiag[3] = -1.0*dwdx3;
    JDiag[4] = -1.0*dwdx4;
    JDiag[5] = -1.0*dwdx5;
    JDiag[6] = -1.0*dwdx6;
    JDiag[7] = -1.0*dwdx7;
    JDiag[8] = -1.0*dwdx8 - 1.0*dwdx9;
    JDiag[9] = -1.0*dwdx10 - 1.0*dwdx11;
    JDiag[11] = -1.0*dwdx13;
    JDiag[12] = -1.0*dwdx14;
    JDiag[13] = -1.0*dwdx15;
    JDiag[14] = -1.0*dwdx16;
    JDiag[15] = -1.0*dwdx17;
    JDiag[16] = -1.0*dwdx18;
    JDiag[17] = -1.0*dwdx19;
    JDiag[18] = -1.0*dwdx20;
    JDiag[19] = -1.0*dwdx23;
    JDiag[20] = -1.0*dwdx24 - 1.0*dwdx26 - 1.0*dwdx27;
    JDiag[21] = -1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30;
    JDiag[22] = -1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33;
    JDiag[23] = -1.0*dwdx34;
    JDiag[24] = -1.0*dwdx37;
    JDiag[25] = -1.0*dwdx39 - 1.0*dwdx40;
    JDiag[26] = -1.0*dwdx41 - 1.0*dwdx42;
    JDiag[27] = -1.0*dwdx43;
    JDiag[28] = -1.0*dwdx48;
    JDiag[29] = -1.0*dwdx49;
    JDiag[30] = -1.0*dwdx51;
    JDiag[31] = -1.0*dwdx53;
    JDiag[32] = -1.0*dwdx58;
}