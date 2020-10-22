#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Leber2015(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 + 1.0*dwdx3 - 1.0*dwdx7;
    JDiag[1] = -1.0*dwdx9;
    JDiag[2] = 1.0*dwdx12 - 1.0*dwdx15;
    JDiag[3] = -1.0*dwdx16;
    JDiag[4] = -1.0*dwdx17 - 1.0*dwdx18;
    JDiag[5] = -1.0*dwdx21;
    JDiag[6] = -1.0*dwdx22;
    JDiag[7] = -0.25*dwdx27 - 0.25*dwdx28;
    JDiag[8] = -0.25*dwdx32;
    JDiag[10] = -0.25*dwdx33 - 0.25*dwdx37;
    JDiag[11] = -0.25*dwdx41;
    JDiag[12] = -14.285714285714285*dwdx42;
    JDiag[15] = -14.285714285714285*dwdx46 - 14.285714285714285*dwdx49;
    JDiag[16] = -14.285714285714285*dwdx52;
    JDiag[17] = -14.285714285714285*dwdx53 + 14.285714285714285*dwdx55;
    JDiag[18] = -1.0*dwdx57 - 1.0*dwdx58 - 1.0*dwdx59;
    JDiag[19] = -1.0*dwdx60;
    JDiag[21] = -1.0*dwdx61;
    JDiag[22] = -1.0*dwdx62;
}