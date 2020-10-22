#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_Yang2007(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx5 - 1.0*dwdx6;
    JDiag[1] = 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[2] = 1.0*dwdx13 - 1.0*dwdx16;
    JDiag[3] = 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx20;
    JDiag[4] = 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx26;
    JDiag[5] = 1.0*dwdx27;
    JDiag[7] = -1.0*dwdx30;
    JDiag[8] = -1.0*dwdx32;
    JDiag[11] = 1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx40;
    JDiag[12] = -1.0*dwdx43 - 1.0*dwdx44;
    JDiag[13] = 1.0*dwdx47 - 1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[14] = -1.0*dwdx52;
    JDiag[18] = 1.0*dwdx57 - 1.0*dwdx62;
    JDiag[19] = 1.0*dwdx64 - 1.0*dwdx65;
    JDiag[20] = 1.0*dwdx66;
    JDiag[21] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[22] = 1.0*dwdx75;
    JDiag[23] = 1.0*dwdx79 - 1.0*dwdx81;
    JDiag[24] = -1.0*dwdx82;
}