#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_xie1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2;
    JDiag[1] = -1.0*dwdx7 - 1.0*dwdx8;
    JDiag[2] = -1.0*dwdx10 - 1.0*dwdx9;
    JDiag[4] = -1.0*dwdx12;
    JDiag[5] = -1.0*dwdx14 - 1.0*dwdx15;
    JDiag[6] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    JDiag[7] = -1.0*dwdx19 - 1.0*dwdx20;
    JDiag[8] = -1.0*dwdx21;
    JDiag[9] = -1.0*dwdx24;
    JDiag[11] = -1.0*dwdx26;
    JDiag[13] = -1.0*dwdx29;
    JDiag[15] = -1.0*dwdx33 + 1.0*dwdx34;
    JDiag[16] = 1.0*dwdx36 - 1.0*dwdx37;
    JDiag[17] = 1.0*dwdx39 - 1.0*dwdx40;
    JDiag[18] = -1.0*dwdx42 + 1.0*dwdx43;
    JDiag[19] = 1.0*dwdx46 - 1.0*dwdx47;
    JDiag[20] = 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[21] = -1.0*dwdx52;
    JDiag[23] = -1.0*dwdx55;
}