#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_leloup2_Fig2A(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx1 - 1.0*dwdx2;
    JDiag[1] = -1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    JDiag[2] = 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[3] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    JDiag[4] = 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    JDiag[5] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    JDiag[6] = -1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22;
    JDiag[7] = 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25;
    JDiag[8] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JDiag[9] = -1.0*dwdx31 - 1.0*dwdx32 + 1.0*dwdx33 - 1.0*dwdx34;
    JDiag[10] = -1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38;
    JDiag[11] = -1.0*dwdx39 - 1.0*dwdx41;
    JDiag[12] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    JDiag[13] = -1.0*dwdx46 - 1.0*dwdx47;
    JDiag[14] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    JDiag[15] = -1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53;
}