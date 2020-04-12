#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_ODea2007(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    JDiag[2] = -1.0*dwdx6 - 1.0*dwdx7 + 1.0*dwdx8;
    JDiag[3] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx9;
    JDiag[4] = 1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14;
    JDiag[5] = 1.0*dwdx15 - 1.0*dwdx16 - 1.0*dwdx17;
    JDiag[6] = 1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[7] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26 - 1.0*dwdx27;
    JDiag[8] = -1.0*dwdx28 - 1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31 - 1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34;
    JDiag[9] = -1.0*dwdx35 - 1.0*dwdx36 + 1.0*dwdx37 - 1.0*dwdx38;
    JDiag[10] = -1.0*dwdx40 - 1.0*dwdx41 + 1.0*dwdx42;
    JDiag[11] = -1.0*dwdx43 + 1.0*dwdx44 + 1.0*dwdx45;
    JDiag[12] = -1.0*dwdx46 - 1.0*dwdx47 + 1.0*dwdx48;
    JDiag[13] = -1.0*dwdx49 - 1.0*dwdx50 + 1.0*dwdx51;
    JDiag[14] = 1.0*dwdx52 - 1.0*dwdx53 - 1.0*dwdx54;
    JDiag[15] = -1.0*dwdx55 - 1.0*dwdx56 - 1.0*dwdx57 - 1.0*dwdx58;
    JDiag[16] = -1.0*dwdx60;
    JDiag[17] = -1.0*dwdx61;
    JDiag[18] = -1.0*dwdx63 - 1.0*dwdx64 - 1.0*dwdx65 - 1.0*dwdx66;
    JDiag[19] = -1.0*dwdx67 - 1.0*dwdx68 + 1.0*dwdx69;
    JDiag[20] = 1.0*dwdx70 - 1.0*dwdx71 - 1.0*dwdx72;
    JDiag[21] = 1.0*dwdx73 - 1.0*dwdx74 - 1.0*dwdx75;
    JDiag[22] = 1.0*dwdx76 + 1.0*dwdx77 - 1.0*dwdx78;
    JDiag[23] = 1.0*dwdx79 - 1.0*dwdx80 - 1.0*dwdx81;
}