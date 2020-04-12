#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model5_levchenko2(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0;
    JDiag[1] = -1.0*dwdx6;
    JDiag[2] = -1.0*dwdx12;
    JDiag[3] = -1.0*dwdx17;
    JDiag[4] = -1.0*dwdx22;
    JDiag[5] = -1.0*dwdx28;
    JDiag[6] = -1.0*dwdx34;
    JDiag[7] = -1.0*dwdx40;
    JDiag[8] = -1.0*dwdx46;
    JDiag[9] = -1.0*dwdx49 - 1.0*dwdx56;
    JDiag[10] = -1.0*dwdx57 - 1.0*dwdx58;
    JDiag[11] = -1.0*dwdx59 - 1.0*dwdx60;
    JDiag[12] = -1.0*dwdx61 - 1.0*dwdx62;
    JDiag[13] = -1.0*dwdx63 - 1.0*dwdx64;
    JDiag[14] = -1.0*dwdx65 - 1.0*dwdx66;
    JDiag[15] = -1.0*dwdx67;
    JDiag[16] = -1.0*dwdx68 - 1.0*dwdx69;
    JDiag[17] = -1.0*dwdx70 - 1.0*dwdx77;
    JDiag[18] = -1.0*dwdx78 - 1.0*dwdx79;
    JDiag[19] = -1.0*dwdx80 - 1.0*dwdx81;
    JDiag[20] = -1.0*dwdx82 - 1.0*dwdx83;
    JDiag[21] = -1.0*dwdx84 - 1.0*dwdx85;
    JDiag[22] = -1.0*dwdx86 - 1.0*dwdx87;
    JDiag[23] = -1.0*dwdx88 - 1.0*dwdx89 - 1.0*dwdx90;
    JDiag[24] = -1.0*dwdx91 - 1.0*dwdx92;
    JDiag[25] = -1.0*dwdx93;
    JDiag[26] = -1.0*dwdx94;
    JDiag[27] = -1.0*dwdx95;
    JDiag[28] = -1.0*dwdx96 - 1.0*dwdx97;
    JDiag[29] = -1.0*dwdx100 - 1.0*dwdx98 - 1.0*dwdx99;
    JDiag[30] = -1.0*dwdx107 - 1.0*dwdx108;
}