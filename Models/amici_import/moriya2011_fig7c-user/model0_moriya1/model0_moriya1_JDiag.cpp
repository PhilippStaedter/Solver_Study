#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_moriya1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 + 1.0*dwdx3;
    JDiag[2] = 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7;
    JDiag[3] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx12 - 1.0*dwdx13 - 1.0*dwdx14 + 1.0*dwdx15;
    JDiag[6] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19 - 1.0*dwdx20;
    JDiag[9] = 1.0*dwdx23 - 1.0*dwdx25;
    JDiag[10] = 1.0*dwdx39 - 1.0*dwdx40;
    JDiag[12] = 1.0*dwdx41 - 1.0*dwdx42;
    JDiag[16] = 1.0*dwdx45 - 1.0*dwdx53 - 1.0*dwdx58 - 1.0*dwdx59;
    JDiag[18] = 1.0*dwdx63 - 1.0*dwdx68 - 1.0*dwdx69 - 1.0*dwdx70;
    JDiag[20] = -1.0*dwdx74 - 1.0*dwdx75 + 1.0*dwdx77 - 1.0*dwdx79;
    JDiag[24] = -1.0*dwdx86 + 1.0*dwdx87 - 1.0*dwdx88 - 1.0*dwdx90 + 1.0*dwdx92;
    JDiag[26] = -1.0*dwdx102 + 1.0*dwdx103;
    JDiag[27] = -1.0*dwdx107;
    JDiag[30] = -1.0*dwdx115;
    JDiag[35] = 1.0*dwdx119 - 1.0*dwdx120;
    JDiag[36] = -1.0*dwdx125;
    JDiag[38] = 1.0*dwdx129 - 1.0*dwdx130;
    JDiag[39] = -1.0*dwdx135;
    JDiag[44] = 1.0*dwdx140 - 1.0*dwdx141;
    JDiag[45] = -1.0*dwdx143;
}