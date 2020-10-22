#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_bachmann(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -2.5*dwdx1;
    JDiag[1] = -2.5*dwdx2;
    JDiag[2] = -3.6363636363636362*dwdx4;
    JDiag[3] = -3.6363636363636362*dwdx5;
    JDiag[4] = -3.6363636363636362*dwdx6;
    JDiag[5] = -3.6363636363636362*dwdx7;
    JDiag[6] = -3.6363636363636362*dwdx8;
    JDiag[8] = -2.5*dwdx10;
    JDiag[9] = -2.5*dwdx11;
    JDiag[10] = -2.5*dwdx16 - 2.5*dwdx17 - 2.5*dwdx18;
    JDiag[11] = -2.5*dwdx19;
    JDiag[12] = -2.5*dwdx20;
    JDiag[13] = -2.5*dwdx29;
    JDiag[14] = -2.5*dwdx33;
    JDiag[15] = -3.6363636363636362*dwdx35;
    JDiag[16] = -3.6363636363636362*dwdx36;
    JDiag[17] = -3.6363636363636362*dwdx37;
    JDiag[18] = -3.6363636363636362*dwdx38;
    JDiag[19] = -3.6363636363636362*dwdx39;
    JDiag[20] = -2.5*dwdx40 - 2.5*dwdx41;
    JDiag[21] = -3.6363636363636362*dwdx42;
    JDiag[22] = -2.5*dwdx49;
    JDiag[23] = -2.5*dwdx54 - 2.5*dwdx55;
    JDiag[24] = -2.5*dwdx58 - 2.5*dwdx59;
    JDiag[25] = -2.5*dwdx60;
}