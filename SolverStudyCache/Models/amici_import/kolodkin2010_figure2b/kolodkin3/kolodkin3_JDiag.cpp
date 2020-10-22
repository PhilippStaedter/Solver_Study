#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_kolodkin3(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = 1.0*dwdx2 - 1.0*dwdx3 - 3.4444444444400002*dwdx4;
    JDiag[2] = -1.0*dwdx5 + 1.0*dwdx6;
    JDiag[3] = 1.0*dwdx7 + 1.0*dwdx8;
    JDiag[4] = 1.0*dwdx10 - 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx11 - 1.0*dwdx12;
    JDiag[6] = -1.0*dwdx13 + 1.0*dwdx14;
    JDiag[7] = -1.0*dwdx15;
    JDiag[8] = -1.0*dwdx16;
    JDiag[9] = 1.0*dwdx17;
}