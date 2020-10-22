#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_kolodkin5(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[1] = 1.0*dwdx2 - 1.0*dwdx3;
    JDiag[2] = -1.0*dwdx4 + 3.4444444444444402*dwdx5;
    JDiag[3] = 1.0*dwdx6 - 1.0*dwdx7;
    JDiag[4] = 3.4444444444444402*dwdx10 - 1.0*dwdx8 + 1.0*dwdx9;
    JDiag[5] = -1.0*dwdx11 - 1.0*dwdx12;
    JDiag[6] = -1.0*dwdx13;
    JDiag[7] = 1.0*dwdx14;
}