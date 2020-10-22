#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_essunger3(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx3;
    JDiag[1] = -1.0*dwdx6 - 1.0*dwdx8;
    JDiag[2] = -1.0*dwdx10 - 1.0*dwdx9;
    JDiag[3] = -1.0*dwdx11 - 1.0*dwdx12;
    JDiag[4] = 1.0*dwdx13 + 1.0*dwdx16;
    JDiag[5] = -1.0*dwdx17 - 1.0*dwdx18;
    JDiag[6] = -1.0*dwdx19 - 1.0*dwdx20;
}