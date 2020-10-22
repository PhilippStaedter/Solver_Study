#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_essunger1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx4 + 1.0*dwdx5;
    JDiag[1] = 1.0*dwdx10 - 1.0*dwdx8 - 1.0*dwdx9;
    JDiag[2] = -1.0*dwdx13 - 1.0*dwdx14;
    JDiag[3] = -1.0*dwdx15 - 1.0*dwdx16;
    JDiag[4] = -1.0*dwdx17 + 1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20;
    JDiag[5] = -1.0*dwdx21 - 1.0*dwdx22;
    JDiag[6] = -1.0*dwdx23 - 1.0*dwdx24;
}