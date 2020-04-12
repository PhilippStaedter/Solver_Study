#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model0_orfao1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx1;
    JDiag[1] = -1.0*dwdx3 - 1.0*dwdx4;
    JDiag[4] = -1.0*dwdx5;
    JDiag[5] = 1.0*dwdx6;
    JDiag[7] = -1.0*dwdx9;
    JDiag[8] = -1.0*dwdx10;
    JDiag[9] = -1.0*dwdx11;
    JDiag[10] = -1.0*dwdx12 - 1.0*dwdx13;
}