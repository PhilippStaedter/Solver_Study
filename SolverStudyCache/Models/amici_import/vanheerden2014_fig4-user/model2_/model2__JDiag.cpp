#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[3] = -2.0*dwdx0 + 4.0*dwdx1 - 1.0*dwdx2;
    JDiag[4] = -1.0*dwdx3 - 1.0*dwdx4;
    JDiag[6] = -2.0*dwdx5 + 1.0*dwdx6;
}