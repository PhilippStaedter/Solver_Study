#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_dupont2_Fig4B(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = 1.0*dwdx0;
    JDiag[2] = -1.0*dwdx2 - 1.0*dwdx3;
    JDiag[3] = -1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7;
}