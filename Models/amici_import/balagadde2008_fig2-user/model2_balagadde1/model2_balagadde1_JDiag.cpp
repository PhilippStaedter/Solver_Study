#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_balagadde1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx1;
    JDiag[1] = -1.0*dwdx3;
    JDiag[2] = 1.0*dwdx4 - 1.0*dwdx5;
    JDiag[3] = -1.0*dwdx10 + 1.0*dwdx9;
}