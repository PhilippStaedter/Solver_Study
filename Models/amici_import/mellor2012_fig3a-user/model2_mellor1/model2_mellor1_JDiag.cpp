#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JDiag_model2_mellor1(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -0.001*dwdx0 - 0.001*dwdx1 - 0.001*dwdx2;
    JDiag[7] = -0.001*dwdx3;
    JDiag[8] = -0.001*dwdx4;
}