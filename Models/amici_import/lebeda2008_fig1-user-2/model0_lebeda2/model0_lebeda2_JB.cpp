#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_lebeda2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[3] = -1.0*dwdx0;
    JB[4] = -1.0*dwdx1;
    JB[5] = 1.0*dwdx1;
    JB[14] = -1.0*dwdx2;
    JB[15] = 1.0*dwdx2;
}