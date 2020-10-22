#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model1_valero(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[2] = -1.0*dwdx0;
    JB[5] = -1.0*dwdx0;
    JB[6] = -2.0*dwdx1;
    JB[7] = 1.0*dwdx1;
    JB[8] = 1.0*dwdx1;
    JB[12] = -2.0*dwdx3;
    JB[13] = -1.0*dwdx2 + 1.0*dwdx3;
    JB[14] = 1.0*dwdx2 + 1.0*dwdx3;
    JB[34] = 1.0*dwdx4;
    JB[35] = 1.0*dwdx4;
}