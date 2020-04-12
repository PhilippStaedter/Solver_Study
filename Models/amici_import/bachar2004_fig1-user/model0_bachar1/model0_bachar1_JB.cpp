#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_bachar1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[5] = -1.0*dwdx1 + 1.0*dwdx2 + 1.0*dwdx3;
    JB[6] = -1.0*dwdx3;
    JB[7] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[9] = -1.0*dwdx5 - 1.0*dwdx6;
    JB[10] = 1.0*dwdx6 + 1.0*dwdx7 + 1.0*dwdx8;
    JB[11] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[13] = -1.0*dwdx10;
    JB[15] = 1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx9;
}