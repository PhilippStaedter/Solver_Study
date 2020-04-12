#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_kolodkin2(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[1] = -1.0*dwdx0;
    JB[7] = -1.0*dwdx1 + 1.0*dwdx2;
    JB[8] = -1.0*dwdx2;
    JB[9] = 1.0*dwdx2;
    JB[13] = 1.0*dwdx4;
    JB[14] = 1.0*dwdx3 - 1.0*dwdx4;
    JB[15] = 1.0*dwdx4;
    JB[16] = 1.0*dwdx3;
    JB[17] = -1.0*dwdx3;
    JB[19] = 1.0*dwdx5;
    JB[20] = -1.0*dwdx5;
    JB[21] = 1.0*dwdx5;
    JB[26] = 1.0*dwdx6;
    JB[28] = 1.0*dwdx6;
    JB[29] = -1.0*dwdx6;
    JB[32] = 1.0*dwdx7;
    JB[34] = 1.0*dwdx7;
    JB[35] = -1.0*dwdx7;
}