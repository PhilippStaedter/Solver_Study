#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_kolodkin2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[6] = 1.0*dwdx0;
    J[7] = 1.0*dwdx1 - 1.0*dwdx2;
    J[8] = -1.0*dwdx4;
    J[9] = -1.0*dwdx5;
    J[13] = 1.0*dwdx2;
    J[14] = -1.0*dwdx3 + 1.0*dwdx4;
    J[15] = 1.0*dwdx5;
    J[16] = -1.0*dwdx6;
    J[17] = -1.0*dwdx7;
    J[19] = -1.0*dwdx2;
    J[20] = -1.0*dwdx4;
    J[21] = -1.0*dwdx5;
    J[26] = -1.0*dwdx3;
    J[28] = -1.0*dwdx6;
    J[29] = -1.0*dwdx7;
    J[32] = 1.0*dwdx3;
    J[34] = 1.0*dwdx6;
    J[35] = 1.0*dwdx7;
}