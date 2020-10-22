#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_marhl(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -0.25*dwdx0 - 0.25*dwdx1;
    J[3] = -0.25*dwdx4 - 0.25*dwdx6 + 0.25*dwdx7;
    J[6] = -0.25*dwdx2;
    J[8] = -0.25*dwdx8 + 0.25*dwdx9;
    J[12] = -1.0*dwdx3;
    J[13] = 1.0*dwdx5;
    J[14] = 1.0*dwdx10;
    J[15] = 1.0*dwdx0 + 1.0*dwdx1;
    J[16] = 1.0*dwdx2;
    J[17] = 1.0*dwdx3;
    J[18] = 1.0*dwdx4 - 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7 + 1.0*dwdx8 - 1.0*dwdx9;
    J[19] = -1.0*dwdx10;
    J[22] = 1.0*dwdx3;
    J[23] = -1.0*dwdx5;
    J[24] = -1.0*dwdx10;
}