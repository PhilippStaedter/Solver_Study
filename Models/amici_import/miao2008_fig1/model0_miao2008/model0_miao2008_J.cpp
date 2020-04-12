#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_miao2008(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    J[1] = -1.0*dwdx4;
    J[2] = -1.0*dwdx8;
    J[3] = -1.0*dwdx10;
    J[4] = 1.0*dwdx1 + 0.25*dwdx3;
    J[5] = 1.0*dwdx4 + 1.0*dwdx5 - 1.0*dwdx6;
    J[6] = 0.25*dwdx8;
    J[7] = -1.0*dwdx11;
    J[8] = 0.5*dwdx3;
    J[9] = 1.0*dwdx6 + 1.0*dwdx7;
    J[10] = 0.5*dwdx8 + 1.0*dwdx9;
    J[11] = 1.0*dwdx11 + 1.0*dwdx13;
    J[12] = 1.0*dwdx2 + 0.25*dwdx3;
    J[13] = -1.0*dwdx7;
    J[14] = 0.25*dwdx8;
    J[15] = 1.0*dwdx10 + 1.0*dwdx12 - 1.0*dwdx13;
}