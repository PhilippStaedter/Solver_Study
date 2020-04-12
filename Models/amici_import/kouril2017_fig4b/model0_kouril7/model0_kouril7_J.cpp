#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_kouril7(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1;
    J[1] = 1.0*dwdx2;
    J[2] = 1.0*dwdx3;
    J[7] = 1.0*dwdx9;
    J[9] = -1.0*dwdx11;
    J[11] = -1.0*dwdx0 + 1.0*dwdx1;
    J[12] = -1.0*dwdx2;
    J[13] = -1.0*dwdx3;
    J[18] = -1.0*dwdx9;
    J[20] = 1.0*dwdx11;
    J[22] = 1.0*dwdx0;
    J[23] = 1.0*dwdx2;
    J[24] = 1.0*dwdx3 - 1.0*dwdx4;
    J[27] = -1.0*dwdx6;
    J[28] = -1.0*dwdx8;
    J[29] = 1.0*dwdx9;
    J[30] = -1.0*dwdx10;
    J[48] = -1.0*dwdx5;
    J[49] = -1.0*dwdx7;
    J[57] = 1.0*dwdx4;
    J[59] = -1.0*dwdx5;
    J[60] = 1.0*dwdx6 - 1.0*dwdx7;
    J[61] = 1.0*dwdx8;
    J[63] = 1.0*dwdx10;
    J[68] = -1.0*dwdx4;
    J[70] = 1.0*dwdx5;
    J[71] = -1.0*dwdx6 + 1.0*dwdx7;
    J[72] = -1.0*dwdx8;
    J[74] = -1.0*dwdx10;
    J[77] = -1.0*dwdx0;
    J[78] = -1.0*dwdx2;
    J[79] = -1.0*dwdx3;
    J[84] = -1.0*dwdx9;
    J[90] = 1.0*dwdx4;
    J[93] = 1.0*dwdx6;
    J[94] = 1.0*dwdx8;
    J[96] = 1.0*dwdx10;
    J[99] = -1.0*dwdx1;
    J[108] = -1.0*dwdx11;
    J[110] = 1.0*dwdx1;
    J[119] = 1.0*dwdx11;
}