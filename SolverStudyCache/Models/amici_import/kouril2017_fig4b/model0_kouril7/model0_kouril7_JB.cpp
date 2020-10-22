#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_kouril7(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = -1.0*dwdx0 + 1.0*dwdx1;
    JB[1] = 1.0*dwdx0 - 1.0*dwdx1;
    JB[2] = -1.0*dwdx0;
    JB[7] = 1.0*dwdx0;
    JB[9] = 1.0*dwdx1;
    JB[10] = -1.0*dwdx1;
    JB[11] = -1.0*dwdx2;
    JB[12] = 1.0*dwdx2;
    JB[13] = -1.0*dwdx2;
    JB[18] = 1.0*dwdx2;
    JB[22] = -1.0*dwdx3;
    JB[23] = 1.0*dwdx3;
    JB[24] = -1.0*dwdx3 + 1.0*dwdx4;
    JB[27] = -1.0*dwdx4;
    JB[28] = 1.0*dwdx4;
    JB[29] = 1.0*dwdx3;
    JB[30] = -1.0*dwdx4;
    JB[48] = 1.0*dwdx5;
    JB[49] = 1.0*dwdx5;
    JB[50] = -1.0*dwdx5;
    JB[57] = 1.0*dwdx6;
    JB[59] = 1.0*dwdx7;
    JB[60] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[61] = 1.0*dwdx6 - 1.0*dwdx7;
    JB[63] = -1.0*dwdx6;
    JB[68] = 1.0*dwdx8;
    JB[71] = -1.0*dwdx8;
    JB[72] = 1.0*dwdx8;
    JB[74] = -1.0*dwdx8;
    JB[77] = -1.0*dwdx9;
    JB[78] = 1.0*dwdx9;
    JB[79] = -1.0*dwdx9;
    JB[84] = 1.0*dwdx9;
    JB[90] = 1.0*dwdx10;
    JB[93] = -1.0*dwdx10;
    JB[94] = 1.0*dwdx10;
    JB[96] = -1.0*dwdx10;
    JB[99] = 1.0*dwdx11;
    JB[100] = -1.0*dwdx11;
    JB[108] = 1.0*dwdx11;
    JB[109] = -1.0*dwdx11;
}