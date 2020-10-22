#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_kouril8(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = 1.0*dwdx0 - 1.0*dwdx1;
    J[1] = 1.0*dwdx2;
    J[2] = 1.0*dwdx3;
    J[4] = 1.0*dwdx6;
    J[9] = -1.0*dwdx13;
    J[13] = -1.0*dwdx0 + 1.0*dwdx1;
    J[14] = -1.0*dwdx2;
    J[15] = -1.0*dwdx3;
    J[17] = -1.0*dwdx6;
    J[22] = 1.0*dwdx13;
    J[26] = 1.0*dwdx0;
    J[27] = 1.0*dwdx2;
    J[28] = 1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    J[30] = 1.0*dwdx6;
    J[31] = -1.0*dwdx7;
    J[33] = -1.0*dwdx10;
    J[34] = -1.0*dwdx12;
    J[36] = -1.0*dwdx14;
    J[52] = -1.0*dwdx0;
    J[53] = -1.0*dwdx2;
    J[54] = -1.0*dwdx3 + 1.0*dwdx4;
    J[56] = -1.0*dwdx6;
    J[67] = 1.0*dwdx5;
    J[70] = 1.0*dwdx7 - 1.0*dwdx8;
    J[72] = 1.0*dwdx10;
    J[73] = 1.0*dwdx12;
    J[75] = 1.0*dwdx14;
    J[84] = -1.0*dwdx9;
    J[85] = -1.0*dwdx11;
    J[93] = 1.0*dwdx5;
    J[96] = 1.0*dwdx7;
    J[97] = -1.0*dwdx9;
    J[98] = 1.0*dwdx10 - 1.0*dwdx11;
    J[99] = 1.0*dwdx12;
    J[101] = 1.0*dwdx14;
    J[106] = -1.0*dwdx5;
    J[109] = -1.0*dwdx7;
    J[110] = 1.0*dwdx9;
    J[111] = -1.0*dwdx10 + 1.0*dwdx11;
    J[112] = -1.0*dwdx12;
    J[114] = -1.0*dwdx14;
    J[117] = -1.0*dwdx1;
    J[126] = -1.0*dwdx13;
    J[143] = 1.0*dwdx1;
    J[152] = 1.0*dwdx13;
}