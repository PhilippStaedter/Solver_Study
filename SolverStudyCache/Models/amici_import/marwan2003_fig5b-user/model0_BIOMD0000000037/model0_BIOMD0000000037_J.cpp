#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_BIOMD0000000037(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[13] = -1.0*dwdx1;
    J[15] = 1.0*dwdx3;
    J[27] = 1.0*dwdx2;
    J[37] = 1.0*dwdx1;
    J[39] = -1.0*dwdx2 - 1.0*dwdx3;
    J[52] = -1.0*dwdx5;
    J[53] = 1.0*dwdx7;
    J[56] = 1.0*dwdx13;
    J[58] = 1.0*dwdx14;
    J[64] = 1.0*dwdx6;
    J[65] = -1.0*dwdx8;
    J[75] = 1.0*dwdx4;
    J[78] = -1.0*dwdx9;
    J[79] = 1.0*dwdx11;
    J[87] = -1.0*dwdx4;
    J[90] = 1.0*dwdx9;
    J[91] = -1.0*dwdx11;
    J[96] = -1.0*dwdx0;
    J[104] = -1.0*dwdx12;
    J[108] = 1.0*dwdx0;
    J[116] = 1.0*dwdx12;
    J[126] = 1.0*dwdx10;
    J[128] = -1.0*dwdx13;
    J[130] = -1.0*dwdx14;
    J[131] = 1.0*dwdx15;
    J[138] = -1.0*dwdx10;
    J[143] = -1.0*dwdx15;
}