#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_martins1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[1] = 1.0*dwdx1;
    J[12] = 1.0*dwdx14;
    J[15] = -1.0*dwdx0 - 1.0*dwdx1;
    J[16] = 1.0*dwdx4;
    J[22] = -1.0*dwdx11;
    J[26] = 1.0*dwdx13;
    J[27] = 1.0*dwdx15;
    J[30] = -1.0*dwdx2 - 1.0*dwdx3 - 1.0*dwdx4;
    J[44] = 1.0*dwdx2;
    J[45] = -1.0*dwdx5 - 1.0*dwdx6 - 1.0*dwdx7;
    J[58] = 1.0*dwdx3;
    J[60] = -1.0*dwdx8 - 1.0*dwdx9;
    J[71] = 1.0*dwdx1;
    J[83] = 1.0*dwdx16;
    J[88] = 1.0*dwdx8;
    J[101] = 1.0*dwdx6;
    J[105] = -1.0*dwdx10;
    J[108] = 1.0*dwdx12;
    J[113] = -1.0*dwdx0;
    J[114] = 1.0*dwdx4;
    J[115] = 1.0*dwdx5 + 1.0*dwdx6 + 1.0*dwdx7;
    J[116] = 1.0*dwdx8 + 1.0*dwdx9;
    J[120] = -1.0*dwdx11;
    J[127] = 1.0*dwdx1;
    J[143] = 1.0*dwdx5;
    J[150] = -1.0*dwdx12;
    J[155] = 1.0*dwdx0;
    J[162] = 1.0*dwdx11;
    J[172] = 1.0*dwdx9;
    J[180] = -1.0*dwdx13 - 1.0*dwdx14;
    J[185] = 1.0*dwdx7;
    J[189] = 1.0*dwdx10;
    J[195] = -1.0*dwdx15 - 1.0*dwdx16;
}