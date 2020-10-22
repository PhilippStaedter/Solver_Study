#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_piedrafita1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[12] = -1.0*dwdx1 + 1.0*dwdx2 - 1.0*dwdx3;
    J[13] = 1.0*dwdx6;
    J[15] = 1.0*dwdx13;
    J[17] = -1.0*dwdx17;
    J[18] = -1.0*dwdx20;
    J[22] = -1.0*dwdx0;
    J[23] = 1.0*dwdx2;
    J[24] = -1.0*dwdx4 + 1.0*dwdx5 + 1.0*dwdx6 - 1.0*dwdx7 + 1.0*dwdx8;
    J[25] = -1.0*dwdx9;
    J[26] = 1.0*dwdx13;
    J[27] = 1.0*dwdx14;
    J[28] = 1.0*dwdx16 + 1.0*dwdx18;
    J[30] = 1.0*dwdx23;
    J[33] = 1.0*dwdx0;
    J[35] = 1.0*dwdx4;
    J[36] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx9;
    J[37] = -1.0*dwdx12;
    J[38] = -1.0*dwdx15;
    J[42] = -1.0*dwdx24;
    J[43] = -1.0*dwdx26;
    J[45] = -1.0*dwdx2;
    J[46] = -1.0*dwdx6;
    J[47] = 1.0*dwdx10;
    J[48] = 1.0*dwdx12 - 1.0*dwdx13;
    J[53] = 1.0*dwdx24;
    J[57] = -1.0*dwdx5;
    J[58] = 1.0*dwdx11;
    J[60] = -1.0*dwdx14 + 1.0*dwdx15;
    J[61] = -1.0*dwdx16;
    J[65] = 1.0*dwdx26;
    J[67] = -1.0*dwdx3;
    J[68] = 1.0*dwdx5 + 1.0*dwdx8;
    J[71] = 1.0*dwdx14;
    J[72] = 1.0*dwdx16 - 1.0*dwdx17 + 1.0*dwdx18 - 1.0*dwdx19;
    J[73] = -1.0*dwdx20;
    J[74] = 1.0*dwdx23;
    J[78] = 1.0*dwdx3;
    J[83] = 1.0*dwdx17;
    J[84] = 1.0*dwdx20 - 1.0*dwdx21;
    J[85] = -1.0*dwdx22;
    J[87] = -1.0*dwdx25;
    J[90] = -1.0*dwdx8;
    J[94] = -1.0*dwdx18;
    J[95] = 1.0*dwdx21;
    J[96] = 1.0*dwdx22 - 1.0*dwdx23;
    J[98] = 1.0*dwdx25;
}