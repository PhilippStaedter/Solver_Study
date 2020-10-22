#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_sarma6(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx3;
    J[8] = -1.0*dwdx22;
    J[9] = 1.0*dwdx27;
    J[11] = 1.0*dwdx0;
    J[12] = -1.0*dwdx3;
    J[19] = 1.0*dwdx22;
    J[20] = -1.0*dwdx27;
    J[34] = -1.0*dwdx1;
    J[36] = -1.0*dwdx6;
    J[37] = 1.0*dwdx11 - 1.0*dwdx8;
    J[38] = 1.0*dwdx15;
    J[41] = -1.0*dwdx23;
    J[43] = 1.0*dwdx29;
    J[45] = 1.0*dwdx1 - 1.0*dwdx2;
    J[47] = 1.0*dwdx6 - 1.0*dwdx7;
    J[48] = 1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx8 - 1.0*dwdx9;
    J[49] = 1.0*dwdx14 - 1.0*dwdx15;
    J[52] = 1.0*dwdx23 - 1.0*dwdx24;
    J[54] = 1.0*dwdx28 - 1.0*dwdx29;
    J[56] = 1.0*dwdx2;
    J[58] = 1.0*dwdx7;
    J[59] = -1.0*dwdx10 + 1.0*dwdx9;
    J[60] = -1.0*dwdx14;
    J[63] = 1.0*dwdx24;
    J[65] = -1.0*dwdx28;
    J[68] = 1.0*dwdx5;
    J[71] = -1.0*dwdx12;
    J[72] = -1.0*dwdx16;
    J[73] = -1.0*dwdx18 + 1.0*dwdx21;
    J[74] = 1.0*dwdx26;
    J[79] = 1.0*dwdx4 - 1.0*dwdx5;
    J[82] = 1.0*dwdx12 - 1.0*dwdx13;
    J[83] = 1.0*dwdx16 - 1.0*dwdx17;
    J[84] = 1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20 - 1.0*dwdx21;
    J[85] = 1.0*dwdx25 - 1.0*dwdx26;
    J[90] = -1.0*dwdx4;
    J[93] = 1.0*dwdx13;
    J[94] = 1.0*dwdx17;
    J[95] = 1.0*dwdx19 - 1.0*dwdx20;
    J[96] = -1.0*dwdx25;
}