#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_aguda1(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    J[1] = 1.0*dwdx5;
    J[5] = -1.0*dwdx17;
    J[8] = -1.0*dwdx25;
    J[12] = 1.0*dwdx3;
    J[13] = -1.0*dwdx5;
    J[20] = 1.0*dwdx25;
    J[26] = -1.0*dwdx6;
    J[28] = 1.0*dwdx13;
    J[32] = 1.0*dwdx24;
    J[36] = 1.0*dwdx0;
    J[37] = 1.0*dwdx4;
    J[39] = -1.0*dwdx10 + 1.0*dwdx7 - 1.0*dwdx8;
    J[40] = 1.0*dwdx11;
    J[45] = -1.0*dwdx27;
    J[46] = 1.0*dwdx30;
    J[50] = 1.0*dwdx6;
    J[52] = -1.0*dwdx13 - 1.0*dwdx14 + 1.0*dwdx15 - 1.0*dwdx16;
    J[55] = 1.0*dwdx21;
    J[56] = -1.0*dwdx24;
    J[60] = -1.0*dwdx2;
    J[65] = -1.0*dwdx17 - 1.0*dwdx18;
    J[69] = 1.0*dwdx28;
    J[72] = 1.0*dwdx0;
    J[73] = 1.0*dwdx4;
    J[76] = 1.0*dwdx11;
    J[78] = -1.0*dwdx20;
    J[82] = 1.0*dwdx30;
    J[87] = 1.0*dwdx9;
    J[88] = -1.0*dwdx15 + 1.0*dwdx16;
    J[91] = -1.0*dwdx21 - 1.0*dwdx22;
    J[96] = -1.0*dwdx3;
    J[97] = 1.0*dwdx5;
    J[98] = 1.0*dwdx6;
    J[100] = -1.0*dwdx12 - 1.0*dwdx13;
    J[104] = -1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26;
    J[111] = -1.0*dwdx8;
    J[113] = 1.0*dwdx19;
    J[114] = 1.0*dwdx20;
    J[117] = -1.0*dwdx27 - 1.0*dwdx29;
    J[120] = -1.0*dwdx0;
    J[121] = -1.0*dwdx4;
    J[123] = 1.0*dwdx8;
    J[124] = -1.0*dwdx11;
    J[129] = 1.0*dwdx27;
    J[130] = -1.0*dwdx30;
}