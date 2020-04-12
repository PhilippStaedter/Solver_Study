#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_panteleev2(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1;
    J[3] = -1.0*dwdx9;
    J[4] = -1.0*dwdx10;
    J[6] = -1.0*dwdx15;
    J[7] = -1.0*dwdx16;
    J[10] = -1.0*dwdx2 + 1.0*dwdx3 - 1.0*dwdx4 + 1.0*dwdx5;
    J[11] = -1.0*dwdx6;
    J[12] = 1.0*dwdx8;
    J[13] = 1.0*dwdx11;
    J[14] = -1.0*dwdx13;
    J[15] = 1.0*dwdx14;
    J[16] = -1.0*dwdx17 + 1.0*dwdx18;
    J[17] = -1.0*dwdx19;
    J[19] = 1.0*dwdx2;
    J[20] = 1.0*dwdx6 - 1.0*dwdx7;
    J[23] = 1.0*dwdx13;
    J[27] = -1.0*dwdx1;
    J[28] = -1.0*dwdx3;
    J[29] = 1.0*dwdx7;
    J[30] = -1.0*dwdx8 - 1.0*dwdx9;
    J[31] = -1.0*dwdx10;
    J[33] = -1.0*dwdx14;
    J[36] = 1.0*dwdx1;
    J[37] = -1.0*dwdx5;
    J[39] = 1.0*dwdx9;
    J[40] = 1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    J[43] = -1.0*dwdx18;
    J[44] = -1.0*dwdx20;
    J[46] = -1.0*dwdx2;
    J[47] = -1.0*dwdx6;
    J[50] = -1.0*dwdx13;
    J[54] = -1.0*dwdx0;
    J[55] = 1.0*dwdx3;
    J[57] = 1.0*dwdx8;
    J[60] = 1.0*dwdx14 - 1.0*dwdx15;
    J[61] = -1.0*dwdx16;
    J[63] = 1.0*dwdx0;
    J[64] = -1.0*dwdx4 + 1.0*dwdx5;
    J[67] = 1.0*dwdx11;
    J[69] = 1.0*dwdx15;
    J[70] = 1.0*dwdx16 - 1.0*dwdx17 + 1.0*dwdx18;
    J[71] = -1.0*dwdx19;
    J[73] = 1.0*dwdx4;
    J[76] = 1.0*dwdx12;
    J[79] = 1.0*dwdx17;
    J[80] = 1.0*dwdx19 + 1.0*dwdx20;
}