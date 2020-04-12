#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model0_essunger3(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0 - 1.0*dwdx1 - 1.0*dwdx3;
    J[2] = 1.0*dwdx9;
    J[4] = -1.0*dwdx13;
    J[5] = 1.0*dwdx17;
    J[6] = -1.0*dwdx19;
    J[7] = 1.0*dwdx0 - 1.0*dwdx2;
    J[8] = -1.0*dwdx6 - 1.0*dwdx8;
    J[10] = 1.0*dwdx11;
    J[11] = -1.0*dwdx16;
    J[13] = 1.0*dwdx19 - 1.0*dwdx21;
    J[14] = 1.0*dwdx4;
    J[16] = -1.0*dwdx10 - 1.0*dwdx9;
    J[18] = 1.0*dwdx14;
    J[21] = 1.0*dwdx2;
    J[22] = 1.0*dwdx7;
    J[24] = -1.0*dwdx11 - 1.0*dwdx12;
    J[25] = 1.0*dwdx15;
    J[27] = 1.0*dwdx21;
    J[28] = -1.0*dwdx1 + 1.0*dwdx3;
    J[29] = -1.0*dwdx6 + 1.0*dwdx8;
    J[30] = -1.0*dwdx10;
    J[31] = -1.0*dwdx12;
    J[32] = 1.0*dwdx13 + 1.0*dwdx16;
    J[33] = -1.0*dwdx18;
    J[34] = 1.0*dwdx22;
    J[40] = -1.0*dwdx17 - 1.0*dwdx18;
    J[41] = 1.0*dwdx22;
    J[42] = -1.0*dwdx0;
    J[43] = 1.0*dwdx5;
    J[48] = -1.0*dwdx19 - 1.0*dwdx20;
}