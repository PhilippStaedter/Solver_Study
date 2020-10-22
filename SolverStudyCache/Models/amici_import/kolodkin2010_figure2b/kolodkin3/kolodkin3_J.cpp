#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_kolodkin3(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[10] = 1.0*dwdx0;
    J[11] = 1.0*dwdx2 - 1.0*dwdx3 - 3.4444444444400002*dwdx4;
    J[13] = -3.4444444444400002*dwdx7;
    J[14] = -1.0*dwdx10;
    J[16] = -3.4444444444400002*dwdx13;
    J[17] = -1.0*dwdx15;
    J[20] = 1.0*dwdx1;
    J[22] = -1.0*dwdx5 + 1.0*dwdx6;
    J[23] = -1.0*dwdx8;
    J[25] = 1.0*dwdx12;
    J[31] = 1.0*dwdx4;
    J[32] = 1.0*dwdx5;
    J[33] = 1.0*dwdx7 + 1.0*dwdx8;
    J[36] = 1.0*dwdx13;
    J[41] = 1.0*dwdx3;
    J[44] = 1.0*dwdx10 - 1.0*dwdx9;
    J[47] = 1.0*dwdx15;
    J[48] = -1.0*dwdx16;
    J[49] = -1.0*dwdx17;
    J[50] = -1.0*dwdx1;
    J[52] = -1.0*dwdx6;
    J[55] = -1.0*dwdx11 - 1.0*dwdx12;
    J[56] = -1.0*dwdx14;
    J[61] = -1.0*dwdx4;
    J[63] = -1.0*dwdx7;
    J[65] = 1.0*dwdx11;
    J[66] = -1.0*dwdx13 + 1.0*dwdx14;
    J[71] = -1.0*dwdx3;
    J[74] = -1.0*dwdx10;
    J[77] = -1.0*dwdx15;
    J[84] = -1.0*dwdx9;
    J[88] = -1.0*dwdx16;
    J[89] = -1.0*dwdx17;
    J[94] = 1.0*dwdx9;
    J[98] = 1.0*dwdx16;
    J[99] = 1.0*dwdx17;
}