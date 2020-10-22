#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_Ueda2001(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[14] = -1.0*dwdx0 + 1.0*dwdx1 - 1.0*dwdx2 - 1.0*dwdx3;
    J[15] = 1.0*dwdx7;
    J[16] = 1.0*dwdx10;
    J[24] = 1.0*dwdx38;
    J[27] = 1.0*dwdx0;
    J[28] = -1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    J[40] = -1.0*dwdx1;
    J[42] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    J[43] = 1.0*dwdx14;
    J[50] = -1.0*dwdx38;
    J[54] = 1.0*dwdx6;
    J[56] = -1.0*dwdx13 - 1.0*dwdx15;
    J[60] = 1.0*dwdx28;
    J[70] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    J[71] = 1.0*dwdx20;
    J[72] = -1.0*dwdx23;
    J[74] = -1.0*dwdx32;
    J[77] = -1.0*dwdx39;
    J[80] = 1.0*dwdx4;
    J[84] = -1.0*dwdx19 - 1.0*dwdx21;
    J[86] = 1.0*dwdx26;
    J[96] = 1.0*dwdx16;
    J[98] = -1.0*dwdx22 + 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25;
    J[99] = 1.0*dwdx29;
    J[100] = 1.0*dwdx32;
    J[111] = 1.0*dwdx22;
    J[112] = -1.0*dwdx29 - 1.0*dwdx30 - 1.0*dwdx31;
    J[122] = -1.0*dwdx16;
    J[124] = -1.0*dwdx23;
    J[126] = -1.0*dwdx32 - 1.0*dwdx33 - 1.0*dwdx34;
    J[127] = 1.0*dwdx36;
    J[132] = 1.0*dwdx5;
    J[138] = 1.0*dwdx27;
    J[140] = -1.0*dwdx35 - 1.0*dwdx37;
}