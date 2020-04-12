#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_sarma5(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx0;
    JB[11] = -1.0*dwdx3;
    JB[12] = 1.0*dwdx3;
    JB[14] = 1.0*dwdx1;
    JB[15] = -1.0*dwdx1 + 1.0*dwdx2;
    JB[16] = -1.0*dwdx2;
    JB[28] = -1.0*dwdx5;
    JB[29] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[30] = 1.0*dwdx4;
    JB[36] = 1.0*dwdx6;
    JB[37] = -1.0*dwdx6 + 1.0*dwdx7;
    JB[38] = -1.0*dwdx7;
    JB[47] = -1.0*dwdx11 + 1.0*dwdx8;
    JB[48] = -1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx8 + 1.0*dwdx9;
    JB[49] = 1.0*dwdx10 - 1.0*dwdx9;
    JB[58] = -1.0*dwdx15;
    JB[59] = -1.0*dwdx14 + 1.0*dwdx15;
    JB[60] = 1.0*dwdx14;
    JB[61] = 1.0*dwdx12;
    JB[62] = -1.0*dwdx12 + 1.0*dwdx13;
    JB[63] = -1.0*dwdx13;
    JB[72] = 1.0*dwdx16;
    JB[73] = -1.0*dwdx16 + 1.0*dwdx17;
    JB[74] = -1.0*dwdx17;
    JB[83] = 1.0*dwdx18 - 1.0*dwdx21;
    JB[84] = -1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20 + 1.0*dwdx21;
    JB[85] = -1.0*dwdx19 + 1.0*dwdx20;
    JB[88] = 1.0*dwdx22;
    JB[89] = -1.0*dwdx22;
    JB[91] = 1.0*dwdx23;
    JB[92] = -1.0*dwdx23 + 1.0*dwdx24;
    JB[93] = -1.0*dwdx24;
    JB[94] = -1.0*dwdx26;
    JB[95] = -1.0*dwdx25 + 1.0*dwdx26;
    JB[96] = 1.0*dwdx25;
    JB[99] = -1.0*dwdx27;
    JB[100] = 1.0*dwdx27;
    JB[113] = -1.0*dwdx29;
    JB[114] = -1.0*dwdx28 + 1.0*dwdx29;
    JB[115] = 1.0*dwdx28;
}