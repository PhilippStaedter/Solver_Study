#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void J_model2_sarma7(realtype *J, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    J[0] = -1.0*dwdx0;
    J[1] = 1.0*dwdx3;
    J[10] = -1.0*dwdx24;
    J[11] = 1.0*dwdx30;
    J[13] = 1.0*dwdx0;
    J[14] = -1.0*dwdx3;
    J[23] = 1.0*dwdx24;
    J[24] = -1.0*dwdx30;
    J[42] = -1.0*dwdx6;
    J[43] = 1.0*dwdx7;
    J[49] = -1.0*dwdx26;
    J[55] = 1.0*dwdx6;
    J[56] = -1.0*dwdx7;
    J[62] = 1.0*dwdx26;
    J[66] = -1.0*dwdx1;
    J[70] = -1.0*dwdx8;
    J[71] = -1.0*dwdx10 + 1.0*dwdx13;
    J[72] = 1.0*dwdx17;
    J[75] = -1.0*dwdx25;
    J[77] = 1.0*dwdx32;
    J[79] = 1.0*dwdx1 - 1.0*dwdx2;
    J[83] = 1.0*dwdx8 - 1.0*dwdx9;
    J[84] = 1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12 - 1.0*dwdx13;
    J[85] = 1.0*dwdx16 - 1.0*dwdx17;
    J[88] = 1.0*dwdx25 - 1.0*dwdx27;
    J[90] = 1.0*dwdx31 - 1.0*dwdx32;
    J[92] = 1.0*dwdx2;
    J[96] = 1.0*dwdx9;
    J[97] = 1.0*dwdx11 - 1.0*dwdx12;
    J[98] = -1.0*dwdx16;
    J[101] = 1.0*dwdx27;
    J[103] = -1.0*dwdx31;
    J[106] = 1.0*dwdx5;
    J[111] = -1.0*dwdx14;
    J[112] = -1.0*dwdx18;
    J[113] = -1.0*dwdx20 + 1.0*dwdx23;
    J[114] = 1.0*dwdx29;
    J[119] = 1.0*dwdx4 - 1.0*dwdx5;
    J[124] = 1.0*dwdx14 - 1.0*dwdx15;
    J[125] = 1.0*dwdx18 - 1.0*dwdx19;
    J[126] = 1.0*dwdx20 - 1.0*dwdx21 + 1.0*dwdx22 - 1.0*dwdx23;
    J[127] = 1.0*dwdx28 - 1.0*dwdx29;
    J[132] = -1.0*dwdx4;
    J[137] = 1.0*dwdx15;
    J[138] = 1.0*dwdx19;
    J[139] = 1.0*dwdx21 - 1.0*dwdx22;
    J[140] = -1.0*dwdx28;
}