#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model1_sarma7(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0;
    JB[1] = -1.0*dwdx0;
    JB[13] = -1.0*dwdx3;
    JB[14] = 1.0*dwdx3;
    JB[18] = 1.0*dwdx1;
    JB[19] = -1.0*dwdx1 + 1.0*dwdx2;
    JB[20] = -1.0*dwdx2;
    JB[34] = -1.0*dwdx5;
    JB[35] = -1.0*dwdx4 + 1.0*dwdx5;
    JB[36] = 1.0*dwdx4;
    JB[42] = 1.0*dwdx6;
    JB[43] = -1.0*dwdx6;
    JB[55] = -1.0*dwdx7;
    JB[56] = 1.0*dwdx7;
    JB[70] = 1.0*dwdx8;
    JB[71] = -1.0*dwdx8 + 1.0*dwdx9;
    JB[72] = -1.0*dwdx9;
    JB[83] = 1.0*dwdx10 - 1.0*dwdx13;
    JB[84] = -1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx12 + 1.0*dwdx13;
    JB[85] = -1.0*dwdx11 + 1.0*dwdx12;
    JB[96] = -1.0*dwdx17;
    JB[97] = -1.0*dwdx16 + 1.0*dwdx17;
    JB[98] = 1.0*dwdx16;
    JB[99] = 1.0*dwdx14;
    JB[100] = -1.0*dwdx14 + 1.0*dwdx15;
    JB[101] = -1.0*dwdx15;
    JB[112] = 1.0*dwdx18;
    JB[113] = -1.0*dwdx18 + 1.0*dwdx19;
    JB[114] = -1.0*dwdx19;
    JB[125] = 1.0*dwdx20 - 1.0*dwdx23;
    JB[126] = -1.0*dwdx20 + 1.0*dwdx21 - 1.0*dwdx22 + 1.0*dwdx23;
    JB[127] = -1.0*dwdx21 + 1.0*dwdx22;
    JB[130] = 1.0*dwdx24;
    JB[131] = -1.0*dwdx24;
    JB[133] = 1.0*dwdx26;
    JB[134] = -1.0*dwdx26;
    JB[135] = 1.0*dwdx25;
    JB[136] = -1.0*dwdx25 + 1.0*dwdx27;
    JB[137] = -1.0*dwdx27;
    JB[138] = -1.0*dwdx29;
    JB[139] = -1.0*dwdx28 + 1.0*dwdx29;
    JB[140] = 1.0*dwdx28;
    JB[143] = -1.0*dwdx30;
    JB[144] = 1.0*dwdx30;
    JB[161] = -1.0*dwdx32;
    JB[162] = -1.0*dwdx31 + 1.0*dwdx32;
    JB[163] = 1.0*dwdx31;
}