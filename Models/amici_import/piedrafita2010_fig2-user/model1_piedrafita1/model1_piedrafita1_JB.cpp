#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model1_piedrafita1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[2] = 1.0*dwdx0;
    JB[3] = -1.0*dwdx0;
    JB[12] = 1.0*dwdx1 - 1.0*dwdx2 + 1.0*dwdx3;
    JB[13] = -1.0*dwdx2;
    JB[15] = 1.0*dwdx2;
    JB[17] = 1.0*dwdx3;
    JB[18] = -1.0*dwdx3;
    JB[23] = -1.0*dwdx6;
    JB[24] = 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6 + 1.0*dwdx7 - 1.0*dwdx8;
    JB[25] = -1.0*dwdx4;
    JB[26] = 1.0*dwdx6;
    JB[27] = 1.0*dwdx5;
    JB[28] = -1.0*dwdx5 - 1.0*dwdx8;
    JB[30] = 1.0*dwdx8;
    JB[35] = 1.0*dwdx9;
    JB[36] = 1.0*dwdx10 + 1.0*dwdx11 - 1.0*dwdx9;
    JB[37] = -1.0*dwdx10;
    JB[38] = -1.0*dwdx11;
    JB[45] = -1.0*dwdx13;
    JB[46] = -1.0*dwdx13;
    JB[47] = 1.0*dwdx12;
    JB[48] = -1.0*dwdx12 + 1.0*dwdx13;
    JB[57] = -1.0*dwdx14;
    JB[58] = 1.0*dwdx15;
    JB[60] = 1.0*dwdx14 - 1.0*dwdx15;
    JB[61] = -1.0*dwdx14;
    JB[67] = 1.0*dwdx17;
    JB[68] = -1.0*dwdx16 - 1.0*dwdx18;
    JB[71] = 1.0*dwdx16;
    JB[72] = -1.0*dwdx16 + 1.0*dwdx17 - 1.0*dwdx18 + 1.0*dwdx19;
    JB[73] = -1.0*dwdx17;
    JB[74] = 1.0*dwdx18;
    JB[78] = 1.0*dwdx20;
    JB[83] = 1.0*dwdx20;
    JB[84] = -1.0*dwdx20 + 1.0*dwdx21;
    JB[85] = -1.0*dwdx21;
    JB[90] = -1.0*dwdx23;
    JB[94] = -1.0*dwdx23;
    JB[95] = 1.0*dwdx22;
    JB[96] = -1.0*dwdx22 + 1.0*dwdx23;
    JB[102] = 1.0*dwdx24;
    JB[103] = -1.0*dwdx24;
    JB[113] = 1.0*dwdx26;
    JB[115] = -1.0*dwdx26;
    JB[117] = 1.0*dwdx25;
    JB[118] = -1.0*dwdx25;
}