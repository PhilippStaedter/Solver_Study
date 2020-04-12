#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JB_model0_essunger1(realtype *JB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JB[0] = 1.0*dwdx0 + 1.0*dwdx1 + 1.0*dwdx4 - 1.0*dwdx5;
    JB[1] = -1.0*dwdx0 + 1.0*dwdx2;
    JB[2] = 1.0*dwdx3 - 1.0*dwdx6;
    JB[3] = -1.0*dwdx2;
    JB[4] = 1.0*dwdx1 + 1.0*dwdx3 + 1.0*dwdx4 - 1.0*dwdx5 - 1.0*dwdx6;
    JB[6] = 1.0*dwdx0;
    JB[8] = -1.0*dwdx10 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[10] = -1.0*dwdx11 + 1.0*dwdx12;
    JB[11] = -1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx12 + 1.0*dwdx8 + 1.0*dwdx9;
    JB[13] = -1.0*dwdx7;
    JB[14] = -1.0*dwdx13;
    JB[16] = 1.0*dwdx13 + 1.0*dwdx14;
    JB[18] = 1.0*dwdx14;
    JB[22] = -1.0*dwdx15;
    JB[24] = 1.0*dwdx15 + 1.0*dwdx16;
    JB[25] = 1.0*dwdx16;
    JB[28] = -1.0*dwdx20;
    JB[29] = -1.0*dwdx18;
    JB[30] = 1.0*dwdx17;
    JB[31] = 1.0*dwdx19;
    JB[32] = 1.0*dwdx17 - 1.0*dwdx18 + 1.0*dwdx19 - 1.0*dwdx20;
    JB[35] = -1.0*dwdx21;
    JB[39] = 1.0*dwdx22;
    JB[40] = 1.0*dwdx21 + 1.0*dwdx22;
    JB[42] = 1.0*dwdx23;
    JB[43] = -1.0*dwdx23 + 1.0*dwdx25;
    JB[45] = -1.0*dwdx25;
    JB[48] = 1.0*dwdx23 + 1.0*dwdx24;
}