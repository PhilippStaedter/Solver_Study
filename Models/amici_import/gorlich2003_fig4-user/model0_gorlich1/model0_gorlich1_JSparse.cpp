#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_gorlich1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 83333333333.333328*dwdx0;
    JSparse[1] = -83333333333.333328*dwdx0;
    JSparse[2] = -83333333333.333328*dwdx1;
    JSparse[3] = 83333333333.333328*dwdx1;
    JSparse[4] = -83333333333.333328*dwdx2 + 83333333333.333328*dwdx3;
    JSparse[5] = 83333333333.333328*dwdx2;
    JSparse[6] = -83333333333.333328*dwdx3;
    JSparse[7] = -83333333333.333328*dwdx2;
    JSparse[8] = 83333333333.333328*dwdx3;
    JSparse[9] = 83333333333.333328*dwdx4 - 83333333333.333328*dwdx5;
    JSparse[10] = -83333333333.333328*dwdx4;
    JSparse[11] = 83333333333.333328*dwdx5;
    JSparse[12] = -83333333333.333328*dwdx7;
    JSparse[13] = 83333333333.333328*dwdx6;
    JSparse[14] = -83333333333.333328*dwdx6 + 83333333333.333328*dwdx7;
    JSparse[15] = -83333333333.333328*dwdx7;
    JSparse[16] = 83333333333.333328*dwdx9;
    JSparse[17] = -83333333333.333328*dwdx8;
    JSparse[18] = 83333333333.333328*dwdx8 - 83333333333.333328*dwdx9;
    JSparse[19] = 83333333333.333328*dwdx9;
    JSparse[20] = -55555555555.555557*dwdx10;
    JSparse[21] = 55555555555.555557*dwdx10;
    JSparse[22] = -55555555555.555557*dwdx10;
    JSparse[23] = 55555555555.555557*dwdx11;
    JSparse[24] = 55555555555.555557*dwdx11 + 55555555555.555557*dwdx12;
    JSparse[25] = -55555555555.555557*dwdx11;
    JSparse[26] = -55555555555.555557*dwdx12;
    JSparse[27] = 55555555555.555557*dwdx13;
    JSparse[28] = -83333333333.333328*dwdx13;
    JSparse[29] = -83333333333.333328*dwdx15;
    JSparse[30] = 83333333333.333328*dwdx15;
    JSparse[31] = 55555555555.555557*dwdx14;
    JSparse[32] = -83333333333.333328*dwdx14 - 83333333333.333328*dwdx15;
    JSparse[33] = 55555555555.555557*dwdx16 - 55555555555.555557*dwdx17;
    JSparse[34] = 55555555555.555557*dwdx16;
    JSparse[35] = -55555555555.555557*dwdx16 + 55555555555.555557*dwdx17;
    JSparse[36] = -55555555555.555557*dwdx17;
    JSparse[37] = -55555555555.555557*dwdx20;
    JSparse[38] = 55555555555.555557*dwdx19;
    JSparse[39] = 55555555555.555557*dwdx20;
    JSparse[40] = 55555555555.555557*dwdx18 - 55555555555.555557*dwdx19 - 55555555555.555557*dwdx20;
    JSparse[41] = -83333333333.333328*dwdx18;
    JSparse[42] = 83333333333.333328*dwdx22;
    JSparse[43] = -83333333333.333328*dwdx22;
    JSparse[44] = 55555555555.555557*dwdx21;
    JSparse[45] = -83333333333.333328*dwdx21 + 83333333333.333328*dwdx22;
}