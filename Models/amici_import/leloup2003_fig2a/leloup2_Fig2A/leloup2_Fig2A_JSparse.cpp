#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_leloup2_Fig2A(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx1 - 1.0*dwdx2;
    JSparse[1] = 1.0*dwdx0;
    JSparse[2] = -1.0*dwdx3 - 1.0*dwdx4 - 1.0*dwdx5;
    JSparse[3] = 1.0*dwdx3;
    JSparse[4] = 1.0*dwdx5;
    JSparse[5] = 1.0*dwdx6 - 1.0*dwdx7 - 1.0*dwdx8 - 1.0*dwdx9;
    JSparse[6] = 1.0*dwdx7;
    JSparse[7] = 1.0*dwdx8;
    JSparse[8] = -1.0*dwdx6;
    JSparse[9] = -1.0*dwdx6;
    JSparse[10] = 1.0*dwdx12;
    JSparse[11] = -1.0*dwdx10 - 1.0*dwdx11 - 1.0*dwdx12;
    JSparse[12] = -1.0*dwdx13;
    JSparse[13] = 1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15 - 1.0*dwdx16;
    JSparse[14] = 1.0*dwdx14;
    JSparse[15] = 1.0*dwdx15;
    JSparse[16] = -1.0*dwdx15;
    JSparse[17] = -1.0*dwdx17 - 1.0*dwdx18 - 1.0*dwdx19;
    JSparse[18] = 1.0*dwdx19;
    JSparse[19] = 1.0*dwdx22;
    JSparse[20] = -1.0*dwdx20 - 1.0*dwdx21 - 1.0*dwdx22;
    JSparse[21] = -1.0*dwdx23;
    JSparse[22] = 1.0*dwdx23 - 1.0*dwdx24 - 1.0*dwdx25;
    JSparse[23] = -1.0*dwdx23;
    JSparse[24] = 1.0*dwdx28;
    JSparse[25] = -1.0*dwdx26 - 1.0*dwdx27 - 1.0*dwdx28;
    JSparse[26] = 1.0*dwdx29;
    JSparse[27] = -1.0*dwdx33;
    JSparse[28] = -1.0*dwdx31;
    JSparse[29] = 1.0*dwdx32;
    JSparse[30] = 1.0*dwdx31;
    JSparse[31] = -1.0*dwdx31 - 1.0*dwdx32 + 1.0*dwdx33 - 1.0*dwdx34;
    JSparse[32] = 1.0*dwdx30;
    JSparse[33] = 1.0*dwdx35;
    JSparse[34] = 1.0*dwdx36;
    JSparse[35] = -1.0*dwdx36 - 1.0*dwdx37 - 1.0*dwdx38;
    JSparse[36] = 1.0*dwdx38;
    JSparse[37] = -1.0*dwdx36;
    JSparse[38] = 1.0*dwdx40;
    JSparse[39] = -1.0*dwdx39 - 1.0*dwdx41;
    JSparse[40] = 1.0*dwdx43;
    JSparse[41] = -1.0*dwdx42 - 1.0*dwdx43 - 1.0*dwdx44;
    JSparse[42] = -1.0*dwdx46 - 1.0*dwdx47;
    JSparse[43] = 1.0*dwdx45;
    JSparse[44] = 1.0*dwdx49;
    JSparse[45] = -1.0*dwdx49;
    JSparse[46] = -1.0*dwdx48 - 1.0*dwdx49 - 1.0*dwdx50;
    JSparse[47] = 1.0*dwdx48;
    JSparse[48] = 1.0*dwdx53;
    JSparse[49] = -1.0*dwdx51 - 1.0*dwdx52 - 1.0*dwdx53;
}