#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_fisher1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = -3717472118959.1074*dwdx0 - 3717472118959.1074*dwdx1 - 3717472118959.1074*dwdx2 - 3717472118959.1074*dwdx3;
    JSparseB[1] = -3717472118959.1074*dwdx6;
    JSparseB[2] = -3717472118959.1074*dwdx8;
    JSparseB[3] = -3717472118959.1074*dwdx12;
    JSparseB[4] = -3717472118959.1074*dwdx16;
    JSparseB[5] = -3717472118959.1074*dwdx22;
    JSparseB[6] = -3717472118959.1074*dwdx30;
    JSparseB[7] = -3717472118959.1074*dwdx36;
    JSparseB[8] = 8849557522123.8945*dwdx2;
    JSparseB[9] = -8849557522123.8945*dwdx4 + 8849557522123.8945*dwdx5 + 8849557522123.8945*dwdx6 - 8849557522123.8945*dwdx7;
    JSparseB[10] = -8849557522123.8945*dwdx10;
    JSparseB[11] = -8849557522123.8945*dwdx14;
    JSparseB[12] = 8849557522123.8945*dwdx19;
    JSparseB[13] = 8849557522123.8945*dwdx26;
    JSparseB[14] = -8849557522123.8945*dwdx33;
    JSparseB[15] = -8849557522123.8945*dwdx39;
    JSparseB[16] = 11152416356877.322*dwdx1;
    JSparseB[17] = 11152416356877.322*dwdx8 + 3717472118959.1074*dwdx9;
    JSparseB[18] = 3717472118959.1074*dwdx11;
    JSparseB[19] = 11152416356877.322*dwdx12;
    JSparseB[20] = 26548672566371.684*dwdx4;
    JSparseB[21] = -8849557522123.8945*dwdx9;
    JSparseB[22] = 26548672566371.684*dwdx10 - 8849557522123.8945*dwdx11;
    JSparseB[23] = 26548672566371.684*dwdx14;
    JSparseB[24] = 3717472118959.1074*dwdx1;
    JSparseB[25] = 3717472118959.1074*dwdx8;
    JSparseB[26] = 3717472118959.1074*dwdx12 + 3717472118959.1074*dwdx13;
    JSparseB[27] = 3717472118959.1074*dwdx15;
    JSparseB[28] = 8849557522123.8945*dwdx4;
    JSparseB[29] = 8849557522123.8945*dwdx10;
    JSparseB[30] = -8849557522123.8945*dwdx13;
    JSparseB[31] = 8849557522123.8945*dwdx14 - 8849557522123.8945*dwdx15;
    JSparseB[32] = 3717472118959.1074*dwdx0;
    JSparseB[33] = 3717472118959.1074*dwdx16 - 3717472118959.1074*dwdx17 + 3717472118959.1074*dwdx18;
    JSparseB[34] = -3717472118959.1074*dwdx21;
    JSparseB[35] = 3717472118959.1074*dwdx22;
    JSparseB[36] = 3717472118959.1074*dwdx29;
    JSparseB[37] = -8849557522123.8945*dwdx5;
    JSparseB[38] = 8849557522123.8945*dwdx17;
    JSparseB[39] = -8849557522123.8945*dwdx19 + 8849557522123.8945*dwdx20 + 8849557522123.8945*dwdx21;
    JSparseB[40] = -8849557522123.8945*dwdx26;
    JSparseB[41] = 8849557522123.8945*dwdx32;
    JSparseB[42] = -3717472118959.1074*dwdx0;
    JSparseB[43] = -3717472118959.1074*dwdx16;
    JSparseB[44] = -3717472118959.1074*dwdx22 - 3717472118959.1074*dwdx23 - 3717472118959.1074*dwdx24;
    JSparseB[45] = -3717472118959.1074*dwdx27;
    JSparseB[46] = -3717472118959.1074*dwdx35;
    JSparseB[47] = 8849557522123.8945*dwdx5;
    JSparseB[48] = 8849557522123.8945*dwdx19;
    JSparseB[49] = 8849557522123.8945*dwdx24;
    JSparseB[50] = -8849557522123.8945*dwdx25 + 8849557522123.8945*dwdx26 + 8849557522123.8945*dwdx27;
    JSparseB[51] = -8849557522123.8945*dwdx37;
    JSparseB[52] = 3717472118959.1074*dwdx3;
    JSparseB[53] = -3717472118959.1074*dwdx18;
    JSparseB[54] = 3717472118959.1074*dwdx28 - 3717472118959.1074*dwdx29 + 3717472118959.1074*dwdx30;
    JSparseB[55] = 3717472118959.1074*dwdx31;
    JSparseB[56] = 3717472118959.1074*dwdx36;
    JSparseB[57] = 8849557522123.8945*dwdx7;
    JSparseB[58] = -8849557522123.8945*dwdx20;
    JSparseB[59] = -8849557522123.8945*dwdx28;
    JSparseB[60] = -8849557522123.8945*dwdx31 - 8849557522123.8945*dwdx32 + 8849557522123.8945*dwdx33;
    JSparseB[61] = 8849557522123.8945*dwdx39;
    JSparseB[62] = -3717472118959.1074*dwdx3;
    JSparseB[63] = 3717472118959.1074*dwdx23;
    JSparseB[64] = -3717472118959.1074*dwdx30;
    JSparseB[65] = 3717472118959.1074*dwdx34 + 3717472118959.1074*dwdx35 - 3717472118959.1074*dwdx36;
    JSparseB[66] = 3717472118959.1074*dwdx38;
    JSparseB[67] = -8849557522123.8945*dwdx7;
    JSparseB[68] = 8849557522123.8945*dwdx25;
    JSparseB[69] = -8849557522123.8945*dwdx33;
    JSparseB[70] = -8849557522123.8945*dwdx34;
    JSparseB[71] = 8849557522123.8945*dwdx37 - 8849557522123.8945*dwdx38 - 8849557522123.8945*dwdx39;
}