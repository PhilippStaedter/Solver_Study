#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_bucher1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = -70.422535211267601*dwdx0;
    JSparseB[1] = -70.422535211267601*dwdx1;
    JSparseB[2] = 70.422535211267601*dwdx0;
    JSparseB[3] = 70.422535211267601*dwdx1 + 70.422535211267601*dwdx2 + 70.422535211267601*dwdx3 + 70.422535211267601*dwdx6 + 70.422535211267601*dwdx7;
    JSparseB[4] = -70.422535211267601*dwdx8;
    JSparseB[5] = 70.422535211267601*dwdx24 + 70.422535211267601*dwdx25 - 70.422535211267601*dwdx29;
    JSparseB[6] = -0.5*dwdx6;
    JSparseB[7] = 0.5*dwdx8 + 0.5*dwdx9;
    JSparseB[8] = -70.422535211267601*dwdx10;
    JSparseB[9] = -70.422535211267601*dwdx11;
    JSparseB[10] = -70.422535211267601*dwdx2;
    JSparseB[11] = 70.422535211267601*dwdx10;
    JSparseB[12] = 70.422535211267601*dwdx11 + 70.422535211267601*dwdx12 + 70.422535211267601*dwdx13;
    JSparseB[13] = -70.422535211267601*dwdx14;
    JSparseB[14] = -70.422535211267601*dwdx24;
    JSparseB[15] = -0.5*dwdx13;
    JSparseB[16] = 0.5*dwdx14 + 0.5*dwdx15;
    JSparseB[17] = -70.422535211267601*dwdx16;
    JSparseB[18] = -70.422535211267601*dwdx17;
    JSparseB[19] = -70.422535211267601*dwdx3;
    JSparseB[20] = 70.422535211267601*dwdx16;
    JSparseB[21] = 70.422535211267601*dwdx17 + 70.422535211267601*dwdx18 + 70.422535211267601*dwdx19;
    JSparseB[22] = -70.422535211267601*dwdx20;
    JSparseB[23] = -70.422535211267601*dwdx25;
    JSparseB[24] = -0.5*dwdx19;
    JSparseB[25] = 0.5*dwdx20 + 0.5*dwdx21;
    JSparseB[26] = -70.422535211267601*dwdx22;
    JSparseB[27] = -70.422535211267601*dwdx23;
    JSparseB[28] = 70.422535211267601*dwdx4 + 70.422535211267601*dwdx5 - 70.422535211267601*dwdx7;
    JSparseB[29] = 70.422535211267601*dwdx22;
    JSparseB[30] = 70.422535211267601*dwdx23 + 70.422535211267601*dwdx26 + 70.422535211267601*dwdx27 + 70.422535211267601*dwdx28 + 70.422535211267601*dwdx29;
    JSparseB[31] = -70.422535211267601*dwdx30;
    JSparseB[32] = -0.5*dwdx9;
    JSparseB[33] = -0.5*dwdx28;
    JSparseB[34] = 0.5*dwdx30;
    JSparseB[35] = -70.422535211267601*dwdx31;
    JSparseB[36] = -70.422535211267601*dwdx32;
    JSparseB[37] = -70.422535211267601*dwdx4;
    JSparseB[38] = -70.422535211267601*dwdx12;
    JSparseB[39] = -70.422535211267601*dwdx26;
    JSparseB[40] = 70.422535211267601*dwdx31;
    JSparseB[41] = 70.422535211267601*dwdx32 + 70.422535211267601*dwdx33;
    JSparseB[42] = -70.422535211267601*dwdx34;
    JSparseB[43] = -0.5*dwdx15;
    JSparseB[44] = -0.5*dwdx33;
    JSparseB[45] = 0.5*dwdx34;
    JSparseB[46] = -70.422535211267601*dwdx35;
    JSparseB[47] = -70.422535211267601*dwdx36;
    JSparseB[48] = -70.422535211267601*dwdx5;
    JSparseB[49] = -70.422535211267601*dwdx18;
    JSparseB[50] = -70.422535211267601*dwdx27;
    JSparseB[51] = 70.422535211267601*dwdx35;
    JSparseB[52] = 70.422535211267601*dwdx36 + 70.422535211267601*dwdx37;
    JSparseB[53] = -70.422535211267601*dwdx38;
    JSparseB[54] = -0.5*dwdx21;
    JSparseB[55] = -0.5*dwdx37;
    JSparseB[56] = 0.5*dwdx38;
}