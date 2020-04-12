#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_bucher1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 70.422535211267601*dwdx0;
    JSparse[1] = -70.422535211267601*dwdx0;
    JSparse[2] = 70.422535211267601*dwdx1;
    JSparse[3] = -70.422535211267601*dwdx1 - 70.422535211267601*dwdx2 - 70.422535211267601*dwdx3 - 70.422535211267601*dwdx6 - 70.422535211267601*dwdx7;
    JSparse[4] = 0.5*dwdx6;
    JSparse[5] = 70.422535211267601*dwdx2;
    JSparse[6] = 70.422535211267601*dwdx3;
    JSparse[7] = -70.422535211267601*dwdx4 - 70.422535211267601*dwdx5 + 70.422535211267601*dwdx7;
    JSparse[8] = 70.422535211267601*dwdx4;
    JSparse[9] = 70.422535211267601*dwdx5;
    JSparse[10] = 70.422535211267601*dwdx8;
    JSparse[11] = -0.5*dwdx8 - 0.5*dwdx9;
    JSparse[12] = 0.5*dwdx9;
    JSparse[13] = 70.422535211267601*dwdx10;
    JSparse[14] = -70.422535211267601*dwdx10;
    JSparse[15] = 70.422535211267601*dwdx11;
    JSparse[16] = -70.422535211267601*dwdx11 - 70.422535211267601*dwdx12 - 70.422535211267601*dwdx13;
    JSparse[17] = 0.5*dwdx13;
    JSparse[18] = 70.422535211267601*dwdx12;
    JSparse[19] = 70.422535211267601*dwdx14;
    JSparse[20] = -0.5*dwdx14 - 0.5*dwdx15;
    JSparse[21] = 0.5*dwdx15;
    JSparse[22] = 70.422535211267601*dwdx16;
    JSparse[23] = -70.422535211267601*dwdx16;
    JSparse[24] = 70.422535211267601*dwdx17;
    JSparse[25] = -70.422535211267601*dwdx17 - 70.422535211267601*dwdx18 - 70.422535211267601*dwdx19;
    JSparse[26] = 0.5*dwdx19;
    JSparse[27] = 70.422535211267601*dwdx18;
    JSparse[28] = 70.422535211267601*dwdx20;
    JSparse[29] = -0.5*dwdx20 - 0.5*dwdx21;
    JSparse[30] = 0.5*dwdx21;
    JSparse[31] = 70.422535211267601*dwdx22;
    JSparse[32] = -70.422535211267601*dwdx22;
    JSparse[33] = -70.422535211267601*dwdx24 - 70.422535211267601*dwdx25 + 70.422535211267601*dwdx29;
    JSparse[34] = 70.422535211267601*dwdx24;
    JSparse[35] = 70.422535211267601*dwdx25;
    JSparse[36] = 70.422535211267601*dwdx23;
    JSparse[37] = -70.422535211267601*dwdx23 - 70.422535211267601*dwdx26 - 70.422535211267601*dwdx27 - 70.422535211267601*dwdx28 - 70.422535211267601*dwdx29;
    JSparse[38] = 0.5*dwdx28;
    JSparse[39] = 70.422535211267601*dwdx26;
    JSparse[40] = 70.422535211267601*dwdx27;
    JSparse[41] = 70.422535211267601*dwdx30;
    JSparse[42] = -0.5*dwdx30;
    JSparse[43] = 70.422535211267601*dwdx31;
    JSparse[44] = -70.422535211267601*dwdx31;
    JSparse[45] = 70.422535211267601*dwdx32;
    JSparse[46] = -70.422535211267601*dwdx32 - 70.422535211267601*dwdx33;
    JSparse[47] = 0.5*dwdx33;
    JSparse[48] = 70.422535211267601*dwdx34;
    JSparse[49] = -0.5*dwdx34;
    JSparse[50] = 70.422535211267601*dwdx35;
    JSparse[51] = -70.422535211267601*dwdx35;
    JSparse[52] = 70.422535211267601*dwdx36;
    JSparse[53] = -70.422535211267601*dwdx36 - 70.422535211267601*dwdx37;
    JSparse[54] = 0.5*dwdx37;
    JSparse[55] = 70.422535211267601*dwdx38;
    JSparse[56] = -0.5*dwdx38;
}