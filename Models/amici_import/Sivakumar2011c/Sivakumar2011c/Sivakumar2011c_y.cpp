#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Sivakumar2011c(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = s5;
    y[1] = s16;
    y[2] = s1;
    y[3] = s27;
    y[4] = s28;
    y[5] = s30;
    y[6] = s31;
    y[7] = s32;
    y[8] = s33;
    y[9] = s37;
    y[10] = s46;
    y[11] = s75;
    y[12] = s101;
    y[13] = s102;
    y[14] = s107;
    y[15] = s121;
    y[16] = s155;
    y[17] = s164;
    y[18] = s171;
    y[19] = s172;
    y[20] = s173;
    y[21] = s170;
    y[22] = s195;
    y[23] = s174;
    y[24] = s239;
    y[25] = s36;
    y[26] = s123;
    y[27] = s129;
    y[28] = s159;
    y[29] = s232;
    y[30] = s176;
    y[31] = s179;
    y[32] = s183;
    y[33] = s188;
    y[34] = s245;
    y[35] = s252;
    y[36] = s268;
    y[37] = s260;
    y[38] = s270;
    y[39] = s275;
    y[40] = s278;
    y[41] = s286;
    y[42] = s288;
    y[43] = s292;
    y[44] = s61;
    y[45] = s259;
    y[46] = s266;
    y[47] = s267;
    y[48] = s304;
    y[49] = s305;
}