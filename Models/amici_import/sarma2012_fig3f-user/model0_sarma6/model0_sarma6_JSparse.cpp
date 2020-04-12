#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_sarma6(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0;
    JSparse[1] = 1.0*dwdx0;
    JSparse[2] = 1.0*dwdx3;
    JSparse[3] = -1.0*dwdx3;
    JSparse[4] = -1.0*dwdx1;
    JSparse[5] = 1.0*dwdx1 - 1.0*dwdx2;
    JSparse[6] = 1.0*dwdx2;
    JSparse[7] = 1.0*dwdx5;
    JSparse[8] = 1.0*dwdx4 - 1.0*dwdx5;
    JSparse[9] = -1.0*dwdx4;
    JSparse[10] = -1.0*dwdx6;
    JSparse[11] = 1.0*dwdx6 - 1.0*dwdx7;
    JSparse[12] = 1.0*dwdx7;
    JSparse[13] = 1.0*dwdx11 - 1.0*dwdx8;
    JSparse[14] = 1.0*dwdx10 - 1.0*dwdx11 + 1.0*dwdx8 - 1.0*dwdx9;
    JSparse[15] = -1.0*dwdx10 + 1.0*dwdx9;
    JSparse[16] = 1.0*dwdx15;
    JSparse[17] = 1.0*dwdx14 - 1.0*dwdx15;
    JSparse[18] = -1.0*dwdx14;
    JSparse[19] = -1.0*dwdx12;
    JSparse[20] = 1.0*dwdx12 - 1.0*dwdx13;
    JSparse[21] = 1.0*dwdx13;
    JSparse[22] = -1.0*dwdx16;
    JSparse[23] = 1.0*dwdx16 - 1.0*dwdx17;
    JSparse[24] = 1.0*dwdx17;
    JSparse[25] = -1.0*dwdx18 + 1.0*dwdx21;
    JSparse[26] = 1.0*dwdx18 - 1.0*dwdx19 + 1.0*dwdx20 - 1.0*dwdx21;
    JSparse[27] = 1.0*dwdx19 - 1.0*dwdx20;
    JSparse[28] = -1.0*dwdx22;
    JSparse[29] = 1.0*dwdx22;
    JSparse[30] = -1.0*dwdx23;
    JSparse[31] = 1.0*dwdx23 - 1.0*dwdx24;
    JSparse[32] = 1.0*dwdx24;
    JSparse[33] = 1.0*dwdx26;
    JSparse[34] = 1.0*dwdx25 - 1.0*dwdx26;
    JSparse[35] = -1.0*dwdx25;
    JSparse[36] = 1.0*dwdx27;
    JSparse[37] = -1.0*dwdx27;
    JSparse[38] = 1.0*dwdx29;
    JSparse[39] = 1.0*dwdx28 - 1.0*dwdx29;
    JSparse[40] = -1.0*dwdx28;
}