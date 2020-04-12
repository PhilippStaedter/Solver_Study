#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_leloup1_Fig2A(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -1.0*dwdx0 + 1.0*dwdx1 - 1.0*dwdx2;
    JSparse[1] = 1.0*dwdx2;
    JSparse[2] = -1.0*dwdx1;
    JSparse[3] = -1.0*dwdx1;
    JSparse[4] = -1.0*dwdx5;
    JSparse[5] = 1.0*dwdx5 - 1.0*dwdx6;
    JSparse[6] = 1.0*dwdx3;
    JSparse[7] = 1.0*dwdx4;
    JSparse[8] = -1.0*dwdx7;
    JSparse[9] = 1.0*dwdx8;
    JSparse[10] = -1.0*dwdx9;
    JSparse[11] = 1.0*dwdx10;
    JSparse[12] = -1.0*dwdx11 - 1.0*dwdx12;
    JSparse[13] = 1.0*dwdx12;
    JSparse[14] = 1.0*dwdx14;
    JSparse[15] = -1.0*dwdx13 - 1.0*dwdx14 - 1.0*dwdx15;
    JSparse[16] = 1.0*dwdx15;
    JSparse[17] = 1.0*dwdx18;
    JSparse[18] = 1.0*dwdx17;
    JSparse[19] = -1.0*dwdx16 - 1.0*dwdx17 - 1.0*dwdx18;
    JSparse[20] = -1.0*dwdx18;
    JSparse[21] = -1.0*dwdx19 - 1.0*dwdx20;
    JSparse[22] = 1.0*dwdx20;
    JSparse[23] = 1.0*dwdx22;
    JSparse[24] = -1.0*dwdx21 - 1.0*dwdx22 - 1.0*dwdx23;
    JSparse[25] = 1.0*dwdx23;
    JSparse[26] = 1.0*dwdx24;
    JSparse[27] = -1.0*dwdx24;
    JSparse[28] = 1.0*dwdx26;
    JSparse[29] = -1.0*dwdx24 - 1.0*dwdx25 - 1.0*dwdx26;
}