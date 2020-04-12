#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_kolodkin5(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = 1.0*dwdx0;
    JSparse[1] = 1.0*dwdx1;
    JSparse[2] = -1.0*dwdx1;
    JSparse[3] = 1.0*dwdx2 - 1.0*dwdx3;
    JSparse[4] = -1.0*dwdx3;
    JSparse[5] = 1.0*dwdx3;
    JSparse[6] = -1.0*dwdx4;
    JSparse[7] = -1.0*dwdx4 + 3.4444444444444402*dwdx5;
    JSparse[8] = 1.0*dwdx4;
    JSparse[9] = -1.0*dwdx5;
    JSparse[10] = 1.0*dwdx6 - 1.0*dwdx7;
    JSparse[11] = 3.4444444444444402*dwdx7;
    JSparse[12] = -1.0*dwdx6;
    JSparse[13] = -1.0*dwdx9;
    JSparse[14] = -1.0*dwdx9;
    JSparse[15] = -1.0*dwdx10;
    JSparse[16] = 3.4444444444444402*dwdx10 - 1.0*dwdx8 + 1.0*dwdx9;
    JSparse[17] = -1.0*dwdx8;
    JSparse[18] = 1.0*dwdx8;
    JSparse[19] = 3.4444444444444402*dwdx12;
    JSparse[20] = 1.0*dwdx11;
    JSparse[21] = -1.0*dwdx11 - 1.0*dwdx12;
    JSparse[22] = -1.0*dwdx13;
    JSparse[23] = -1.0*dwdx13;
    JSparse[24] = 1.0*dwdx13;
    JSparse[25] = -1.0*dwdx14;
    JSparse[26] = -1.0*dwdx14;
    JSparse[27] = 1.0*dwdx14;
}