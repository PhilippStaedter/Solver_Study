#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparse_model0_zi1(realtype *JSparse, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JSparse[0] = -952.38095238095241*dwdx0 - 952.38095238095241*dwdx1;
    JSparse[1] = 952.38095238095241*dwdx0;
    JSparse[2] = 952.38095238095241*dwdx0;
    JSparse[3] = 1.0*dwdx0;
    JSparse[4] = -952.38095238095241*dwdx2 - 952.38095238095241*dwdx3;
    JSparse[5] = -952.38095238095241*dwdx4;
    JSparse[6] = -952.38095238095241*dwdx4;
    JSparse[7] = 952.38095238095241*dwdx4;
    JSparse[8] = 952.38095238095241*dwdx2;
    JSparse[9] = 952.38095238095241*dwdx2;
    JSparse[10] = 1.0*dwdx2;
    JSparse[11] = 952.38095238095241*dwdx5;
    JSparse[12] = 952.38095238095241*dwdx6;
    JSparse[13] = -952.38095238095241*dwdx5 - 952.38095238095241*dwdx6;
    JSparse[14] = -952.38095238095241*dwdx7 - 952.38095238095241*dwdx8;
    JSparse[15] = 2857.1428571428573*dwdx7;
    JSparse[16] = -952.38095238095241*dwdx8;
    JSparse[17] = 952.38095238095241*dwdx8;
    JSparse[18] = 952.38095238095241*dwdx9;
    JSparse[19] = -2857.1428571428573*dwdx9;
    JSparse[20] = -952.38095238095241*dwdx10;
    JSparse[21] = -952.38095238095241*dwdx10 - 952.38095238095241*dwdx11;
    JSparse[22] = 2857.1428571428573*dwdx11;
    JSparse[23] = 952.38095238095241*dwdx10;
    JSparse[24] = 952.38095238095241*dwdx12;
    JSparse[25] = -2857.1428571428573*dwdx12;
    JSparse[26] = -952.38095238095241*dwdx13;
    JSparse[27] = 2857.1428571428573*dwdx13;
    JSparse[28] = -952.38095238095241*dwdx15;
    JSparse[29] = 2857.1428571428573*dwdx14;
    JSparse[30] = 2857.1428571428573*dwdx14;
    JSparse[31] = -2857.1428571428573*dwdx14;
    JSparse[32] = -952.38095238095241*dwdx16;
    JSparse[33] = 952.38095238095241*dwdx16;
    JSparse[34] = -952.38095238095241*dwdx17 - 952.38095238095241*dwdx18;
    JSparse[35] = 952.38095238095241*dwdx18;
    JSparse[36] = 952.38095238095241*dwdx19;
    JSparse[37] = 952.38095238095241*dwdx20;
    JSparse[38] = 952.38095238095241*dwdx21;
    JSparse[39] = -952.38095238095241*dwdx19 - 952.38095238095241*dwdx20 - 952.38095238095241*dwdx21;
    JSparse[40] = -952.38095238095241*dwdx19;
    JSparse[41] = -1.0*dwdx19;
    JSparse[42] = -952.38095238095241*dwdx22;
    JSparse[43] = 952.38095238095241*dwdx22;
    JSparse[44] = -952.38095238095241*dwdx23 - 952.38095238095241*dwdx24;
    JSparse[45] = 952.38095238095241*dwdx23;
    JSparse[46] = 1.0*dwdx23;
    JSparse[47] = 952.38095238095241*dwdx27;
    JSparse[48] = -952.38095238095241*dwdx27;
    JSparse[49] = 952.38095238095241*dwdx25;
    JSparse[50] = 952.38095238095241*dwdx26;
    JSparse[51] = -952.38095238095241*dwdx25 - 952.38095238095241*dwdx26 - 952.38095238095241*dwdx27;
    JSparse[52] = -1.0*dwdx27;
    JSparse[53] = 952.38095238095241*dwdx28;
    JSparse[54] = -952.38095238095241*dwdx28;
    JSparse[55] = -952.38095238095241*dwdx28;
    JSparse[56] = -1.0*dwdx28;
}