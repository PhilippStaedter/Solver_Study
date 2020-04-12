#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"
#include "dwdx.h"

void JSparseB_model0_zi1(realtype *JSparseB, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *xB, const realtype *w, const realtype *dwdx){
    JSparseB[0] = 952.38095238095241*dwdx0 + 952.38095238095241*dwdx1;
    JSparseB[1] = -952.38095238095241*dwdx5;
    JSparseB[2] = 952.38095238095241*dwdx15;
    JSparseB[3] = 952.38095238095241*dwdx2 + 952.38095238095241*dwdx3;
    JSparseB[4] = -952.38095238095241*dwdx6;
    JSparseB[5] = 952.38095238095241*dwdx5 + 952.38095238095241*dwdx6;
    JSparseB[6] = -952.38095238095241*dwdx19;
    JSparseB[7] = -952.38095238095241*dwdx27;
    JSparseB[8] = -952.38095238095241*dwdx28;
    JSparseB[9] = 952.38095238095241*dwdx4;
    JSparseB[10] = 952.38095238095241*dwdx7 + 952.38095238095241*dwdx8;
    JSparseB[11] = -952.38095238095241*dwdx9;
    JSparseB[12] = 952.38095238095241*dwdx10;
    JSparseB[13] = -2857.1428571428573*dwdx7;
    JSparseB[14] = 2857.1428571428573*dwdx9;
    JSparseB[15] = -2857.1428571428573*dwdx14;
    JSparseB[16] = 952.38095238095241*dwdx4;
    JSparseB[17] = 952.38095238095241*dwdx8;
    JSparseB[18] = 952.38095238095241*dwdx10 + 952.38095238095241*dwdx11;
    JSparseB[19] = -952.38095238095241*dwdx12;
    JSparseB[20] = -2857.1428571428573*dwdx11;
    JSparseB[21] = 2857.1428571428573*dwdx12;
    JSparseB[22] = -2857.1428571428573*dwdx14;
    JSparseB[23] = -952.38095238095241*dwdx4;
    JSparseB[24] = -952.38095238095241*dwdx8;
    JSparseB[25] = -952.38095238095241*dwdx10;
    JSparseB[26] = 952.38095238095241*dwdx13;
    JSparseB[27] = -2857.1428571428573*dwdx13;
    JSparseB[28] = 2857.1428571428573*dwdx14;
    JSparseB[29] = 952.38095238095241*dwdx16;
    JSparseB[30] = -952.38095238095241*dwdx20;
    JSparseB[31] = 952.38095238095241*dwdx17 + 952.38095238095241*dwdx18;
    JSparseB[32] = -952.38095238095241*dwdx21;
    JSparseB[33] = -952.38095238095241*dwdx0;
    JSparseB[34] = -952.38095238095241*dwdx2;
    JSparseB[35] = -952.38095238095241*dwdx16;
    JSparseB[36] = -952.38095238095241*dwdx18;
    JSparseB[37] = 952.38095238095241*dwdx19 + 952.38095238095241*dwdx20 + 952.38095238095241*dwdx21;
    JSparseB[38] = 952.38095238095241*dwdx27;
    JSparseB[39] = 952.38095238095241*dwdx28;
    JSparseB[40] = 952.38095238095241*dwdx22;
    JSparseB[41] = -952.38095238095241*dwdx25;
    JSparseB[42] = 952.38095238095241*dwdx23 + 952.38095238095241*dwdx24;
    JSparseB[43] = -952.38095238095241*dwdx26;
    JSparseB[44] = -952.38095238095241*dwdx0;
    JSparseB[45] = -952.38095238095241*dwdx2;
    JSparseB[46] = 952.38095238095241*dwdx19;
    JSparseB[47] = -952.38095238095241*dwdx22;
    JSparseB[48] = -952.38095238095241*dwdx23;
    JSparseB[49] = 952.38095238095241*dwdx25 + 952.38095238095241*dwdx26 + 952.38095238095241*dwdx27;
    JSparseB[50] = 952.38095238095241*dwdx28;
    JSparseB[51] = -1.0*dwdx0;
    JSparseB[52] = -1.0*dwdx2;
    JSparseB[53] = 1.0*dwdx19;
    JSparseB[54] = -1.0*dwdx23;
    JSparseB[55] = 1.0*dwdx27;
    JSparseB[56] = 1.0*dwdx28;
}