#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_model2_mcauley1(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = 1.0;
    dxdotdw[1] = 1.0;
    dxdotdw[2] = -1.0;
    dxdotdw[3] = 1.0;
    dxdotdw[4] = 1.0;
    dxdotdw[5] = -1.0;
    dxdotdw[6] = -1.0;
    dxdotdw[7] = 1.0;
    dxdotdw[8] = 1.0;
    dxdotdw[9] = 1.0;
    dxdotdw[10] = -1.0;
    dxdotdw[11] = 1.0;
    dxdotdw[12] = -1.0;
    dxdotdw[13] = 1.0;
    dxdotdw[14] = -1.0;
    dxdotdw[15] = 1.0;
    dxdotdw[16] = -1.0;
    dxdotdw[17] = 1.0;
    dxdotdw[18] = 1.0;
    dxdotdw[19] = -1.0;
    dxdotdw[20] = 1.0;
    dxdotdw[21] = -1.0;
    dxdotdw[22] = 1.0;
    dxdotdw[23] = -1.0;
    dxdotdw[24] = 1.0;
    dxdotdw[25] = -1.0;
    dxdotdw[26] = 1.0;
    dxdotdw[27] = 1.0;
    dxdotdw[28] = -1.0;
    dxdotdw[29] = 1.0;
    dxdotdw[30] = -1.0;
    dxdotdw[31] = 1.0;
    dxdotdw[32] = -1.0;
    dxdotdw[33] = 1.0;
    dxdotdw[34] = -1.0;
    dxdotdw[35] = 1.0;
    dxdotdw[36] = 1.0;
    dxdotdw[37] = -1.0;
    dxdotdw[38] = -1.0;
    dxdotdw[39] = 1.0;
    dxdotdw[40] = -1.0;
    dxdotdw[41] = 1.0;
    dxdotdw[42] = -1.0;
    dxdotdw[43] = -1.0;
    dxdotdw[44] = 1.0;
    dxdotdw[45] = 1.0;
    dxdotdw[46] = 1.0;
    dxdotdw[47] = -1.0;
    dxdotdw[48] = 1.0;
    dxdotdw[49] = -1.0;
    dxdotdw[50] = -1.0;
    dxdotdw[51] = 1.0;
    dxdotdw[52] = 1.0;
    dxdotdw[53] = -1.0;
    dxdotdw[54] = -1.0;
    dxdotdw[55] = 1.0;
    dxdotdw[56] = 1.0;
    dxdotdw[57] = -1.0;
    dxdotdw[58] = -1.0;
    dxdotdw[59] = 1.0;
    dxdotdw[60] = -1.0;
    dxdotdw[61] = 1.0;
    dxdotdw[62] = 1.0;
}