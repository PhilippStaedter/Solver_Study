#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_bachmann(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -2.5;
    dxdotdw[1] = 2.5;
    dxdotdw[2] = -2.5;
    dxdotdw[3] = -2.5;
    dxdotdw[4] = 2.5;
    dxdotdw[5] = 2.5;
    dxdotdw[6] = -2.5;
    dxdotdw[7] = -2.5;
    dxdotdw[8] = 2.5;
    dxdotdw[9] = -2.5;
    dxdotdw[10] = 2.5;
    dxdotdw[11] = 3.6363636363636362;
    dxdotdw[12] = -2.5;
    dxdotdw[13] = 2.5;
    dxdotdw[14] = -3.6363636363636362;
    dxdotdw[15] = 3.6363636363636362;
    dxdotdw[16] = -3.6363636363636362;
    dxdotdw[17] = 3.6363636363636362;
    dxdotdw[18] = -3.6363636363636362;
    dxdotdw[19] = 3.6363636363636362;
    dxdotdw[20] = 2.5;
    dxdotdw[21] = -2.5;
    dxdotdw[22] = -3.6363636363636362;
    dxdotdw[23] = 3.6363636363636362;
    dxdotdw[24] = -3.6363636363636362;
    dxdotdw[25] = 3.6363636363636362;
    dxdotdw[26] = 2.5;
    dxdotdw[27] = -3.6363636363636362;
    dxdotdw[28] = -2.5;
    dxdotdw[29] = 2.5;
    dxdotdw[30] = -2.5;
    dxdotdw[31] = 2.5;
    dxdotdw[32] = 3.6363636363636362;
    dxdotdw[33] = -3.6363636363636362;
    dxdotdw[34] = 3.6363636363636362;
    dxdotdw[35] = -3.6363636363636362;
    dxdotdw[36] = 3.6363636363636362;
    dxdotdw[37] = -2.5;
    dxdotdw[38] = 2.5;
    dxdotdw[39] = -3.6363636363636362;
    dxdotdw[40] = 3.6363636363636362;
    dxdotdw[41] = -3.6363636363636362;
    dxdotdw[42] = 3.6363636363636362;
    dxdotdw[43] = 2.5;
    dxdotdw[44] = -3.6363636363636362;
    dxdotdw[45] = -2.5;
    dxdotdw[46] = 2.5;
    dxdotdw[47] = -2.5;
    dxdotdw[48] = 2.5;
    dxdotdw[49] = -2.5;
    dxdotdw[50] = 2.5;
    dxdotdw[51] = 2.5;
    dxdotdw[52] = -2.5;
    dxdotdw[53] = 2.5;
    dxdotdw[54] = -2.5;
    dxdotdw[55] = 2.5;
    dxdotdw[56] = -2.5;
    dxdotdw[57] = 2.5;
    dxdotdw[58] = -2.5;
    dxdotdw[59] = 2.5;
    dxdotdw[60] = -2.5;
}