#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dxdotdw_Leber2015(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdotdw[0] = -14.285714285714285;
    dxdotdw[1] = -1.0;
    dxdotdw[2] = -14.285714285714285;
    dxdotdw[3] = -14.285714285714285;
    dxdotdw[4] = -1.0;
    dxdotdw[5] = -0.25;
    dxdotdw[6] = 0.25;
    dxdotdw[7] = -14.285714285714285;
    dxdotdw[8] = 1.0;
    dxdotdw[9] = -1.0;
    dxdotdw[10] = 14.285714285714285;
    dxdotdw[11] = -1.0;
    dxdotdw[12] = 1.0;
    dxdotdw[13] = 1.0;
    dxdotdw[14] = 14.285714285714285;
    dxdotdw[15] = -1.0;
    dxdotdw[16] = 14.285714285714285;
    dxdotdw[17] = -1.0;
    dxdotdw[18] = -14.285714285714285;
    dxdotdw[19] = 14.285714285714285;
    dxdotdw[20] = 14.285714285714285;
    dxdotdw[21] = -1.0;
    dxdotdw[22] = -0.25;
    dxdotdw[23] = 0.25;
    dxdotdw[24] = 0.25;
    dxdotdw[25] = -0.25;
    dxdotdw[26] = 0.25;
    dxdotdw[27] = -0.25;
    dxdotdw[28] = -1.0;
    dxdotdw[29] = 1.0;
    dxdotdw[30] = 0.25;
    dxdotdw[31] = -0.25;
    dxdotdw[32] = -1.0;
    dxdotdw[33] = 1.0;
    dxdotdw[34] = -1.0;
    dxdotdw[35] = 1.0;
    dxdotdw[36] = -1.0;
    dxdotdw[37] = -1.0;
    dxdotdw[38] = 1.0;
    dxdotdw[39] = -1.0;
    dxdotdw[40] = 1.0;
    dxdotdw[41] = -1.0;
    dxdotdw[42] = 1.0;
    dxdotdw[43] = -1.0;
    dxdotdw[44] = -1.0;
    dxdotdw[45] = 0.25;
    dxdotdw[46] = -0.25;
}