#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Levchenko2000a(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = MAPKP;
    y[1] = MEKP;
    y[2] = RAFK;
    y[3] = RAFP;
    y[4] = K_1_0;
    y[5] = K_1_1;
    y[6] = K_1_2;
    y[7] = K_2_0;
    y[8] = K_2_1;
    y[9] = K_2_2;
    y[10] = K_3_0;
    y[11] = K_3_1;
    y[12] = K_K_1_0_2_2;
    y[13] = K_K_1_1_2_2;
    y[14] = K_K_2_0_3_1;
    y[15] = K_K_2_1_3_1;
    y[16] = K_MAPKP_1_1;
    y[17] = K_MAPKP_1_2;
    y[18] = K_MEKP_2_1;
    y[19] = K_MEKP_2_2;
    y[20] = K_RAFK_3_0;
    y[21] = K_RAFP_3_1;
    y[22] = S_m1_m1_m1;
    y[23] = S_m1_m1_0;
    y[24] = S_m1_m1_1;
    y[25] = S_m1_0_m1;
    y[26] = S_m1_0_0;
    y[27] = S_m1_0_1;
    y[28] = S_m1_1_m1;
    y[29] = S_m1_1_0;
    y[30] = S_m1_1_1;
    y[31] = S_m1_2_m1;
    y[32] = S_m1_2_0;
    y[33] = S_m1_2_1;
    y[34] = S_0_m1_m1;
    y[35] = S_0_m1_0;
    y[36] = S_0_m1_1;
    y[37] = S_0_0_m1;
    y[38] = S_0_0_0;
    y[39] = S_0_0_1;
    y[40] = S_0_1_m1;
    y[41] = S_0_1_0;
    y[42] = S_0_1_1;
    y[43] = S_0_2_m1;
    y[44] = S_0_2_0;
    y[45] = S_0_2_1;
    y[46] = S_1_m1_m1;
    y[47] = S_1_m1_0;
    y[48] = S_1_m1_1;
    y[49] = S_1_0_m1;
    y[50] = S_1_0_0;
    y[51] = S_1_0_1;
    y[52] = S_1_1_m1;
    y[53] = S_1_1_0;
    y[54] = S_1_1_1;
    y[55] = S_1_2_m1;
    y[56] = S_1_2_0;
    y[57] = S_1_2_1;
    y[58] = S_2_m1_m1;
    y[59] = S_2_m1_0;
    y[60] = S_2_m1_1;
    y[61] = S_2_0_m1;
    y[62] = S_2_0_0;
    y[63] = S_2_0_1;
    y[64] = S_2_1_m1;
    y[65] = S_2_1_0;
    y[66] = S_2_1_1;
    y[67] = S_2_2_m1;
    y[68] = S_2_2_0;
    y[69] = S_2_2_1;
    y[70] = S_RAFK_m1_m1_0;
    y[71] = S_RAFK_m1_0_0;
    y[72] = S_RAFK_m1_1_0;
    y[73] = S_RAFK_m1_2_0;
    y[74] = S_RAFK_0_m1_0;
    y[75] = S_RAFK_0_0_0;
    y[76] = S_RAFK_0_1_0;
    y[77] = S_RAFK_0_2_0;
    y[78] = S_RAFK_1_m1_0;
    y[79] = S_RAFK_1_0_0;
    y[80] = S_RAFK_1_1_0;
    y[81] = S_RAFK_1_2_0;
    y[82] = S_RAFK_2_m1_0;
    y[83] = S_RAFK_2_0_0;
    y[84] = S_RAFK_2_1_0;
    y[85] = S_RAFK_2_2_0;
}