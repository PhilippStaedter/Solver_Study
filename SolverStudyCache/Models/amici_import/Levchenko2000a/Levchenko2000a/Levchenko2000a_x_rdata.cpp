#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Levchenko2000a(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = MAPKP;
    x_rdata[1] = MEKP;
    x_rdata[2] = RAFK;
    x_rdata[3] = RAFP;
    x_rdata[4] = K_1_0;
    x_rdata[5] = K_1_1;
    x_rdata[6] = K_1_2;
    x_rdata[7] = K_2_0;
    x_rdata[8] = K_2_1;
    x_rdata[9] = K_2_2;
    x_rdata[10] = K_3_0;
    x_rdata[11] = K_3_1;
    x_rdata[12] = K_K_1_0_2_2;
    x_rdata[13] = K_K_1_1_2_2;
    x_rdata[14] = K_K_2_0_3_1;
    x_rdata[15] = K_K_2_1_3_1;
    x_rdata[16] = K_MAPKP_1_1;
    x_rdata[17] = K_MAPKP_1_2;
    x_rdata[18] = K_MEKP_2_1;
    x_rdata[19] = K_MEKP_2_2;
    x_rdata[20] = K_RAFK_3_0;
    x_rdata[21] = K_RAFP_3_1;
    x_rdata[22] = S_m1_m1_m1;
    x_rdata[23] = S_m1_m1_0;
    x_rdata[24] = S_m1_m1_1;
    x_rdata[25] = S_m1_0_m1;
    x_rdata[26] = S_m1_0_0;
    x_rdata[27] = S_m1_0_1;
    x_rdata[28] = S_m1_1_m1;
    x_rdata[29] = S_m1_1_0;
    x_rdata[30] = S_m1_1_1;
    x_rdata[31] = S_m1_2_m1;
    x_rdata[32] = S_m1_2_0;
    x_rdata[33] = S_m1_2_1;
    x_rdata[34] = S_0_m1_m1;
    x_rdata[35] = S_0_m1_0;
    x_rdata[36] = S_0_m1_1;
    x_rdata[37] = S_0_0_m1;
    x_rdata[38] = S_0_0_0;
    x_rdata[39] = S_0_0_1;
    x_rdata[40] = S_0_1_m1;
    x_rdata[41] = S_0_1_0;
    x_rdata[42] = S_0_1_1;
    x_rdata[43] = S_0_2_m1;
    x_rdata[44] = S_0_2_0;
    x_rdata[45] = S_0_2_1;
    x_rdata[46] = S_1_m1_m1;
    x_rdata[47] = S_1_m1_0;
    x_rdata[48] = S_1_m1_1;
    x_rdata[49] = S_1_0_m1;
    x_rdata[50] = S_1_0_0;
    x_rdata[51] = S_1_0_1;
    x_rdata[52] = S_1_1_m1;
    x_rdata[53] = S_1_1_0;
    x_rdata[54] = S_1_1_1;
    x_rdata[55] = S_1_2_m1;
    x_rdata[56] = S_1_2_0;
    x_rdata[57] = S_1_2_1;
    x_rdata[58] = S_2_m1_m1;
    x_rdata[59] = S_2_m1_0;
    x_rdata[60] = S_2_m1_1;
    x_rdata[61] = S_2_0_m1;
    x_rdata[62] = S_2_0_0;
    x_rdata[63] = S_2_0_1;
    x_rdata[64] = S_2_1_m1;
    x_rdata[65] = S_2_1_0;
    x_rdata[66] = S_2_1_1;
    x_rdata[67] = S_2_2_m1;
    x_rdata[68] = S_2_2_0;
    x_rdata[69] = S_2_2_1;
    x_rdata[70] = S_RAFK_m1_m1_0;
    x_rdata[71] = S_RAFK_m1_0_0;
    x_rdata[72] = S_RAFK_m1_1_0;
    x_rdata[73] = S_RAFK_m1_2_0;
    x_rdata[74] = S_RAFK_0_m1_0;
    x_rdata[75] = S_RAFK_0_0_0;
    x_rdata[76] = S_RAFK_0_1_0;
    x_rdata[77] = S_RAFK_0_2_0;
    x_rdata[78] = S_RAFK_1_m1_0;
    x_rdata[79] = S_RAFK_1_0_0;
    x_rdata[80] = S_RAFK_1_1_0;
    x_rdata[81] = S_RAFK_1_2_0;
    x_rdata[82] = S_RAFK_2_m1_0;
    x_rdata[83] = S_RAFK_2_0_0;
    x_rdata[84] = S_RAFK_2_1_0;
    x_rdata[85] = S_RAFK_2_2_0;
}