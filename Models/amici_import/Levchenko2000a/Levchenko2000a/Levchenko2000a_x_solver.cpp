#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Levchenko2000a(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = MAPKP;
    x_solver[1] = MEKP;
    x_solver[2] = RAFK;
    x_solver[3] = RAFP;
    x_solver[4] = K_1_0;
    x_solver[5] = K_1_1;
    x_solver[6] = K_1_2;
    x_solver[7] = K_2_0;
    x_solver[8] = K_2_1;
    x_solver[9] = K_2_2;
    x_solver[10] = K_3_0;
    x_solver[11] = K_3_1;
    x_solver[12] = K_K_1_0_2_2;
    x_solver[13] = K_K_1_1_2_2;
    x_solver[14] = K_K_2_0_3_1;
    x_solver[15] = K_K_2_1_3_1;
    x_solver[16] = K_MAPKP_1_1;
    x_solver[17] = K_MAPKP_1_2;
    x_solver[18] = K_MEKP_2_1;
    x_solver[19] = K_MEKP_2_2;
    x_solver[20] = K_RAFK_3_0;
    x_solver[21] = K_RAFP_3_1;
    x_solver[22] = S_m1_m1_m1;
    x_solver[23] = S_m1_m1_0;
    x_solver[24] = S_m1_m1_1;
    x_solver[25] = S_m1_0_m1;
    x_solver[26] = S_m1_0_0;
    x_solver[27] = S_m1_0_1;
    x_solver[28] = S_m1_1_m1;
    x_solver[29] = S_m1_1_0;
    x_solver[30] = S_m1_1_1;
    x_solver[31] = S_m1_2_m1;
    x_solver[32] = S_m1_2_0;
    x_solver[33] = S_m1_2_1;
    x_solver[34] = S_0_m1_m1;
    x_solver[35] = S_0_m1_0;
    x_solver[36] = S_0_m1_1;
    x_solver[37] = S_0_0_m1;
    x_solver[38] = S_0_0_0;
    x_solver[39] = S_0_0_1;
    x_solver[40] = S_0_1_m1;
    x_solver[41] = S_0_1_0;
    x_solver[42] = S_0_1_1;
    x_solver[43] = S_0_2_m1;
    x_solver[44] = S_0_2_0;
    x_solver[45] = S_0_2_1;
    x_solver[46] = S_1_m1_m1;
    x_solver[47] = S_1_m1_0;
    x_solver[48] = S_1_m1_1;
    x_solver[49] = S_1_0_m1;
    x_solver[50] = S_1_0_0;
    x_solver[51] = S_1_0_1;
    x_solver[52] = S_1_1_m1;
    x_solver[53] = S_1_1_0;
    x_solver[54] = S_1_1_1;
    x_solver[55] = S_1_2_m1;
    x_solver[56] = S_1_2_0;
    x_solver[57] = S_1_2_1;
    x_solver[58] = S_2_m1_m1;
    x_solver[59] = S_2_m1_0;
    x_solver[60] = S_2_m1_1;
    x_solver[61] = S_2_0_m1;
    x_solver[62] = S_2_0_0;
    x_solver[63] = S_2_0_1;
    x_solver[64] = S_2_1_m1;
    x_solver[65] = S_2_1_0;
    x_solver[66] = S_2_1_1;
    x_solver[67] = S_2_2_m1;
    x_solver[68] = S_2_2_0;
    x_solver[69] = S_2_2_1;
    x_solver[70] = S_RAFK_m1_m1_0;
    x_solver[71] = S_RAFK_m1_0_0;
    x_solver[72] = S_RAFK_m1_1_0;
    x_solver[73] = S_RAFK_m1_2_0;
    x_solver[74] = S_RAFK_0_m1_0;
    x_solver[75] = S_RAFK_0_0_0;
    x_solver[76] = S_RAFK_0_1_0;
    x_solver[77] = S_RAFK_0_2_0;
    x_solver[78] = S_RAFK_1_m1_0;
    x_solver[79] = S_RAFK_1_0_0;
    x_solver[80] = S_RAFK_1_1_0;
    x_solver[81] = S_RAFK_1_2_0;
    x_solver[82] = S_RAFK_2_m1_0;
    x_solver[83] = S_RAFK_2_0_0;
    x_solver[84] = S_RAFK_2_1_0;
    x_solver[85] = S_RAFK_2_2_0;
}