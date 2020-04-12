#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_sarma1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*reaction_1_k1*species_2;
    dwdx[1] = -1.0*reaction_9_k2*species_13;
    dwdx[2] = 1.0*reaction_12_k1*species_8;
    dwdx[3] = 1.0*reaction_16_k1*species_13;
    dwdx[4] = 1.0*reaction_26_k1*species_20;
    dwdx[5] = -1.0*reaction_12_k2;
    dwdx[6] = 1.0*reaction_13_k1;
    dwdx[7] = -1.0*reaction_14_k2;
    dwdx[8] = 1.0*reaction_15_k1;
    dwdx[9] = 1.0*reaction_14_k1*species_2;
    dwdx[10] = 1.0*reaction_16_k1*species_10;
    dwdx[11] = -1.0*reaction_18_k2*species_7;
    dwdx[12] = 1.0*reaction_5_k1*species_6;
    dwdx[13] = 1.0*reaction_7_k1*species_4;
    dwdx[14] = -1.0*reaction_9_k2*species_1;
    dwdx[15] = -1.0*reaction_16_k2;
    dwdx[16] = 1.0*reaction_17_k1;
    dwdx[17] = 1.0*reaction_18_k1;
    dwdx[18] = 1.0*reaction_19_k1*species_18;
    dwdx[19] = -1.0*reaction_24_k2*species_20;
    dwdx[20] = -1.0*reaction_19_k2;
    dwdx[21] = 1.0*reaction_20_k1;
    dwdx[22] = 1.0*reaction_19_k1*species_16;
    dwdx[23] = -1.0*reaction_21_k2;
    dwdx[24] = 1.0*reaction_22_k1;
    dwdx[25] = 1.0*reaction_1_k1*species_1;
    dwdx[26] = 1.0*reaction_14_k1*species_13;
    dwdx[27] = 1.0*reaction_23_k1*species_20;
    dwdx[28] = 1.0*reaction_3_k1*species_4;
    dwdx[29] = 1.0*reaction_21_k1*species_8;
    dwdx[30] = 1.0*reaction_23_k1*species_2;
    dwdx[31] = -1.0*reaction_24_k2*species_16;
    dwdx[32] = 1.0*reaction_26_k1*species_10;
    dwdx[33] = -1.0*reaction_28_k2*species_7;
    dwdx[34] = -1.0*reaction_5_k2;
    dwdx[35] = 1.0*reaction_6_k1;
    dwdx[36] = -1.0*reaction_7_k2;
    dwdx[37] = 1.0*reaction_8_k1;
    dwdx[38] = 1.0*reaction_9_k1;
    dwdx[39] = -1.0*reaction_23_k2;
    dwdx[40] = 1.0*reaction_25_k1;
    dwdx[41] = -1.0*reaction_26_k2;
    dwdx[42] = 1.0*reaction_27_k1;
    dwdx[43] = 1.0*reaction_28_k1;
    dwdx[44] = 1.0*reaction_24_k1;
    dwdx[45] = -1.0*reaction_1_k2;
    dwdx[46] = 1.0*reaction_2_k1;
    dwdx[47] = 1.0*reaction_3_k1*species_2;
    dwdx[48] = 1.0*reaction_7_k1*species_13;
    dwdx[49] = -1.0*reaction_3_k2;
    dwdx[50] = 1.0*reaction_4_k1;
    dwdx[51] = 1.0*reaction_5_k1*species_13;
    dwdx[52] = 1.0*reaction_10_k1*species_8;
    dwdx[53] = -1.0*reaction_18_k2*species_13;
    dwdx[54] = -1.0*reaction_28_k2*species_20;
    dwdx[55] = 1.0*reaction_10_k1*species_7;
    dwdx[56] = 1.0*reaction_12_k1*species_10;
    dwdx[57] = 1.0*reaction_21_k1*species_20;
    dwdx[58] = -1.0*reaction_10_k2;
    dwdx[59] = 1.0*reaction_11_k1;
}