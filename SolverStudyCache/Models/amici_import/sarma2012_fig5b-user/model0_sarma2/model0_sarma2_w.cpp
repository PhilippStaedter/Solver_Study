#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_sarma2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*reaction_1_k1*species_1*species_2 - 1.0*reaction_1_k2*species_3;
    w[1] = 1.0*reaction_10_k1*species_7*species_8 - 1.0*reaction_10_k2*species_9;
    w[2] = 1.0*reaction_11_k1*species_9;
    w[3] = 1.0*reaction_12_k1*species_10*species_8 - 1.0*reaction_12_k2*species_11;
    w[4] = 1.0*reaction_13_k1*species_11;
    w[5] = 1.0*reaction_14_k1*species_13*species_2 - 1.0*reaction_14_k2*species_12;
    w[6] = 1.0*reaction_15_k1*species_12;
    w[7] = 1.0*reaction_16_k1*species_10*species_13 - 1.0*reaction_16_k2*species_14;
    w[8] = 1.0*reaction_17_k1*species_14;
    w[9] = 1.0*reaction_18_k1*species_15 - 1.0*reaction_18_k2*species_13*species_7;
    w[10] = 1.0*reaction_19_k1*species_16*species_18 - 1.0*reaction_19_k2*species_17;
    w[11] = 1.0*reaction_2_k1*species_3;
    w[12] = 1.0*reaction_20_k1*species_17;
    w[13] = 1.0*reaction_21_k1*species_20*species_8 - 1.0*reaction_21_k2*species_19;
    w[14] = 1.0*reaction_22_k1*species_19;
    w[15] = 1.0*reaction_23_k1*species_2*species_20 - 1.0*reaction_23_k2*species_24;
    w[16] = 1.0*reaction_24_k1*species_27 - 1.0*reaction_24_k2*species_16*species_20;
    w[17] = 1.0*reaction_25_k1*species_24;
    w[18] = 1.0*reaction_26_k1*species_10*species_20 - 1.0*reaction_26_k2*species_25;
    w[19] = 1.0*reaction_27_k1*species_25;
    w[20] = 1.0*reaction_28_k1*species_26 - 1.0*reaction_28_k2*species_20*species_7;
    w[21] = 1.0*reaction_3_k1*species_2*species_4 - 1.0*reaction_3_k2*species_5;
    w[22] = 1.0*reaction_4_k1*species_5;
    w[23] = 1.0*reaction_5_k1*species_13*species_6 - 1.0*reaction_5_k2*species_21;
    w[24] = 1.0*reaction_6_k1*species_21;
    w[25] = 1.0*reaction_7_k1*species_13*species_4 - 1.0*reaction_7_k2*species_22;
    w[26] = 1.0*reaction_8_k1*species_22;
    w[27] = 1.0*reaction_9_k1*species_23 - 1.0*reaction_9_k2*species_1*species_13;
}