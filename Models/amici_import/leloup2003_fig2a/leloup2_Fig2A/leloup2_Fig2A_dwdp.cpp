#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_leloup2_Fig2A(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*pow(reaction_0_K, reaction_0_m)*reaction_0_vsb*log(reaction_0_K)/(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m)) + 1.0*pow(reaction_0_K, reaction_0_m)*reaction_0_vsb*(-pow(reaction_0_K, reaction_0_m)*log(reaction_0_K) - pow(species_3, reaction_0_m)*log(species_3))/pow(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m), 2);
            break;
        case 1:
            dwdp[0] = -1.0*pow(reaction_0_K, 2*reaction_0_m)*reaction_0_m*reaction_0_vsb/(reaction_0_K*pow(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m), 2)) + 1.0*pow(reaction_0_K, reaction_0_m)*reaction_0_m*reaction_0_vsb/(reaction_0_K*(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m)));
            break;
        case 2:
            dwdp[0] = 1.0*pow(reaction_0_K, reaction_0_m)/(pow(reaction_0_K, reaction_0_m) + pow(species_3, reaction_0_m));
            break;
        case 3:
            dwdp[1] = 1.0*species_0;
            break;
        case 4:
            dwdp[2] = 1.0*species_7;
            break;
        case 5:
            dwdp[3] = -1.0*reaction_11_V*species_9/pow(reaction_11_Km + species_9, 2);
            break;
        case 6:
            dwdp[3] = 1.0*species_9/(reaction_11_Km + species_9);
            break;
        case 7:
            dwdp[4] = -1.0*reaction_12_V*species_8/pow(reaction_12_Km + species_8, 2);
            break;
        case 8:
            dwdp[4] = 1.0*species_8/(reaction_12_Km + species_8);
            break;
        case 9:
            dwdp[5] = -1.0*species_10;
            break;
        case 10:
            dwdp[5] = 1.0*species_4*species_8;
            break;
        case 11:
            dwdp[6] = -1.0*reaction_14_V*species_10/pow(reaction_14_Km + species_10, 2);
            break;
        case 12:
            dwdp[6] = 1.0*species_10/(reaction_14_Km + species_10);
            break;
        case 13:
            dwdp[7] = -1.0*reaction_15_V*species_11/pow(reaction_15_Km + species_11, 2);
            break;
        case 14:
            dwdp[7] = 1.0*species_11/(reaction_15_Km + species_11);
            break;
        case 15:
            dwdp[8] = -1.0*species_12;
            break;
        case 16:
            dwdp[8] = 1.0*species_10;
            break;
        case 17:
            dwdp[9] = 1.0*species_14;
            break;
        case 18:
            dwdp[10] = 1.0*species_2;
            break;
        case 19:
            dwdp[11] = 1.0*species_13;
            break;
        case 20:
            dwdp[12] = 1.0*species_0;
            break;
        case 21:
            dwdp[13] = -1.0*pow(reaction_20_K, reaction_20_n)*reaction_20_Vs*reaction_20_n*pow(species_3, reaction_20_n)/(reaction_20_K*pow(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n), 2));
            break;
        case 22:
            dwdp[13] = 1.0*reaction_20_Vs*pow(species_3, reaction_20_n)*log(species_3)/(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n)) + 1.0*reaction_20_Vs*pow(species_3, reaction_20_n)*(-pow(reaction_20_K, reaction_20_n)*log(reaction_20_K) - pow(species_3, reaction_20_n)*log(species_3))/pow(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n), 2);
            break;
        case 23:
            dwdp[13] = 1.0*pow(species_3, reaction_20_n)/(pow(reaction_20_K, reaction_20_n) + pow(species_3, reaction_20_n));
            break;
        case 24:
            dwdp[14] = -1.0*reaction_21_V*species_12/pow(reaction_21_Km + species_12, 2);
            break;
        case 25:
            dwdp[14] = 1.0*species_12/(reaction_21_Km + species_12);
            break;
        case 26:
            dwdp[15] = 1.0*species_7;
            break;
        case 27:
            dwdp[16] = -1.0*species_15;
            break;
        case 28:
            dwdp[16] = 1.0*species_12*species_3;
            break;
        case 29:
            dwdp[17] = -1.0*reaction_24_V*species_0/pow(reaction_24_Km + species_0, 2);
            break;
        case 30:
            dwdp[17] = 1.0*species_0/(reaction_24_Km + species_0);
            break;
        case 31:
            dwdp[18] = -1.0*reaction_25_V*species_5/pow(reaction_25_Km + species_5, 2);
            break;
        case 32:
            dwdp[18] = 1.0*species_5/(reaction_25_Km + species_5);
            break;
        case 33:
            dwdp[19] = -1.0*reaction_26_V*species_7/pow(reaction_26_Km + species_7, 2);
            break;
        case 34:
            dwdp[19] = 1.0*species_7/(reaction_26_Km + species_7);
            break;
        case 35:
            dwdp[20] = 1.0*species_8;
            break;
        case 36:
            dwdp[21] = 1.0*species_4;
            break;
        case 37:
            dwdp[22] = 1.0*species_9;
            break;
        case 38:
            dwdp[23] = -1.0*reaction_3_V*species_1/pow(reaction_3_Km + species_1, 2);
            break;
        case 39:
            dwdp[23] = 1.0*species_1/(reaction_3_Km + species_1);
            break;
        case 40:
            dwdp[24] = 1.0*species_6;
            break;
        case 41:
            dwdp[25] = 1.0*species_11;
            break;
        case 42:
            dwdp[26] = 1.0*species_10;
            break;
        case 43:
            dwdp[27] = -1.0*reaction_33_V*species_14/pow(reaction_33_Km + species_14, 2);
            break;
        case 44:
            dwdp[27] = 1.0*species_14/(reaction_33_Km + species_14);
            break;
        case 45:
            dwdp[28] = 1.0*species_1;
            break;
        case 46:
            dwdp[29] = -1.0*reaction_35_V*species_2/pow(reaction_35_Km + species_2, 2);
            break;
        case 47:
            dwdp[29] = 1.0*species_2/(reaction_35_Km + species_2);
            break;
        case 48:
            dwdp[30] = -1.0*reaction_36_V*species_3/pow(reaction_36_Km + species_3, 2);
            break;
        case 49:
            dwdp[30] = 1.0*species_3/(reaction_36_Km + species_3);
            break;
        case 50:
            dwdp[31] = -1.0*reaction_37_V*species_13/pow(reaction_37_Km + species_13, 2);
            break;
        case 51:
            dwdp[31] = 1.0*species_13/(reaction_37_Km + species_13);
            break;
        case 52:
            dwdp[32] = 1.0*species_15;
            break;
        case 53:
            dwdp[33] = -1.0*reaction_39_V*species_15/pow(reaction_39_Km + species_15, 2);
            break;
        case 54:
            dwdp[33] = 1.0*species_15/(reaction_39_Km + species_15);
            break;
        case 55:
            dwdp[34] = -1.0*species_3;
            break;
        case 56:
            dwdp[34] = 1.0*species_1;
            break;
        case 57:
            dwdp[35] = 1.0*species_3;
            break;
        case 58:
            dwdp[36] = -1.0*reaction_41_V*species_2/pow(reaction_41_Km + species_2, 2);
            break;
        case 59:
            dwdp[36] = 1.0*species_2/(reaction_41_Km + species_2);
            break;
        case 60:
            dwdp[37] = -1.0*reaction_42_V*species_13/pow(reaction_42_Km + species_13, 2);
            break;
        case 61:
            dwdp[37] = 1.0*species_13/(reaction_42_Km + species_13);
            break;
        case 62:
            dwdp[38] = -1.0*reaction_43_V*species_6/pow(reaction_43_Km + species_6, 2);
            break;
        case 63:
            dwdp[38] = 1.0*species_6/(reaction_43_Km + species_6);
            break;
        case 64:
            dwdp[39] = -1.0*reaction_44_V*species_9/pow(reaction_44_Km + species_9, 2);
            break;
        case 65:
            dwdp[39] = 1.0*species_9/(reaction_44_Km + species_9);
            break;
        case 66:
            dwdp[40] = -1.0*reaction_45_V*species_14/pow(reaction_45_Km + species_14, 2);
            break;
        case 67:
            dwdp[40] = 1.0*species_14/(reaction_45_Km + species_14);
            break;
        case 68:
            dwdp[41] = 1.0*species_12;
            break;
        case 69:
            dwdp[42] = -1.0*reaction_47_V*species_11/pow(reaction_47_Km + species_11, 2);
            break;
        case 70:
            dwdp[42] = 1.0*species_11/(reaction_47_Km + species_11);
            break;
        case 71:
            dwdp[43] = 1.0*species_5;
            break;
        case 72:
            dwdp[44] = 1.0*species_5;
            break;
        case 73:
            dwdp[45] = -1.0*reaction_7_V*species_4/pow(reaction_7_Km + species_4, 2);
            break;
        case 74:
            dwdp[45] = 1.0*species_4/(reaction_7_Km + species_4);
            break;
        case 75:
            dwdp[46] = -1.0*reaction_8_V*species_6/pow(reaction_8_Km + species_6, 2);
            break;
        case 76:
            dwdp[46] = 1.0*species_6/(reaction_8_Km + species_6);
            break;
        case 77:
            dwdp[47] = -1.0*pow(reaction_9_K, reaction_9_n)*reaction_9_Vs*reaction_9_n*pow(species_3, reaction_9_n)/(reaction_9_K*pow(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n), 2));
            break;
        case 78:
            dwdp[47] = 1.0*reaction_9_Vs*pow(species_3, reaction_9_n)*log(species_3)/(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n)) + 1.0*reaction_9_Vs*pow(species_3, reaction_9_n)*(-pow(reaction_9_K, reaction_9_n)*log(reaction_9_K) - pow(species_3, reaction_9_n)*log(species_3))/pow(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n), 2);
            break;
        case 79:
            dwdp[47] = 1.0*pow(species_3, reaction_9_n)/(pow(reaction_9_K, reaction_9_n) + pow(species_3, reaction_9_n));
            break;
    }
}