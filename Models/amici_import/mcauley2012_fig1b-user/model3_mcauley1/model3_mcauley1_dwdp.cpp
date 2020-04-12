#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model3_mcauley1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = reaction_1_k1*species_1;
            break;
        case 1:
            dwdp[0] = multiplier*species_1;
            break;
        case 2:
            dwdp[1] = 1.0/(pow(reaction_10_BCRt/species_7, reaction_10_BS) + 1);
            break;
        case 3:
            dwdp[1] = -reaction_10_BCRmax*reaction_10_BS*pow(reaction_10_BCRt/species_7, reaction_10_BS)/(reaction_10_BCRt*pow(pow(reaction_10_BCRt/species_7, reaction_10_BS) + 1, 2));
            break;
        case 4:
            dwdp[1] = -reaction_10_BCRmax*pow(reaction_10_BCRt/species_7, reaction_10_BS)*log(reaction_10_BCRt/species_7)/pow(pow(reaction_10_BCRt/species_7, reaction_10_BS) + 1, 2);
            break;
        case 5:
            dwdp[2] = 1.0/(pow(species_7/reaction_11_HCSt, reaction_11_HS) + 1);
            break;
        case 6:
            dwdp[2] = reaction_11_HCSmax*reaction_11_HS*pow(species_7/reaction_11_HCSt, reaction_11_HS)/(reaction_11_HCSt*pow(pow(species_7/reaction_11_HCSt, reaction_11_HS) + 1, 2));
            break;
        case 7:
            dwdp[2] = -reaction_11_HCSmax*pow(species_7/reaction_11_HCSt, reaction_11_HS)*log(species_7/reaction_11_HCSt)/pow(pow(species_7/reaction_11_HCSt, reaction_11_HS) + 1, 2);
            break;
        case 8:
            dwdp[3] = species_14*species_7;
            break;
        case 9:
            dwdp[4] = species_13*species_15;
            break;
        case 10:
            dwdp[5] = species_11;
            break;
        case 11:
            dwdp[6] = species_7;
            break;
        case 12:
            dwdp[7] = species_19/species_7;
            break;
        case 13:
            dwdp[8] = species_18;
            break;
        case 14:
            dwdp[9] = species_17;
            break;
        case 15:
            dwdp[10] = species_17*species_22;
            break;
        case 16:
            dwdp[11] = 1.0/(pow(species_2/reaction_2_ICt, reaction_2_IS) + 1);
            break;
        case 17:
            dwdp[11] = reaction_2_ICSmax*reaction_2_IS*pow(species_2/reaction_2_ICt, reaction_2_IS)/(reaction_2_ICt*pow(pow(species_2/reaction_2_ICt, reaction_2_IS) + 1, 2));
            break;
        case 18:
            dwdp[11] = -reaction_2_ICSmax*pow(species_2/reaction_2_ICt, reaction_2_IS)*log(species_2/reaction_2_ICt)/pow(pow(species_2/reaction_2_ICt, reaction_2_IS) + 1, 2);
            break;
        case 19:
            dwdp[12] = species_21;
            break;
        case 20:
            dwdp[13] = species_21*species_24;
            break;
        case 21:
            dwdp[14] = species_18*species_23;
            break;
        case 22:
            dwdp[15] = species_23;
            break;
        case 23:
            dwdp[16] = species_23*species_25;
            break;
        case 24:
            dwdp[17] = species_23;
            break;
        case 25:
            dwdp[18] = species_26/species_11;
            break;
        case 26:
            dwdp[19] = species_25;
            break;
        case 27:
            dwdp[20] = species_11*species_14;
            break;
        case 28:
            dwdp[21] = species_15*species_28;
            break;
        case 29:
            dwdp[22] = species_4;
            break;
        case 30:
            dwdp[23] = species_11;
            break;
        case 31:
            dwdp[24] = species_10*species_11*species_31;
            break;
        case 32:
            dwdp[25] = 1.0/(pow(species_11/reaction_32_PPCt, reaction_32_PCSS) + 1);
            break;
        case 33:
            dwdp[25] = reaction_32_PCSS*reaction_32_PCSmax*pow(species_11/reaction_32_PPCt, reaction_32_PCSS)/(reaction_32_PPCt*pow(pow(species_11/reaction_32_PPCt, reaction_32_PCSS) + 1, 2));
            break;
        case 34:
            dwdp[25] = -reaction_32_PCSmax*pow(species_11/reaction_32_PPCt, reaction_32_PCSS)*log(species_11/reaction_32_PPCt)/pow(pow(species_11/reaction_32_PPCt, reaction_32_PCSS) + 1, 2);
            break;
        case 35:
            dwdp[26] = species_30*species_33;
            break;
        case 36:
            dwdp[27] = species_30*species_33;
            break;
        case 37:
            dwdp[28] = species_30*species_34;
            break;
        case 38:
            dwdp[29] = species_5;
            break;
        case 39:
            dwdp[30] = species_5;
            break;
        case 40:
            dwdp[31] = species_7/species_4;
            break;
        case 41:
            dwdp[32] = species_2*species_5;
            break;
        case 42:
            dwdp[33] = species_2*species_5;
            break;
        case 43:
            dwdp[34] = species_11;
            break;
    }
}