#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model1_sarma1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 1:
            dwdp[0] = 1.0*species_1*species_2;
            break;
        case 2:
            dwdp[0] = -1.0*species_3;
            break;
        case 3:
            dwdp[1] = 1.0*species_7*species_8;
            break;
        case 4:
            dwdp[1] = -1.0*species_9;
            break;
        case 5:
            dwdp[2] = 1.0*species_9;
            break;
        case 6:
            dwdp[3] = 1.0*species_10*species_8;
            break;
        case 7:
            dwdp[3] = -1.0*species_11;
            break;
        case 8:
            dwdp[4] = 1.0*species_11;
            break;
        case 9:
            dwdp[5] = 1.0*species_13*species_2;
            break;
        case 10:
            dwdp[5] = -1.0*species_12;
            break;
        case 11:
            dwdp[6] = 1.0*species_12;
            break;
        case 12:
            dwdp[7] = 1.0*species_10*species_13;
            break;
        case 13:
            dwdp[7] = -1.0*species_14;
            break;
        case 14:
            dwdp[8] = 1.0*species_14;
            break;
        case 15:
            dwdp[9] = 1.0*species_15;
            break;
        case 16:
            dwdp[9] = -1.0*species_13*species_7;
            break;
        case 17:
            dwdp[10] = 1.0*species_16*species_18;
            break;
        case 18:
            dwdp[10] = -1.0*species_17;
            break;
        case 19:
            dwdp[11] = 1.0*species_3;
            break;
        case 20:
            dwdp[12] = 1.0*species_17;
            break;
        case 21:
            dwdp[13] = 1.0*species_20*species_8;
            break;
        case 22:
            dwdp[13] = -1.0*species_19;
            break;
        case 23:
            dwdp[14] = 1.0*species_19;
            break;
        case 24:
            dwdp[15] = 1.0*species_2*species_20;
            break;
        case 25:
            dwdp[15] = -1.0*species_24;
            break;
        case 26:
            dwdp[16] = 1.0*species_27;
            break;
        case 27:
            dwdp[16] = -1.0*species_16*species_20;
            break;
        case 28:
            dwdp[17] = 1.0*species_24;
            break;
        case 29:
            dwdp[18] = 1.0*species_10*species_20;
            break;
        case 30:
            dwdp[18] = -1.0*species_25;
            break;
        case 31:
            dwdp[19] = 1.0*species_25;
            break;
        case 32:
            dwdp[20] = 1.0*species_26;
            break;
        case 33:
            dwdp[20] = -1.0*species_20*species_7;
            break;
        case 34:
            dwdp[21] = 1.0*species_2*species_4;
            break;
        case 35:
            dwdp[21] = -1.0*species_5;
            break;
        case 36:
            dwdp[22] = 1.0*species_5;
            break;
        case 37:
            dwdp[23] = 1.0*species_13*species_6;
            break;
        case 38:
            dwdp[23] = -1.0*species_21;
            break;
        case 39:
            dwdp[24] = 1.0*species_21;
            break;
        case 40:
            dwdp[25] = 1.0*species_13*species_4;
            break;
        case 41:
            dwdp[25] = -1.0*species_22;
            break;
        case 42:
            dwdp[26] = 1.0*species_22;
            break;
        case 43:
            dwdp[27] = 1.0*species_23;
            break;
        case 44:
            dwdp[27] = -1.0*species_1*species_13;
            break;
    }
}