#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_kofahl(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = Bar1aex*alpha;
            break;
        case 1:
            dwdp[1] = Gbc*complexC;
            break;
        case 2:
            dwdp[2] = complexD;
            break;
        case 3:
            dwdp[3] = Ste11*Ste5;
            break;
        case 4:
            dwdp[4] = complexA;
            break;
        case 5:
            dwdp[5] = Fus3*Ste7;
            break;
        case 6:
            dwdp[6] = complexB;
            break;
        case 7:
            dwdp[7] = complexA*complexB;
            break;
        case 8:
            dwdp[8] = complexC;
            break;
        case 9:
            dwdp[9] = Ste20*complexD;
            break;
        case 10:
            dwdp[10] = complexE;
            break;
        case 11:
            dwdp[11] = Ste2*alpha;
            break;
        case 12:
            dwdp[12] = complexE;
            break;
        case 13:
            dwdp[13] = complexE;
            break;
        case 14:
            dwdp[14] = complexF;
            break;
        case 15:
            dwdp[15] = complexF;
            break;
        case 16:
            dwdp[16] = complexG;
            break;
        case 17:
            dwdp[17] = complexG;
            break;
        case 18:
            dwdp[18] = complexH;
            break;
        case 19:
            dwdp[19] = complexH;
            break;
        case 20:
            dwdp[20] = complexI;
            break;
        case 21:
            dwdp[21] = Fus3*complexL;
            break;
        case 22:
            dwdp[22] = Ste2a;
            break;
        case 23:
            dwdp[23] = complexK;
            break;
        case 24:
            dwdp[24] = complexK;
            break;
        case 25:
            dwdp[25] = complexL;
            break;
        case 26:
            dwdp[26] = Fus3PP;
            break;
        case 27:
            dwdp[27] = Fus3PP*Ste12;
            break;
        case 28:
            dwdp[28] = Ste12a;
            break;
        case 29:
            dwdp[29] = Bar1*Ste12a;
            break;
        case 30:
            dwdp[30] = Bar1a;
            break;
        case 31:
            dwdp[31] = Bar1a;
            break;
        case 32:
            dwdp[32] = Far1*pow(Fus3PP, 2)/(pow(Fus3PP, 2) + 10000);
            break;
        case 33:
            dwdp[33] = Ste2a;
            break;
        case 34:
            dwdp[34] = Far1PP;
            break;
        case 35:
            dwdp[35] = Cdc28*Far1;
            break;
        case 36:
            dwdp[36] = Far1PP*Gbc;
            break;
        case 37:
            dwdp[37] = complexM;
            break;
        case 38:
            dwdp[38] = complexN;
            break;
        case 39:
            dwdp[39] = Cdc28*Far1PP;
            break;
        case 40:
            dwdp[40] = pow(Fus3PP, 2)/(pow(Fus3PP, 2) + 16);
            break;
        case 41:
            dwdp[41] = Sst2;
            break;
        case 42:
            dwdp[42] = Ste2;
            break;
        case 43:
            dwdp[43] = Gabc*Ste2a;
            break;
        case 44:
            dwdp[44] = GaGTP;
            break;
        case 45:
            dwdp[45] = GaGTP*Sst2;
            break;
        case 46:
            dwdp[46] = GaGDP*Gbc;
            break;
    }
}