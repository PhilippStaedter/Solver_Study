#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model5_levchenko2(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 8:
            dwdp[32] = 1.0*C2*RAFp;
            dwdp[35] = 1.0*C2*RAFp;
            dwdp[40] = 1.0*C6*RAFp;
            dwdp[43] = 1.0*C6*RAFp;
            dwdp[45] = 1.0*C9*RAFp;
            dwdp[46] = 1.0*C9*RAFp;
            break;
        case 9:
            dwdp[42] = C7;
            dwdp[45] = C7;
            break;
        case 10:
            dwdp[31] = C2;
            dwdp[32] = C2;
            dwdp[37] = C6;
            dwdp[39] = C9;
            dwdp[40] = C6;
            dwdp[46] = C9;
            dwdp[53] = C2 + C6 + C9;
            break;
        case 11:
            dwdp[31] = C4;
            dwdp[33] = C6;
            dwdp[35] = C7;
            dwdp[36] = C4;
            dwdp[40] = C6;
            dwdp[42] = C7;
            dwdp[50] = C4 + C6 + C7;
            break;
        case 12:
            dwdp[31] = C3;
            dwdp[34] = C3;
            dwdp[37] = C7;
            dwdp[39] = C8;
            dwdp[42] = C7;
            dwdp[44] = C8;
            dwdp[51] = C3 + C7 + C8;
            break;
        case 13:
            dwdp[31] = C5;
            dwdp[33] = C9;
            dwdp[35] = C8;
            dwdp[38] = C5;
            dwdp[44] = C8;
            dwdp[46] = C9;
            dwdp[48] = C5 + C8 + C9;
            break;
        case 14:
            dwdp[30] = 1.0*C1*MEK;
            dwdp[33] = 1.0*C1*MEK;
            dwdp[36] = 1.0*C4*MEK;
            dwdp[38] = 1.0*C5*MEK;
            dwdp[41] = 1.0*C4*MEK;
            dwdp[47] = 1.0*C5*MEK;
            dwdp[52] = 1.0*MEK*(C1 + C4 + C5);
            break;
        case 15:
            dwdp[30] = 1.0*C1*MAPK;
            dwdp[32] = 1.0*C2*MAPK;
            dwdp[34] = 1.0*C3*MAPK;
            dwdp[37] = 1.0*C1*MAPK;
            dwdp[41] = 1.0*C2*MAPK;
            dwdp[43] = 1.0*C3*MAPK;
            dwdp[49] = 1.0*MAPK*(C1 + C2 + C3);
            break;
        case 16:
            dwdp[0] = 1.0*RAF*RAFK;
            break;
        case 17:
            dwdp[1] = 1.0*MEKPH*MEKp;
            break;
        case 18:
            dwdp[2] = 1.0*MEKpMEKPH;
            break;
        case 19:
            dwdp[3] = 1.0*MEKpMEKPH;
            break;
        case 20:
            dwdp[4] = 1.0*MEKp*RAFp;
            break;
        case 21:
            dwdp[5] = 1.0*MEKpRAFp;
            break;
        case 22:
            dwdp[6] = 1.0*MEKpRAFp;
            break;
        case 23:
            dwdp[7] = 1.0*MEKPH*MEKpp;
            break;
        case 24:
            dwdp[8] = 1.0*MEKppMEKPH;
            break;
        case 25:
            dwdp[9] = 1.0*MEKppMEKPH;
            break;
        case 26:
            dwdp[10] = 1.0*MAPK*MEKpp;
            break;
        case 27:
            dwdp[11] = 1.0*RAFRAFK;
            break;
        case 28:
            dwdp[12] = 1.0*MAPKMEKpp;
            break;
        case 29:
            dwdp[13] = 1.0*MAPKMEKpp;
            break;
        case 30:
            dwdp[14] = 1.0*MAPKPH*MAPKp;
            break;
        case 31:
            dwdp[15] = 1.0*MAPKpMAPKPH;
            break;
        case 32:
            dwdp[16] = 1.0*MAPKpMAPKPH;
            break;
        case 33:
            dwdp[17] = 1.0*MAPKp*MEKpp;
            break;
        case 34:
            dwdp[18] = 1.0*MAPKpMEKpp;
            break;
        case 35:
            dwdp[19] = 1.0*MAPKpMEKpp;
            break;
        case 36:
            dwdp[20] = 1.0*MAPKPH*MAPKpp;
            break;
        case 37:
            dwdp[21] = 1.0*MAPKppMAPKPH;
            break;
        case 38:
            dwdp[22] = 1.0*RAFRAFK;
            break;
        case 39:
            dwdp[23] = 1.0*MAPKppMAPKPH;
            break;
        case 40:
            dwdp[24] = 1.0*RAFPH*RAFp;
            break;
        case 41:
            dwdp[25] = 1.0*RAFpRAFPH;
            break;
        case 42:
            dwdp[26] = 1.0*RAFpRAFPH;
            break;
        case 43:
            dwdp[27] = 1.0*MEK*RAFp;
            break;
        case 44:
            dwdp[28] = 1.0*MEKRAFp;
            break;
        case 45:
            dwdp[29] = 1.0*MEKRAFp;
            break;
    }
}