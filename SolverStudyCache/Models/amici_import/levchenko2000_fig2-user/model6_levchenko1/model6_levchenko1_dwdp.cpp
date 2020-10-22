#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model6_levchenko1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*RAF*RAFK;
            break;
        case 1:
            dwdp[1] = 1.0*MEKPH*MEKp;
            break;
        case 2:
            dwdp[2] = 1.0*MEKpMEKPH;
            break;
        case 3:
            dwdp[3] = 1.0*MEKpMEKPH;
            break;
        case 4:
            dwdp[4] = 1.0*MEKp*RAFp;
            break;
        case 5:
            dwdp[5] = 1.0*MEKpRAFp;
            break;
        case 6:
            dwdp[6] = 1.0*MEKpRAFp;
            break;
        case 7:
            dwdp[7] = 1.0*MEKPH*MEKpp;
            break;
        case 8:
            dwdp[8] = 1.0*MEKppMEKPH;
            break;
        case 9:
            dwdp[9] = 1.0*MEKppMEKPH;
            break;
        case 10:
            dwdp[10] = 1.0*MAPK*MEKpp;
            break;
        case 11:
            dwdp[11] = 1.0*RAFRAFK;
            break;
        case 12:
            dwdp[12] = 1.0*MAPKMEKpp;
            break;
        case 13:
            dwdp[13] = 1.0*MAPKMEKpp;
            break;
        case 14:
            dwdp[14] = 1.0*MAPKPH*MAPKp;
            break;
        case 15:
            dwdp[15] = 1.0*MAPKpMAPKPH;
            break;
        case 16:
            dwdp[16] = 1.0*MAPKpMAPKPH;
            break;
        case 17:
            dwdp[17] = 1.0*MAPKp*MEKpp;
            break;
        case 18:
            dwdp[18] = 1.0*MAPKpMEKpp;
            break;
        case 19:
            dwdp[19] = 1.0*MAPKpMEKpp;
            break;
        case 20:
            dwdp[20] = 1.0*MAPKPH*MAPKpp;
            break;
        case 21:
            dwdp[21] = 1.0*MAPKppMAPKPH;
            break;
        case 22:
            dwdp[22] = 1.0*RAFRAFK;
            break;
        case 23:
            dwdp[23] = 1.0*MAPKppMAPKPH;
            break;
        case 24:
            dwdp[24] = 1.0*RAFPH*RAFp;
            break;
        case 25:
            dwdp[25] = 1.0*RAFpRAFPH;
            break;
        case 26:
            dwdp[26] = 1.0*RAFpRAFPH;
            break;
        case 27:
            dwdp[27] = 1.0*MEK*RAFp;
            break;
        case 28:
            dwdp[28] = 1.0*MEKRAFp;
            break;
        case 29:
            dwdp[29] = 1.0*MEKRAFp;
            break;
    }
}