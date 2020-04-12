#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model0_schilling1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[11] = 1.0*EpoR*pJAK2;
            break;
        case 1:
            dwdp[20] = 1.0*ERK1*ppMEK2;
            dwdp[23] = 1.0*ERK1*ppMEK1;
            break;
        case 2:
            dwdp[21] = 1.0*ERK2*ppMEK2;
            dwdp[24] = 1.0*ERK2*ppMEK1;
            break;
        case 3:
            dwdp[29] = 1.0*ppERK1;
            dwdp[30] = 1.0*ppERK2;
            break;
        case 4:
            dwdp[13] = 1.0*MEK1*pRaf;
            break;
        case 5:
            dwdp[12] = 1.0*MEK2*pRaf;
            break;
        case 6:
            dwdp[16] = 1.0*ppMEK2;
            dwdp[17] = 1.0*ppMEK1;
            break;
        case 7:
            dwdp[0] = 1.0*Epo*JAK2;
            break;
        case 8:
            dwdp[22] = 1.0*SHP1*pEpoR;
            break;
        case 9:
            dwdp[1] = 1.0*Delay06_mSHP1;
            dwdp[2] = 1.0*Delay07_mSHP1;
            dwdp[3] = 1.0*Delay08_mSHP1;
            dwdp[33] = 1.0*mSHP1;
            dwdp[37] = 1.0*Delay01_mSHP1;
            dwdp[38] = 1.0*Delay02_mSHP1;
            dwdp[39] = 1.0*Delay03_mSHP1;
            dwdp[40] = 1.0*Delay04_mSHP1;
            dwdp[41] = 1.0*Delay05_mSHP1;
            break;
        case 10:
            dwdp[7] = 1.0*SOS*pEpoR;
            break;
        case 11:
            dwdp[25] = 1.0*pERK1*ppMEK2;
            dwdp[27] = 1.0*pERK1*ppMEK1;
            break;
        case 12:
            dwdp[26] = 1.0*pERK2*ppMEK2;
            dwdp[28] = 1.0*pERK2*ppMEK1;
            break;
        case 13:
            dwdp[31] = 1.0*pERK1;
            dwdp[32] = 1.0*pERK2;
            break;
        case 14:
            dwdp[15] = 1.0*pMEK1*pRaf;
            break;
        case 15:
            dwdp[14] = 1.0*pMEK2*pRaf;
            break;
        case 16:
            dwdp[18] = 1.0*pMEK2;
            dwdp[19] = 1.0*pMEK1;
            break;
        case 17:
            dwdp[4] = 1.0*actSHP1;
            break;
        case 18:
            dwdp[9] = 1.0*Raf*mSOS;
            break;
        case 19:
            dwdp[8] = 1.0*mSOS;
            break;
        case 20:
            dwdp[5] = 1.0*actSHP1*pEpoR;
            break;
        case 21:
            dwdp[6] = 1.0*actSHP1*pJAK2;
            break;
        case 22:
            dwdp[10] = 1.0*pRaf;
            break;
        case 23:
            dwdp[36] = 1.0*pSOS;
            break;
        case 24:
            dwdp[34] = 1.0*mSOS*ppERK1;
            dwdp[35] = 1.0*mSOS*ppERK2;
            break;
    }
}