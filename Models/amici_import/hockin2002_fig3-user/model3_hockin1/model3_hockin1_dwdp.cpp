#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_model3_hockin1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = -1.0*TF_VII;
            break;
        case 1:
            dwdp[26] = 1.0*TF_VIIa_X;
            break;
        case 2:
            dwdp[27] = -1.0*TF_VIIa_Xa;
            break;
        case 3:
            dwdp[27] = 1.0*TF_VIIa*Xa;
            break;
        case 4:
            dwdp[28] = -1.0*TF_VIIa_IX;
            break;
        case 5:
            dwdp[28] = 1.0*IX*TF_VIIa;
            break;
        case 6:
            dwdp[29] = 1.0*TF_VIIa_IX;
            break;
        case 7:
            dwdp[30] = 1.0*II*Xa;
            break;
        case 8:
            dwdp[1] = 1.0*IIa*VIII;
            break;
        case 9:
            dwdp[2] = -1.0*IXa_VIIIa;
            break;
        case 10:
            dwdp[2] = 1.0*IXa*VIIIa;
            break;
        case 11:
            dwdp[0] = 1.0*TF*VII;
            break;
        case 12:
            dwdp[3] = -1.0*IXa_VIIIa_X;
            break;
        case 13:
            dwdp[3] = 1.0*IXa_VIIIa*X;
            break;
        case 14:
            dwdp[4] = 1.0*IXa_VIIIa_X;
            break;
        case 15:
            dwdp[5] = -1.0*VIIIa1_L*VIIIa2;
            break;
        case 16:
            dwdp[5] = 1.0*VIIIa;
            break;
        case 17:
            dwdp[6] = 1.0*IXa_VIIIa_X;
            dwdp[7] = 1.0*IXa_VIIIa;
            break;
        case 18:
            dwdp[8] = 1.0*IIa*V;
            break;
        case 19:
            dwdp[9] = -1.0*Xa_Va;
            break;
        case 20:
            dwdp[9] = 1.0*Va*Xa;
            break;
        case 21:
            dwdp[10] = -1.0*Xa_Va_II;
            break;
        case 22:
            dwdp[13] = -1.0*TF_VIIa;
            break;
        case 23:
            dwdp[10] = 1.0*II*Xa_Va;
            break;
        case 24:
            dwdp[11] = 1.0*Xa_Va_II;
            break;
        case 25:
            dwdp[12] = 1.0*Xa_Va*mIIa;
            break;
        case 26:
            dwdp[14] = -1.0*Xa_TFPI;
            break;
        case 27:
            dwdp[14] = 1.0*TFPI*Xa;
            break;
        case 28:
            dwdp[15] = -1.0*TF_VIIa_Xa_TFPI;
            break;
        case 29:
            dwdp[15] = 1.0*TFPI*TF_VIIa_Xa;
            break;
        case 30:
            dwdp[16] = 1.0*TF_VIIa*Xa_TFPI;
            break;
        case 31:
            dwdp[17] = 1.0*ATIII*Xa;
            break;
        case 32:
            dwdp[18] = 1.0*ATIII*mIIa;
            break;
        case 33:
            dwdp[13] = 1.0*TF*VIIa;
            break;
        case 34:
            dwdp[19] = 1.0*ATIII*IXa;
            break;
        case 35:
            dwdp[20] = 1.0*ATIII*IIa;
            break;
        case 36:
            dwdp[21] = 1.0*ATIII*TF_VIIa;
            break;
        case 37:
            dwdp[22] = 1.0*TF_VIIa*VII;
            break;
        case 38:
            dwdp[23] = 1.0*VII*Xa;
            break;
        case 39:
            dwdp[24] = 1.0*IIa*VII;
            break;
        case 40:
            dwdp[25] = -1.0*TF_VIIa_X;
            break;
        case 41:
            dwdp[25] = 1.0*TF_VIIa*X;
            break;
    }
}