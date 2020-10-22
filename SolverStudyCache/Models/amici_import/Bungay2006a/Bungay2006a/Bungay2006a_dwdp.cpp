#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Bungay2006a(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.0*II_f*LIPID/nva;
            break;
        case 1:
            dwdp[0] = -1.0*II_f*LIPID*konII/pow(nva, 2);
            dwdp[1] = -1.0*LIPID*konmIIa*mIIa_f/pow(nva, 2);
            dwdp[2] = -1.0*LIPID*V_f*konV/pow(nva, 2);
            dwdp[3] = -1.0*LIPID*Va_f*konVa/pow(nva, 2);
            dwdp[4] = -1.0*LIPID*VII_f*konVII/pow(nva, 2);
            dwdp[5] = -1.0*LIPID*VIIa_f*konVIIa/pow(nva, 2);
            dwdp[6] = -1.0*LIPID*X_f*konX/pow(nva, 2);
            dwdp[7] = -LIPID*Xa_f*konXa/pow(nva, 2);
            dwdp[8] = -1.0*APC_f*LIPID*konAPC/pow(nva, 2);
            dwdp[9] = -1.0*LIPID*PS_f*konPS/pow(nva, 2);
            dwdp[10] = -1.0*LIPID*Vai_f*konVai/pow(nva, 2);
            dwdp[11] = -1.0*LIPID*PC_f*konPC/pow(nva, 2);
            break;
        case 2:
            dwdp[0] = -1.0*II_l;
            break;
        case 3:
            dwdp[1] = 1.0*LIPID*mIIa_f/nva;
            break;
        case 4:
            dwdp[1] = -1.0*mIIa_l;
            break;
        case 5:
            dwdp[2] = 1.0*LIPID*V_f/nva;
            break;
        case 6:
            dwdp[2] = -1.0*V_l;
            break;
        case 7:
            dwdp[3] = 1.0*LIPID*Va_f/nva;
            break;
        case 8:
            dwdp[3] = -1.0*Va_l;
            break;
        case 9:
            dwdp[4] = 1.0*LIPID*VII_f/nva;
            break;
        case 10:
            dwdp[4] = -1.0*VII_l;
            break;
        case 11:
            dwdp[5] = 1.0*LIPID*VIIa_f/nva;
            break;
        case 12:
            dwdp[5] = -1.0*VIIa_l;
            break;
        case 13:
            dwdp[6] = 1.0*LIPID*X_f/nva;
            break;
        case 14:
            dwdp[6] = -1.0*X_l;
            break;
        case 15:
            dwdp[7] = LIPID*Xa_f/nva;
            break;
        case 16:
            dwdp[7] = -Xa_l;
            break;
        case 17:
            dwdp[8] = 1.0*APC_f*LIPID/nva;
            break;
        case 18:
            dwdp[8] = -1.0*APC_l;
            break;
        case 19:
            dwdp[9] = 1.0*LIPID*PS_f/nva;
            break;
        case 20:
            dwdp[9] = -1.0*PS_l;
            break;
        case 21:
            dwdp[10] = 1.0*LIPID*Vai_f/nva;
            break;
        case 22:
            dwdp[10] = -1.0*Vai_l;
            break;
        case 23:
            dwdp[11] = 1.0*LIPID*PC_f/nva;
            break;
        case 24:
            dwdp[11] = -1.0*PC_l;
            break;
        case 25:
            dwdp[12] = 1.0*TF_l*VIIa_l;
            break;
        case 26:
            dwdp[12] = -1.0*TF_VIIa_l;
            break;
        case 27:
            dwdp[13] = 1.0*TF_l*VII_l;
            break;
        case 28:
            dwdp[13] = -1.0*TF_VII_l;
            break;
        case 29:
            dwdp[14] = 1.0*TF_VIIa_l*X_l;
            break;
        case 30:
            dwdp[14] = -1.0*TF_VIIa_X_l;
            break;
        case 31:
            dwdp[15] = 1.0*TF_VIIa_X_l;
            break;
        case 32:
            dwdp[17] = 1.0*TF_VII_l*Xa_l;
            break;
        case 33:
            dwdp[17] = -1.0*TF_VII_Xa_l;
            break;
        case 34:
            dwdp[18] = 1.0*TF_VII_Xa_l;
            break;
        case 35:
            dwdp[19] = 1.0*Va_l*Xa_l;
            break;
        case 36:
            dwdp[19] = -1.0*Xa_Va_l;
            break;
        case 37:
            dwdp[20] = 1.0*V_l*Xa_l;
            break;
        case 38:
            dwdp[20] = -1.0*V_Xa_l;
            break;
        case 39:
            dwdp[21] = 1.0*V_Xa_l;
            break;
        case 40:
            dwdp[22] = 1.0*IIa_f*V_l;
            break;
        case 41:
            dwdp[22] = -1.0*V_IIa_l;
            break;
        case 42:
            dwdp[23] = 1.0*V_IIa_l;
            break;
        case 43:
            dwdp[24] = 1.0*II_l*Xa_Va_l;
            break;
        case 44:
            dwdp[24] = -1.0*Xa_Va_II_l;
            break;
        case 45:
            dwdp[25] = 1.0*Xa_Va_l*mIIa_l;
            break;
        case 46:
            dwdp[25] = -1.0*Xa_Va_mIIa_l;
            break;
        case 47:
            dwdp[26] = 1.0*Xa_Va_II_l;
            break;
        case 48:
            dwdp[27] = 1.0*Xa_Va_mIIa_l;
            break;
        case 49:
            dwdp[28] = 1.0*VII_l*Xa_l;
            break;
        case 50:
            dwdp[28] = -1.0*VII_Xa_l;
            break;
        case 51:
            dwdp[29] = 1.0*VII_Xa_l;
            break;
        case 52:
            dwdp[30] = 1.0*APC_PS_l*Va_l;
            break;
        case 53:
            dwdp[30] = -1.0*APC_PS_Va_l;
            break;
        case 54:
            dwdp[31] = 1.0*APC_PS_Va_l;
            break;
        case 55:
            dwdp[32] = 1.0*TFPI_f*Xa_f;
            break;
        case 56:
            dwdp[32] = -1.0*TFPI_Xa_l;
            break;
        case 57:
            dwdp[33] = 1.0*TFPI_Xa_l*TF_VIIa_l;
            break;
        case 58:
            dwdp[33] = -1.0*TFPI_Xa_TF_VIIa_l;
            break;
        case 59:
            dwdp[34] = 1.0*AT_f*Xa_f;
            break;
        case 60:
            dwdp[35] = 1.0*AT_f*IIa_f;
            break;
        case 61:
            dwdp[36] = 1.0*V_l*mIIa_l;
            break;
        case 62:
            dwdp[36] = -1.0*V_mIIa_l;
            break;
        case 63:
            dwdp[37] = 1.0*V_mIIa_l;
            break;
        case 64:
            dwdp[38] = 1.0*IIa_f*TM_l;
            break;
        case 65:
            dwdp[38] = -1.0*IIa_TM_l;
            break;
        case 66:
            dwdp[39] = 1.0*IIa_TM_l*PC_l;
            break;
        case 67:
            dwdp[39] = -1.0*IIa_TM_PC_l;
            break;
        case 68:
            dwdp[40] = 1.0*IIa_TM_PC_l;
            break;
        case 69:
            dwdp[41] = 1.0*AT_f*mIIa_f;
            break;
        case 70:
            dwdp[42] = 1.0*APC_l*PS_l;
            break;
        case 71:
            dwdp[42] = -1.0*APC_PS_l;
            break;
        case 72:
            dwdp[16] = 1.0*TF_VIIa_Xa_l;
            break;
        case 73:
            dwdp[43] = 1.0*IIa_f*alpha2M_l;
            break;
        case 74:
            dwdp[44] = 1.0*Xa_f*alpha2M_l;
            break;
    }
}