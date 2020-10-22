#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_Sivakumar2011c(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[1] = s1*s5;
            break;
        case 1:
            dwdp[1] = -s16;
            break;
        case 2:
            dwdp[2] = s16*s28;
            break;
        case 3:
            dwdp[2] = -s27;
            break;
        case 4:
            dwdp[3] = s121*s36;
            break;
        case 5:
            dwdp[3] = -s123;
            break;
        case 6:
            dwdp[4] = s123*s46;
            break;
        case 7:
            dwdp[4] = -s129;
            break;
        case 8:
            dwdp[5] = s123*s75;
            break;
        case 9:
            dwdp[5] = -s159;
            break;
        case 10:
            dwdp[6] = s36;
            break;
        case 11:
            dwdp[6] = -s232;
            break;
        case 12:
            dwdp[7] = s174*s232;
            break;
        case 13:
            dwdp[7] = -s176;
            break;
        case 14:
            dwdp[24] = s170*s176;
            break;
        case 15:
            dwdp[24] = -s179;
            break;
        case 16:
            dwdp[23] = s171*s179;
            break;
        case 17:
            dwdp[23] = -s183;
            break;
        case 18:
            dwdp[8] = s173*s183;
            break;
        case 19:
            dwdp[8] = -s188;
            break;
        case 20:
            dwdp[0] = s305;
            break;
        case 21:
            dwdp[9] = s252*s61;
            break;
        case 22:
            dwdp[9] = -s259;
            break;
        case 23:
            dwdp[10] = s259*s268;
            break;
        case 24:
            dwdp[10] = -s266;
            break;
        case 25:
            dwdp[11] = s266;
            break;
        case 26:
            dwdp[11] = -s155*s267;
            break;
        case 27:
            dwdp[12] = s267;
            break;
        case 28:
            dwdp[12] = -s260*s61;
            break;
        case 29:
            dwdp[13] = s159*s268;
            break;
        case 30:
            dwdp[13] = -s275;
            break;
        case 31:
            dwdp[14] = s275;
            break;
        case 32:
            dwdp[14] = -s101*s278;
            break;
        case 33:
            dwdp[15] = s278;
            break;
        case 34:
            dwdp[15] = -s164*s270;
            break;
        case 35:
            dwdp[16] = s286*s31;
            break;
        case 36:
            dwdp[16] = -s288;
            break;
        case 37:
            dwdp[17] = s102*s288;
            break;
        case 38:
            dwdp[17] = -s292;
            break;
        case 39:
            dwdp[18] = s292;
            break;
        case 40:
            dwdp[18] = -s37;
            break;
        case 41:
            dwdp[19] = s286;
            break;
        case 42:
            dwdp[19] = -s30;
            break;
        case 43:
            dwdp[20] = s239;
            break;
        case 44:
            dwdp[20] = -s5;
            break;
        case 45:
            dwdp[21] = s107*s30*s32;
            break;
        case 46:
            dwdp[21] = -s286*s30*s33;
            break;
        case 47:
            dwdp[22] = s129*s30*s32;
            break;
        case 48:
            dwdp[22] = -s245*s30*s33;
            break;
        case 49:
            dwdp[25] = s260;
            break;
        case 50:
            dwdp[26] = s270;
            break;
        case 51:
            dwdp[28] = kI_r86_s304*s245*pow(s32, 3)*s37/(kI_r86_s304 + s304);
            break;
        case 52:
            dwdp[28] = -kI_r86_s304*s252*pow(s33, 3)*s37/(kI_r86_s304 + s304);
            break;
        case 53:
            dwdp[27] = s172*s188;
            break;
        case 54:
            dwdp[27] = -s305;
            break;
        case 55:
            dwdp[28] = -kI_r86_s304*s37*(kass_r86_s37*s245*pow(s32, 3) - kdiss_r86_s37*s252*pow(s33, 3))/pow(kI_r86_s304 + s304, 2) + s37*(kass_r86_s37*s245*pow(s32, 3) - kdiss_r86_s37*s252*pow(s33, 3))/(kI_r86_s304 + s304);
            break;
    }
}