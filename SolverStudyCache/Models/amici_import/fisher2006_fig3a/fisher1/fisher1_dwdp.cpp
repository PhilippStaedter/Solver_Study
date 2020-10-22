#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_fisher1(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[0] = 1.13e-13*NFAT_Pi_Nuc;
            dwdp[3] = 2.6900000000000001e-13*NFAT_Pi_Cyt;
            break;
        case 1:
            dwdp[14] = 1.13e-13*NFAT_Act_C_Nuc;
            break;
        case 2:
            dwdp[13] = -1.13e-13*Act_C_Nuc*NFAT_Pi_Nuc;
            dwdp[16] = -2.6900000000000001e-13*Act_C_Cyt*NFAT_Pi_Cyt;
            break;
        case 3:
            dwdp[13] = 1.13e-13*NFAT_Pi_Act_C_Nuc;
            dwdp[16] = 2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt;
            break;
        case 4:
            dwdp[12] = -1.13e-13*NFAT_Pi_Act_C_Nuc;
            dwdp[15] = -2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt;
            break;
        case 5:
            dwdp[12] = 1.13e-13*NFAT_Act_C_Nuc;
            dwdp[15] = 2.6900000000000001e-13*NFAT_Act_C_Cyt;
            break;
        case 6:
            dwdp[2] = 2.6900000000000001e-13*NFAT_Act_C_Cyt;
            dwdp[9] = -1.13e-13*NFAT_Act_C_Nuc;
            break;
        case 7:
            dwdp[2] = -2.6900000000000001e-13*Act_C_Cyt*NFAT_Cyt;
            dwdp[9] = 1.13e-13*Act_C_Nuc*NFAT_Nuc;
            break;
        case 8:
            dwdp[10] = -2.6900000000000001e-13*NFAT_Cyt;
            break;
        case 9:
            dwdp[10] = 1.13e-13*NFAT_Nuc;
            break;
        case 10:
            dwdp[4] = 2.6900000000000001e-13*pow(Ca_Cyt, 3)*Inact_C_Cyt;
            dwdp[5] = 1.13e-13*pow(Ca_Nuc, 3)*Inact_C_Nuc;
            break;
        case 11:
            dwdp[0] = -1.13e-13*NFAT_Nuc;
            dwdp[3] = -2.6900000000000001e-13*NFAT_Cyt;
            break;
        case 12:
            dwdp[4] = -2.6900000000000001e-13*Act_C_Cyt;
            dwdp[5] = -1.13e-13*Act_C_Nuc;
            break;
        case 13:
            dwdp[7] = 2.6900000000000001e-13*Ca_Cyt;
            break;
        case 14:
            dwdp[7] = -1.13e-13*Ca_Nuc;
            break;
        case 15:
            dwdp[1] = 2.6900000000000001e-13*NFAT_Pi_Cyt;
            break;
        case 16:
            dwdp[1] = -1.13e-13*NFAT_Pi_Nuc;
            break;
        case 17:
            dwdp[6] = 2.6900000000000001e-13*Inact_C_Cyt;
            dwdp[11] = -2.6900000000000001e-13*Act_C_Cyt;
            break;
        case 18:
            dwdp[6] = -1.13e-13*Inact_C_Nuc;
            dwdp[11] = 1.13e-13*Act_C_Nuc;
            break;
        case 19:
            dwdp[8] = 2.6900000000000001e-13*NFAT_Pi_Act_C_Cyt;
            break;
        case 20:
            dwdp[8] = -1.13e-13*NFAT_Pi_Act_C_Nuc;
            break;
        case 21:
            dwdp[14] = -2.6900000000000001e-13*NFAT_Act_C_Cyt;
            break;
    }
}