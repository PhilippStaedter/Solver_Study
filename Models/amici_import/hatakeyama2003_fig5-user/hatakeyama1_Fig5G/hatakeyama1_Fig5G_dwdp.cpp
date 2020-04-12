#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdp_hatakeyama1_Fig5G(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp, const int ip){
    switch(ip) {
        case 0:
            dwdp[9] = -1.0*ShP*V10/pow(K10 + ShP, 2);
            break;
        case 1:
            dwdp[10] = -1.0*RasGDP*ShGS*k11/pow(K11 + RasGDP, 2);
            break;
        case 2:
            dwdp[11] = -1.0*RasGTP*V12/pow(K12 + RasGTP, 2);
            break;
        case 3:
            dwdp[12] = -1.0*Raf*RasGTP*k13/pow(K13 + Raf, 2);
            break;
        case 4:
            dwdp[13] = -1.0*Rafstar*k14*(AktPIPP + E)/pow(K14 + Rafstar, 2);
            break;
        case 5:
            dwdp[14] = 1.0*MEK*Rafstar*k15*(-1 - MEKP/K17)/pow(K15*(1 + MEKP/K17) + MEK, 2);
            dwdp[16] = 1.0*K17*MEK*MEKP*Rafstar*k17/(pow(K15, 2)*pow(K17*(1 + MEK/K15) + MEKP, 2));
            break;
        case 6:
            dwdp[15] = 1.0*MEKP*PP2A*k16*(-AktPIP/K31 - AktPIPP/K33 - 1 - MEKPP/K18)/pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2);
            dwdp[17] = 1.0*K18*MEKP*MEKPP*PP2A*k18/(pow(K16, 2)*pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2));
            dwdp[30] = 1.0*AktPIP*K31*MEKP*PP2A*k31/(pow(K16, 2)*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
            dwdp[32] = 1.0*AktPIPP*K33*MEKP*PP2A*k33/(pow(K16, 2)*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
            break;
        case 7:
            dwdp[14] = 1.0*K15*MEK*MEKP*Rafstar*k15/(pow(K17, 2)*pow(K15*(1 + MEKP/K17) + MEK, 2));
            dwdp[16] = 1.0*MEKP*Rafstar*k17*(-1 - MEK/K15)/pow(K17*(1 + MEK/K15) + MEKP, 2);
            break;
        case 8:
            dwdp[15] = 1.0*K16*MEKP*MEKPP*PP2A*k16/(pow(K18, 2)*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
            dwdp[17] = 1.0*MEKPP*PP2A*k18*(-AktPIPP/K33 - AktPIPP/K31 - 1 - MEKP/K16)/pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2);
            dwdp[30] = 1.0*AktPIP*K31*MEKPP*PP2A*k31/(pow(K18, 2)*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
            dwdp[32] = 1.0*AktPIPP*K33*MEKPP*PP2A*k33/(pow(K18, 2)*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
            break;
        case 9:
            dwdp[18] = 1.0*ERK*MEKPP*k19*(-ERKP/K21 - 1)/pow(ERK + K19*(ERKP/K21 + 1), 2);
            dwdp[20] = 1.0*ERK*ERKP*K21*MEKPP*k21/(pow(K19, 2)*pow(ERKP + K21*(ERK/K19 + 1), 2));
            break;
        case 10:
            dwdp[19] = 1.0*ERKP*MKP3*k20*(-ERKPP/K22 - 1)/pow(ERKP + K20*(ERKPP/K22 + 1), 2);
            dwdp[21] = 1.0*ERKP*ERKPP*K22*MKP3*k22/(pow(K20, 2)*pow(ERKPP + K22*(ERKP/K20 + 1), 2));
            break;
        case 11:
            dwdp[18] = 1.0*ERK*ERKP*K19*MEKPP*k19/(pow(K21, 2)*pow(ERK + K19*(ERKP/K21 + 1), 2));
            dwdp[20] = 1.0*ERKP*MEKPP*k21*(-ERK/K19 - 1)/pow(ERKP + K21*(ERK/K19 + 1), 2);
            break;
        case 12:
            dwdp[19] = 1.0*ERKP*ERKPP*K20*MKP3*k20/(pow(K22, 2)*pow(ERKP + K20*(ERKPP/K22 + 1), 2));
            dwdp[21] = 1.0*ERKPP*MKP3*k22*(-ERKP/K20 - 1)/pow(ERKPP + K22*(ERKP/K20 + 1), 2);
            break;
        case 13:
            dwdp[25] = -1.0*PI3Kstar*V26/pow(K26 + PI3Kstar, 2);
            break;
        case 14:
            dwdp[26] = -1.0*PI3Kstar*P_I*k27/pow(K27 + P_I, 2);
            break;
        case 15:
            dwdp[27] = -1.0*PIP3*V28/pow(K28 + PIP3, 2);
            break;
        case 16:
            dwdp[29] = 1.0*AktPIP3*V30*(-AktPIP/K32 - 1)/pow(AktPIP3 + K30*(AktPIP/K32 + 1), 2);
            dwdp[31] = 1.0*AktPIP*AktPIP3*K32*V32/(pow(K30, 2)*pow(AktPIP + K32*(AktPIP3/K30 + 1), 2));
            break;
        case 17:
            dwdp[15] = 1.0*AktPIP*K16*MEKP*PP2A*k16/(pow(K31, 2)*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
            dwdp[17] = 1.0*AktPIPP*K18*MEKPP*PP2A*k18/(pow(K31, 2)*pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2));
            dwdp[30] = 1.0*AktPIP*PP2A*k31*(-AktPIPP/K33 - 1 - MEKPP/K18 - MEKP/K16)/pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2);
            dwdp[32] = 1.0*AktPIP*AktPIPP*K33*PP2A*k33/(pow(K31, 2)*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
            break;
        case 18:
            dwdp[29] = 1.0*AktPIP*AktPIP3*K30*V30/(pow(K32, 2)*pow(AktPIP3 + K30*(AktPIP/K32 + 1), 2));
            dwdp[31] = 1.0*AktPIP*V32*(-AktPIP3/K30 - 1)/pow(AktPIP + K32*(AktPIP3/K30 + 1), 2);
            break;
        case 19:
            dwdp[15] = 1.0*AktPIPP*K16*MEKP*PP2A*k16/(pow(K33, 2)*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
            dwdp[17] = 1.0*AktPIPP*K18*MEKPP*PP2A*k18/(pow(K33, 2)*pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2));
            dwdp[30] = 1.0*AktPIP*AktPIPP*K31*PP2A*k31/(pow(K33, 2)*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
            dwdp[32] = 1.0*AktPIPP*PP2A*k33*(-AktPIP/K31 - 1 - MEKPP/K18 - MEKP/K16)/pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2);
            break;
        case 20:
            dwdp[3] = -1.0*RP*V4/pow(K4 + RP, 2);
            break;
        case 21:
            dwdp[9] = 1.0*ShP/(K10 + ShP);
            break;
        case 22:
            dwdp[11] = 1.0*RasGTP/(K12 + RasGTP);
            break;
        case 23:
            dwdp[25] = 1.0*PI3Kstar/(K26 + PI3Kstar);
            break;
        case 24:
            dwdp[27] = 1.0*PIP3/(K28 + PIP3);
            break;
        case 25:
            dwdp[29] = 1.0*AktPIP3/(AktPIP3 + K30*(AktPIP/K32 + 1));
            break;
        case 26:
            dwdp[31] = 1.0*AktPIP/(AktPIP + K32*(AktPIP3/K30 + 1));
            break;
        case 27:
            dwdp[3] = 1.0*RP/(K4 + RP);
            break;
        case 28:
            dwdp[0] = 1.0*HRG*R;
            break;
        case 29:
            dwdp[10] = 1.0*RasGDP*ShGS/(K11 + RasGDP);
            break;
        case 30:
            dwdp[12] = 1.0*Raf*RasGTP/(K13 + Raf);
            break;
        case 31:
            dwdp[13] = 1.0*Rafstar*(AktPIPP + E)/(K14 + Rafstar);
            break;
        case 32:
            dwdp[14] = 1.0*MEK*Rafstar/(K15*(1 + MEKP/K17) + MEK);
            break;
        case 33:
            dwdp[15] = 1.0*MEKP*PP2A/(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP);
            break;
        case 34:
            dwdp[16] = 1.0*MEKP*Rafstar/(K17*(1 + MEK/K15) + MEKP);
            break;
        case 35:
            dwdp[17] = 1.0*MEKPP*PP2A/(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP);
            break;
        case 36:
            dwdp[18] = 1.0*ERK*MEKPP/(ERK + K19*(ERKP/K21 + 1));
            break;
        case 37:
            dwdp[1] = 1.0*pow(RHRG, 2);
            break;
        case 38:
            dwdp[19] = 1.0*ERKP*MKP3/(ERKP + K20*(ERKPP/K22 + 1));
            break;
        case 39:
            dwdp[20] = 1.0*ERKP*MEKPP/(ERKP + K21*(ERK/K19 + 1));
            break;
        case 40:
            dwdp[21] = 1.0*ERKPP*MKP3/(ERKPP + K22*(ERKP/K20 + 1));
            break;
        case 41:
            dwdp[22] = 1.0*PI3K*RP;
            break;
        case 42:
            dwdp[23] = 1.0*RPI3K;
            break;
        case 43:
            dwdp[24] = 1.0*RPI3Kstar;
            break;
        case 44:
            dwdp[26] = 1.0*PI3Kstar*P_I/(K27 + P_I);
            break;
        case 45:
            dwdp[28] = 1.0*Akt*PIP3;
            break;
        case 46:
            dwdp[2] = 1.0*RHRG2;
            break;
        case 47:
            dwdp[30] = 1.0*AktPIP*PP2A/(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16));
            break;
        case 48:
            dwdp[32] = 1.0*AktPIPP*PP2A/(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16));
            break;
        case 49:
            dwdp[33] = 1.0*RP;
            break;
        case 50:
            dwdp[4] = 1.0*RP*Shc;
            break;
        case 51:
            dwdp[5] = 1.0*RShc;
            break;
        case 52:
            dwdp[6] = 1.0*GS*RShP;
            break;
        case 53:
            dwdp[7] = 1.0*RShGS;
            break;
        case 54:
            dwdp[8] = 1.0*ShGS;
            break;
        case 55:
            dwdp[0] = -1.0*RHRG;
            break;
        case 56:
            dwdp[1] = -1.0*RHRG2;
            break;
        case 57:
            dwdp[22] = -1.0*RPI3K;
            break;
        case 58:
            dwdp[23] = -1.0*RPI3Kstar;
            break;
        case 59:
            dwdp[24] = -1.0*PI3Kstar*RP;
            break;
        case 60:
            dwdp[28] = -1.0*AktPIP3;
            break;
        case 61:
            dwdp[2] = -1.0*RP;
            break;
        case 62:
            dwdp[33] = -1.0*internalization;
            break;
        case 63:
            dwdp[4] = -1.0*RShc;
            break;
        case 64:
            dwdp[5] = -1.0*RShP;
            break;
        case 65:
            dwdp[6] = -1.0*RShGS;
            break;
        case 66:
            dwdp[7] = -1.0*RP*ShGS;
            break;
        case 67:
            dwdp[8] = -1.0*GS*ShP;
            break;
    }
}