#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_hatakeyama1_Fig5G(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*PIP3*k29;
    dwdx[1] = -1.0*K16*MEKP*PP2A*k16/(K31*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
    dwdx[2] = -1.0*AktPIP3*K30*V30/(K32*pow(AktPIP3 + K30*(AktPIP/K32 + 1), 2));
    dwdx[3] = -1.0*AktPIP*PP2A*k31/pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2) + 1.0*PP2A*k31/(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16));
    dwdx[4] = -1.0*AktPIP*V32/pow(AktPIP + K32*(AktPIP3/K30 + 1), 2) + 1.0*V32/(AktPIP + K32*(AktPIP3/K30 + 1));
    dwdx[5] = -1.0*AktPIPP*K33*PP2A*k33/(K31*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[6] = -1.0*k_29;
    dwdx[7] = -1.0*AktPIP3*V30/pow(AktPIP3 + K30*(AktPIP/K32 + 1), 2) + 1.0*V30/(AktPIP3 + K30*(AktPIP/K32 + 1));
    dwdx[8] = -1.0*AktPIP*K32*V32/(K30*pow(AktPIP + K32*(AktPIP3/K30 + 1), 2));
    dwdx[9] = 1.0*Rafstar*k14/(K14 + Rafstar);
    dwdx[10] = -1.0*K16*MEKP*PP2A*k16/(K33*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
    dwdx[11] = -1.0*K18*MEKPP*PP2A*k18*(1.0/K33 + 1.0/K31)/pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2);
    dwdx[12] = -1.0*AktPIP*K31*PP2A*k31/(K33*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[13] = -1.0*AktPIPP*PP2A*k33/pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2) + 1.0*PP2A*k33/(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16));
    dwdx[14] = 1.0*Rafstar*k14/(K14 + Rafstar);
    dwdx[15] = -1.0*ERK*MEKPP*k19/pow(ERK + K19*(ERKP/K21 + 1), 2) + 1.0*MEKPP*k19/(ERK + K19*(ERKP/K21 + 1));
    dwdx[16] = -1.0*ERKP*K21*MEKPP*k21/(K19*pow(ERKP + K21*(ERK/K19 + 1), 2));
    dwdx[17] = -1.0*ERK*K19*MEKPP*k19/(K21*pow(ERK + K19*(ERKP/K21 + 1), 2));
    dwdx[18] = -1.0*ERKP*MKP3*k20/pow(ERKP + K20*(ERKPP/K22 + 1), 2) + 1.0*MKP3*k20/(ERKP + K20*(ERKPP/K22 + 1));
    dwdx[19] = -1.0*ERKP*MEKPP*k21/pow(ERKP + K21*(ERK/K19 + 1), 2) + 1.0*MEKPP*k21/(ERKP + K21*(ERK/K19 + 1));
    dwdx[20] = -1.0*ERKPP*K22*MKP3*k22/(K20*pow(ERKPP + K22*(ERKP/K20 + 1), 2));
    dwdx[21] = -1.0*ERKP*K20*MKP3*k20/(K22*pow(ERKP + K20*(ERKPP/K22 + 1), 2));
    dwdx[22] = -1.0*ERKPP*MKP3*k22/pow(ERKPP + K22*(ERKP/K20 + 1), 2) + 1.0*MKP3*k22/(ERKPP + K22*(ERKP/K20 + 1));
    dwdx[23] = 1.0*RShP*k7;
    dwdx[24] = -1.0*ShP*k_9;
    dwdx[25] = 1.0*R*k1;
    dwdx[26] = -1.0*MEK*Rafstar*k15/pow(K15*(1 + MEKP/K17) + MEK, 2) + 1.0*Rafstar*k15/(K15*(1 + MEKP/K17) + MEK);
    dwdx[27] = -1.0*K17*MEKP*Rafstar*k17/(K15*pow(K17*(1 + MEK/K15) + MEKP, 2));
    dwdx[28] = -1.0*K15*MEK*Rafstar*k15/(K17*pow(K15*(1 + MEKP/K17) + MEK, 2));
    dwdx[29] = -1.0*MEKP*PP2A*k16/pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2) + 1.0*PP2A*k16/(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP);
    dwdx[30] = -1.0*MEKP*Rafstar*k17/pow(K17*(1 + MEK/K15) + MEKP, 2) + 1.0*Rafstar*k17/(K17*(1 + MEK/K15) + MEKP);
    dwdx[31] = -1.0*K18*MEKPP*PP2A*k18/(K16*pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2));
    dwdx[32] = -1.0*AktPIP*K31*PP2A*k31/(K16*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[33] = -1.0*AktPIPP*K33*PP2A*k33/(K16*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[34] = -1.0*K16*MEKP*PP2A*k16/(K18*pow(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP, 2));
    dwdx[35] = -1.0*MEKPP*PP2A*k18/pow(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP, 2) + 1.0*PP2A*k18/(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP);
    dwdx[36] = 1.0*ERK*k19/(ERK + K19*(ERKP/K21 + 1));
    dwdx[37] = 1.0*ERKP*k21/(ERKP + K21*(ERK/K19 + 1));
    dwdx[38] = -1.0*AktPIP*K31*PP2A*k31/(K18*pow(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[39] = -1.0*AktPIPP*K33*PP2A*k33/(K18*pow(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16), 2));
    dwdx[40] = 1.0*ERKP*k20/(ERKP + K20*(ERKPP/K22 + 1));
    dwdx[41] = 1.0*ERKPP*k22/(ERKPP + K22*(ERKP/K20 + 1));
    dwdx[42] = 1.0*RP*k23;
    dwdx[43] = -1.0*RP*k_25;
    dwdx[44] = -1.0*PI3Kstar*V26/pow(K26 + PI3Kstar, 2) + 1.0*V26/(K26 + PI3Kstar);
    dwdx[45] = 1.0*P_I*k27/(K27 + P_I);
    dwdx[46] = -1.0*PIP3*V28/pow(K28 + PIP3, 2) + 1.0*V28/(K28 + PIP3);
    dwdx[47] = 1.0*Akt*k29;
    dwdx[48] = 1.0*MEKP*k16/(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP);
    dwdx[49] = 1.0*MEKPP*k18/(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP);
    dwdx[50] = 1.0*AktPIP*k31/(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16));
    dwdx[51] = 1.0*AktPIPP*k33/(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16));
    dwdx[52] = -1.0*PI3Kstar*P_I*k27/pow(K27 + P_I, 2) + 1.0*PI3Kstar*k27/(K27 + P_I);
    dwdx[53] = 1.0*HRG*k1;
    dwdx[54] = -1.0*k_1;
    dwdx[55] = 2.0*RHRG*k2;
    dwdx[56] = -1.0*k_2;
    dwdx[57] = 1.0*k3;
    dwdx[58] = -1.0*k_3;
    dwdx[59] = -1.0*RP*V4/pow(K4 + RP, 2) + 1.0*V4/(K4 + RP);
    dwdx[60] = 1.0*Shc*k5;
    dwdx[61] = -1.0*ShGS*k_8;
    dwdx[62] = 1.0*PI3K*k23;
    dwdx[63] = -1.0*PI3Kstar*k_25;
    dwdx[64] = 1.0*k34;
    dwdx[65] = -1.0*k_23;
    dwdx[66] = 1.0*k24;
    dwdx[67] = -1.0*k_24;
    dwdx[68] = 1.0*k25;
    dwdx[69] = -1.0*k_7;
    dwdx[70] = 1.0*k8;
    dwdx[71] = -1.0*k_6;
    dwdx[72] = 1.0*GS*k7;
    dwdx[73] = -1.0*k_5;
    dwdx[74] = 1.0*k6;
    dwdx[75] = -1.0*Raf*RasGTP*k13/pow(K13 + Raf, 2) + 1.0*RasGTP*k13/(K13 + Raf);
    dwdx[76] = -1.0*Rafstar*k14*(AktPIPP + E)/pow(K14 + Rafstar, 2) + 1.0*k14*(AktPIPP + E)/(K14 + Rafstar);
    dwdx[77] = 1.0*MEK*k15/(K15*(1 + MEKP/K17) + MEK);
    dwdx[78] = 1.0*MEKP*k17/(K17*(1 + MEK/K15) + MEKP);
    dwdx[79] = -1.0*RasGDP*ShGS*k11/pow(K11 + RasGDP, 2) + 1.0*ShGS*k11/(K11 + RasGDP);
    dwdx[80] = -1.0*RasGTP*V12/pow(K12 + RasGTP, 2) + 1.0*V12/(K12 + RasGTP);
    dwdx[81] = 1.0*Raf*k13/(K13 + Raf);
    dwdx[82] = -1.0*RP*k_8;
    dwdx[83] = 1.0*k9;
    dwdx[84] = 1.0*RasGDP*k11/(K11 + RasGDP);
    dwdx[85] = -1.0*GS*k_9;
    dwdx[86] = -1.0*ShP*V10/pow(K10 + ShP, 2) + 1.0*V10/(K10 + ShP);
    dwdx[87] = 1.0*RP*k5;
    dwdx[88] = -1.0*k_34;
}