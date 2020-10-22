#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_hatakeyama1_Fig5G(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*HRG*R*k1 - 1.0*RHRG*k_1;
    w[1] = 1.0*pow(RHRG, 2)*k2 - 1.0*RHRG2*k_2;
    w[2] = 1.0*RHRG2*k3 - 1.0*RP*k_3;
    w[3] = 1.0*RP*V4/(K4 + RP);
    w[4] = 1.0*RP*Shc*k5 - 1.0*RShc*k_5;
    w[5] = -1.0*RShP*k_6 + 1.0*RShc*k6;
    w[6] = 1.0*GS*RShP*k7 - 1.0*RShGS*k_7;
    w[7] = -1.0*RP*ShGS*k_8 + 1.0*RShGS*k8;
    w[8] = -1.0*GS*ShP*k_9 + 1.0*ShGS*k9;
    w[9] = 1.0*ShP*V10/(K10 + ShP);
    w[10] = 1.0*RasGDP*ShGS*k11/(K11 + RasGDP);
    w[11] = 1.0*RasGTP*V12/(K12 + RasGTP);
    w[12] = 1.0*Raf*RasGTP*k13/(K13 + Raf);
    w[13] = 1.0*Rafstar*k14*(AktPIPP + E)/(K14 + Rafstar);
    w[14] = 1.0*MEK*Rafstar*k15/(K15*(1 + MEKP/K17) + MEK);
    w[15] = 1.0*MEKP*PP2A*k16/(K16*(AktPIP/K31 + AktPIPP/K33 + 1 + MEKPP/K18) + MEKP);
    w[16] = 1.0*MEKP*Rafstar*k17/(K17*(1 + MEK/K15) + MEKP);
    w[17] = 1.0*MEKPP*PP2A*k18/(K18*(AktPIPP/K33 + AktPIPP/K31 + 1 + MEKP/K16) + MEKPP);
    w[18] = 1.0*ERK*MEKPP*k19/(ERK + K19*(ERKP/K21 + 1));
    w[19] = 1.0*ERKP*MKP3*k20/(ERKP + K20*(ERKPP/K22 + 1));
    w[20] = 1.0*ERKP*MEKPP*k21/(ERKP + K21*(ERK/K19 + 1));
    w[21] = 1.0*ERKPP*MKP3*k22/(ERKPP + K22*(ERKP/K20 + 1));
    w[22] = 1.0*PI3K*RP*k23 - 1.0*RPI3K*k_23;
    w[23] = 1.0*RPI3K*k24 - 1.0*RPI3Kstar*k_24;
    w[24] = -1.0*PI3Kstar*RP*k_25 + 1.0*RPI3Kstar*k25;
    w[25] = 1.0*PI3Kstar*V26/(K26 + PI3Kstar);
    w[26] = 1.0*PI3Kstar*P_I*k27/(K27 + P_I);
    w[27] = 1.0*PIP3*V28/(K28 + PIP3);
    w[28] = 1.0*Akt*PIP3*k29 - 1.0*AktPIP3*k_29;
    w[29] = 1.0*AktPIP3*V30/(AktPIP3 + K30*(AktPIP/K32 + 1));
    w[30] = 1.0*AktPIP*PP2A*k31/(AktPIP + K31*(AktPIPP/K33 + 1 + MEKPP/K18 + MEKP/K16));
    w[31] = 1.0*AktPIP*V32/(AktPIP + K32*(AktPIP3/K30 + 1));
    w[32] = 1.0*AktPIPP*PP2A*k33/(AktPIPP + K33*(AktPIP/K31 + 1 + MEKPP/K18 + MEKP/K16));
    w[33] = 1.0*RP*k34 - 1.0*internalization*k_34;
}