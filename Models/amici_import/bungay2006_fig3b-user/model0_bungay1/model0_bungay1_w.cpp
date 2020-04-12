#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bungay1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*II_f*LIPID*konII/nva - 1.0*II_l*koffII;
    w[1] = 1.0*IXa_f*LIPID*konIXa/nva - 1.0*IXa_l*koffIXa;
    w[2] = 1.0*LIPID*X_f*konX/nva - 1.0*X_l*koffX;
    w[3] = LIPID*Xa_f*konXa/nva - Xa_l*koffXa;
    w[4] = 1.0*APC_f*LIPID*konAPC/nva - 1.0*APC_l*koffAPC;
    w[5] = 1.0*LIPID*PS_f*konPS/nva - 1.0*PS_l*koffPS;
    w[6] = 1.0*LIPID*VIIIai_f*konVIIIai/nva - 1.0*VIIIai_l*koffVIIIai;
    w[7] = 1.0*LIPID*Vai_f*konVai/nva - 1.0*Vai_l*koffVai;
    w[8] = 1.0*LIPID*PC_f*konPC/nva - 1.0*PC_l*koffPC;
    w[9] = 1.0*LIPID*konmIIa*mIIa_f/nva - 1.0*koffmIIa*mIIa_l;
    w[10] = 1.0*LIPID*V_f*konV/nva - 1.0*V_l*koffV;
    w[11] = 1.0*LIPID*Va_f*konVa/nva - 1.0*Va_l*koffVa;
    w[12] = 1.0*LIPID*VII_f*konVII/nva - 1.0*VII_l*koffVII;
    w[13] = 1.0*LIPID*VIIa_f*konVIIa/nva - 1.0*VIIa_l*koffVIIa;
    w[14] = 1.0*LIPID*VIII_f*konVIII/nva - 1.0*VIII_l*koffVIII;
    w[15] = 1.0*LIPID*VIIIa_f*konVIIIa/nva - 1.0*VIIIa_l*koffVIIIa;
    w[16] = 1.0*IX_f*LIPID*konIX/nva - 1.0*IX_l*koffIX;
    w[17] = -1.0*TF_VIIa_l*k2 + 1.0*TF_l*VIIa_l*k1;
    w[18] = -VIII_Xa_l*k25 + VIII_l*Xa_l*k24;
    w[19] = 1.0*VIII_Xa_l*k26;
    w[20] = 1.0*IIa_f*V_l*k27 - 1.0*V_IIa_l*k28;
    w[21] = 1.0*V_IIa_l*k29;
    w[22] = 1.0*IIa_f*VIII_l*k30 - 1.0*VIII_IIa_l*k31;
    w[23] = 1.0*VIII_IIa_l*k32;
    w[24] = 1.0*II_l*Xa_Va_l*k33 - 1.0*Xa_Va_II_l*k34;
    w[25] = 1.0*Xa_Va_l*k35*mIIa_l - 1.0*Xa_Va_mIIa_l*k36;
    w[26] = 1.0*Xa_Va_II_l*k37;
    w[27] = 1.0*Xa_Va_mIIa_l*k38;
    w[28] = -1.0*VII_Xa_l*k40 + 1.0*VII_l*Xa_l*k39;
    w[29] = 1.0*VII_Xa_l*k41;
    w[30] = 1.0*IIa_f*XI_f*k42 - 1.0*XI_IIa_l*k43;
    w[31] = 1.0*XI_IIa_l*k44;
    w[32] = -1.0*APC_PS_VIIIa_l*k46 + 1.0*APC_PS_l*VIIIa_l*k45;
    w[33] = 1.0*APC_PS_VIIIa_l*k47;
    w[34] = -1.0*APC_PS_Va_l*k49 + 1.0*APC_PS_l*Va_l*k48;
    w[35] = 1.0*APC_PS_Va_l*k50;
    w[36] = -1.0*TF_VII_l*k4 + 1.0*TF_l*VII_l*k3;
    w[37] = -1.0*TFPI_Xa_l*k52 + 1.0*TFPI_f*Xa_f*k51;
    w[38] = -1.0*TFPI_Xa_TF_VIIa_l*k54 + 1.0*TFPI_Xa_l*TF_VIIa_l*k53;
    w[39] = 1.0*AT_f*IXa_f*k55;
    w[40] = 1.0*AT_f*Xa_f*k56;
    w[41] = 1.0*AT_f*IIa_f*k57;
    w[42] = 1.0*V_l*k58*mIIa_l - 1.0*V_mIIa_l*k59;
    w[43] = 1.0*V_mIIa_l*k60;
    w[44] = 1.0*VIII_l*k61*mIIa_l - 1.0*VIII_mIIa_l*k62;
    w[45] = 1.0*VIII_mIIa_l*k63;
    w[46] = -1.0*IIa_TM_l*k65 + 1.0*IIa_f*TM_l*k64;
    w[47] = -1.0*IIa_TM_PC_l*k67 + 1.0*IIa_TM_l*PC_l*k66;
    w[48] = 1.0*IIa_TM_PC_l*k68;
    w[49] = 1.0*AT_f*k69*mIIa_f;
    w[50] = 1.0*IX_l*TF_VIIa_l*k5 - 1.0*TF_VIIa_IX_l*k6;
    w[51] = -1.0*APC_PS_l*k71 + 1.0*APC_l*PS_l*k70;
    w[52] = 1.0*IX_l*XIa_l*k72 - 1.0*XIa_IX_l*k73;
    w[53] = 1.0*XIa_IX_l*k74;
    w[54] = 1.0*AT_f*XIa_l*k76;
    w[55] = 1.0*IIa_f*alpha2M_l*k77;
    w[56] = 1.0*Xa_f*alpha2M_l*k78;
    w[57] = 1.0*TF_VIIa_IX_l*k7;
    w[58] = -1.0*TF_VIIa_X_l*k9 + 1.0*TF_VIIa_l*X_l*k8;
    w[59] = 1.0*TF_VIIa_X_l*k10;
    w[60] = 1.0*TF_VIIa_Xa_l*k75;
    w[61] = -1.0*TF_VII_Xa_l*k12 + 1.0*TF_VII_l*Xa_l*k11;
    w[62] = 1.0*TF_VII_Xa_l*k13;
    w[63] = -1.0*IXa_VIIIa_l*k15 + 1.0*IXa_l*VIIIa_l*k14;
    w[64] = 1.0*Va_l*Xa_l*k16 - 1.0*Xa_Va_l*k17;
    w[65] = -1.0*IXa_VIIIa_X_l*k19 + 1.0*IXa_VIIIa_l*X_l*k18;
    w[66] = 1.0*IXa_VIIIa_X_l*k20;
    w[67] = -1.0*V_Xa_l*k22 + 1.0*V_l*Xa_l*k21;
    w[68] = 1.0*V_Xa_l*k23;
}