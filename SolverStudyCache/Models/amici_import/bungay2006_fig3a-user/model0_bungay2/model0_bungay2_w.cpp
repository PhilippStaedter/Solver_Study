#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_bungay2(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*II_f*LIPID*konII/nva - 1.0*II_l*koffII;
    w[1] = 1.0*LIPID*X_f*konX/nva - 1.0*X_l*koffX;
    w[2] = LIPID*Xa_f*konXa/nva - Xa_l*koffXa;
    w[3] = 1.0*APC_f*LIPID*konAPC/nva - 1.0*APC_l*koffAPC;
    w[4] = 1.0*LIPID*PS_f*konPS/nva - 1.0*PS_l*koffPS;
    w[5] = 1.0*LIPID*Vai_f*konVai/nva - 1.0*Vai_l*koffVai;
    w[6] = 1.0*LIPID*PC_f*konPC/nva - 1.0*PC_l*koffPC;
    w[7] = 1.0*LIPID*konmIIa*mIIa_f/nva - 1.0*koffmIIa*mIIa_l;
    w[8] = 1.0*LIPID*V_f*konV/nva - 1.0*V_l*koffV;
    w[9] = 1.0*LIPID*Va_f*konVa/nva - 1.0*Va_l*koffVa;
    w[10] = 1.0*LIPID*VII_f*konVII/nva - 1.0*VII_l*koffVII;
    w[11] = 1.0*LIPID*VIIa_f*konVIIa/nva - 1.0*VIIa_l*koffVIIa;
    w[12] = -1.0*TF_VIIa_l*k2 + 1.0*TF_l*VIIa_l*k1;
    w[13] = 1.0*IIa_f*V_l*k27 - 1.0*V_IIa_l*k28;
    w[14] = 1.0*V_IIa_l*k29;
    w[15] = 1.0*II_l*Xa_Va_l*k33 - 1.0*Xa_Va_II_l*k34;
    w[16] = 1.0*Xa_Va_l*k35*mIIa_l - 1.0*Xa_Va_mIIa_l*k36;
    w[17] = 1.0*Xa_Va_II_l*k37;
    w[18] = 1.0*Xa_Va_mIIa_l*k38;
    w[19] = -1.0*VII_Xa_l*k40 + 1.0*VII_l*Xa_l*k39;
    w[20] = 1.0*VII_Xa_l*k41;
    w[21] = -1.0*APC_PS_Va_l*k49 + 1.0*APC_PS_l*Va_l*k48;
    w[22] = 1.0*APC_PS_Va_l*k50;
    w[23] = -1.0*TF_VII_l*k4 + 1.0*TF_l*VII_l*k3;
    w[24] = -1.0*TFPI_Xa_l*k52 + 1.0*TFPI_f*Xa_f*k51;
    w[25] = -1.0*TFPI_Xa_TF_VIIa_l*k54 + 1.0*TFPI_Xa_l*TF_VIIa_l*k53;
    w[26] = 1.0*AT_f*Xa_f*k56;
    w[27] = 1.0*AT_f*IIa_f*k57;
    w[28] = 1.0*V_l*k58*mIIa_l - 1.0*V_mIIa_l*k59;
    w[29] = 1.0*V_mIIa_l*k60;
    w[30] = -1.0*IIa_TM_l*k65 + 1.0*IIa_f*TM_l*k64;
    w[31] = -1.0*IIa_TM_PC_l*k67 + 1.0*IIa_TM_l*PC_l*k66;
    w[32] = 1.0*IIa_TM_PC_l*k68;
    w[33] = 1.0*AT_f*k69*mIIa_f;
    w[34] = -1.0*APC_PS_l*k71 + 1.0*APC_l*PS_l*k70;
    w[35] = 1.0*IIa_f*alpha2M_l*k77;
    w[36] = 1.0*Xa_f*alpha2M_l*k78;
    w[37] = -1.0*TF_VIIa_X_l*k9 + 1.0*TF_VIIa_l*X_l*k8;
    w[38] = 1.0*TF_VIIa_X_l*k10;
    w[39] = 1.0*TF_VIIa_Xa_l*k75;
    w[40] = -1.0*TF_VII_Xa_l*k12 + 1.0*TF_VII_l*Xa_l*k11;
    w[41] = 1.0*TF_VII_Xa_l*k13;
    w[42] = 1.0*Va_l*Xa_l*k16 - 1.0*Xa_Va_l*k17;
    w[43] = -1.0*V_Xa_l*k22 + 1.0*V_l*Xa_l*k21;
    w[44] = 1.0*V_Xa_l*k23;
}