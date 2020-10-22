#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_Bungay2006a(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*LIPID*konII/nva;
    dwdx[1] = -1.0*koffII;
    dwdx[2] = 1.0*Xa_Va_l*k33;
    dwdx[3] = 1.0*LIPID*konmIIa/nva;
    dwdx[4] = 1.0*AT_f*k69;
    dwdx[5] = -1.0*koffmIIa;
    dwdx[6] = 1.0*Xa_Va_l*k35;
    dwdx[7] = 1.0*V_l*k58;
    dwdx[8] = 1.0*LIPID*konV/nva;
    dwdx[9] = -1.0*koffV;
    dwdx[10] = 1.0*Xa_l*k21;
    dwdx[11] = 1.0*IIa_f*k27;
    dwdx[12] = 1.0*k58*mIIa_l;
    dwdx[13] = 1.0*LIPID*konVa/nva;
    dwdx[14] = -1.0*koffVa;
    dwdx[15] = 1.0*Xa_l*k16;
    dwdx[16] = 1.0*APC_PS_l*k48;
    dwdx[17] = 1.0*LIPID*konVII/nva;
    dwdx[18] = -1.0*koffVII;
    dwdx[19] = 1.0*TF_l*k3;
    dwdx[20] = 1.0*Xa_l*k39;
    dwdx[21] = 1.0*LIPID*konVIIa/nva;
    dwdx[22] = -1.0*koffVIIa;
    dwdx[23] = 1.0*TF_l*k1;
    dwdx[24] = 1.0*LIPID*konX/nva;
    dwdx[25] = -1.0*koffX;
    dwdx[26] = 1.0*TF_VIIa_l*k8;
    dwdx[27] = LIPID*konXa/nva;
    dwdx[28] = 1.0*TFPI_f*k51;
    dwdx[29] = 1.0*AT_f*k56;
    dwdx[30] = 1.0*alpha2M_l*k78;
    dwdx[31] = -koffXa;
    dwdx[32] = 1.0*TF_VII_l*k11;
    dwdx[33] = 1.0*Va_l*k16;
    dwdx[34] = 1.0*V_l*k21;
    dwdx[35] = 1.0*VII_l*k39;
    dwdx[36] = 1.0*LIPID*konAPC/nva;
    dwdx[37] = -1.0*koffAPC;
    dwdx[38] = 1.0*PS_l*k70;
    dwdx[39] = 1.0*LIPID*konPS/nva;
    dwdx[40] = -1.0*koffPS;
    dwdx[41] = 1.0*APC_l*k70;
    dwdx[42] = 1.0*LIPID*konVai/nva;
    dwdx[43] = -1.0*koffVai;
    dwdx[44] = 1.0*LIPID*konPC/nva;
    dwdx[45] = -1.0*koffPC;
    dwdx[46] = 1.0*IIa_TM_l*k66;
    dwdx[47] = 1.0*VIIa_l*k1;
    dwdx[48] = 1.0*VII_l*k3;
    dwdx[49] = -1.0*k2;
    dwdx[50] = 1.0*X_l*k8;
    dwdx[51] = 1.0*TFPI_Xa_l*k53;
    dwdx[52] = -1.0*k4;
    dwdx[53] = 1.0*Xa_l*k11;
    dwdx[54] = -1.0*k9;
    dwdx[55] = 1.0*k10;
    dwdx[56] = 1.0*k75;
    dwdx[57] = -1.0*k12;
    dwdx[58] = 1.0*k13;
    dwdx[59] = -1.0*k17;
    dwdx[60] = 1.0*II_l*k33;
    dwdx[61] = 1.0*k35*mIIa_l;
    dwdx[62] = -1.0*k22;
    dwdx[63] = 1.0*k23;
    dwdx[64] = 1.0*V_l*k27;
    dwdx[65] = 1.0*AT_f*k57;
    dwdx[66] = 1.0*TM_l*k64;
    dwdx[67] = 1.0*alpha2M_l*k77;
    dwdx[68] = -1.0*k28;
    dwdx[69] = 1.0*k29;
    dwdx[70] = -1.0*k34;
    dwdx[71] = 1.0*k37;
    dwdx[72] = -1.0*k36;
    dwdx[73] = 1.0*k38;
    dwdx[74] = 1.0*Va_l*k48;
    dwdx[75] = -1.0*k71;
    dwdx[76] = 1.0*Xa_f*k51;
    dwdx[77] = 1.0*Xa_f*k56;
    dwdx[78] = 1.0*IIa_f*k57;
    dwdx[79] = 1.0*k69*mIIa_f;
    dwdx[80] = -1.0*k52;
    dwdx[81] = 1.0*TF_VIIa_l*k53;
    dwdx[82] = -1.0*k54;
    dwdx[83] = -1.0*k49;
    dwdx[84] = 1.0*k50;
    dwdx[85] = -1.0*k40;
    dwdx[86] = 1.0*k41;
    dwdx[87] = -1.0*k59;
    dwdx[88] = 1.0*k60;
    dwdx[89] = 1.0*IIa_f*k64;
    dwdx[90] = -1.0*k65;
    dwdx[91] = 1.0*PC_l*k66;
    dwdx[92] = -1.0*k67;
    dwdx[93] = 1.0*k68;
    dwdx[94] = 1.0*II_f*konII/nva;
    dwdx[95] = 1.0*konmIIa*mIIa_f/nva;
    dwdx[96] = 1.0*V_f*konV/nva;
    dwdx[97] = 1.0*Va_f*konVa/nva;
    dwdx[98] = 1.0*VII_f*konVII/nva;
    dwdx[99] = 1.0*VIIa_f*konVIIa/nva;
    dwdx[100] = 1.0*X_f*konX/nva;
    dwdx[101] = Xa_f*konXa/nva;
    dwdx[102] = 1.0*APC_f*konAPC/nva;
    dwdx[103] = 1.0*PS_f*konPS/nva;
    dwdx[104] = 1.0*Vai_f*konVai/nva;
    dwdx[105] = 1.0*PC_f*konPC/nva;
    dwdx[106] = 1.0*IIa_f*k77;
    dwdx[107] = 1.0*Xa_f*k78;
}