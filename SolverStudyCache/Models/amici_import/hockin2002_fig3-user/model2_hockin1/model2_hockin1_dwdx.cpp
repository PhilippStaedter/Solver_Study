#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model2_hockin1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = 1.0*Xa*k38;
    dwdx[1] = 1.0*k39*mIIa;
    dwdx[2] = 1.0*IXa*k40;
    dwdx[3] = 1.0*IIa*k41;
    dwdx[4] = 1.0*TF_VIIa*k42;
    dwdx[5] = 1.0*Xa_Va*k30;
    dwdx[6] = 1.0*Xa*k16;
    dwdx[7] = 1.0*VIII*k17;
    dwdx[8] = 1.0*V*k26;
    dwdx[9] = 1.0*ATIII*k41;
    dwdx[10] = 1.0*VII*k7;
    dwdx[11] = 1.0*TF_VIIa*k14;
    dwdx[12] = 1.0*VIIIa*k19;
    dwdx[13] = 1.0*ATIII*k40;
    dwdx[14] = -1.0*k18;
    dwdx[15] = 1.0*X*k21;
    dwdx[16] = 1.0*k25;
    dwdx[17] = -1.0*k20;
    dwdx[18] = 1.0*k22;
    dwdx[19] = 1.0*k25;
    dwdx[20] = 1.0*VII*k2;
    dwdx[21] = 1.0*VIIa*k4;
    dwdx[22] = 1.0*Xa*k34;
    dwdx[23] = 1.0*TF_VIIa_Xa*k36;
    dwdx[24] = -1.0*k1;
    dwdx[25] = -1.0*k3;
    dwdx[26] = 1.0*Xa_TFPI*k37;
    dwdx[27] = 1.0*ATIII*k42;
    dwdx[28] = 1.0*VII*k5;
    dwdx[29] = 1.0*X*k9;
    dwdx[30] = 1.0*Xa*k12;
    dwdx[31] = 1.0*IX*k14;
    dwdx[32] = -1.0*k13;
    dwdx[33] = 1.0*k15;
    dwdx[34] = -1.0*k8;
    dwdx[35] = 1.0*k10;
    dwdx[36] = 1.0*TFPI*k36;
    dwdx[37] = -1.0*k11;
    dwdx[38] = -1.0*k35;
    dwdx[39] = 1.0*IIa*k26;
    dwdx[40] = 1.0*TF*k2;
    dwdx[41] = 1.0*TF_VIIa*k5;
    dwdx[42] = 1.0*Xa*k6;
    dwdx[43] = 1.0*IIa*k7;
    dwdx[44] = 1.0*IIa*k17;
    dwdx[45] = 1.0*IXa*k19;
    dwdx[46] = 1.0*k24;
    dwdx[47] = -1.0*VIIIa2*k23;
    dwdx[48] = -1.0*VIIIa1_L*k23;
    dwdx[49] = 1.0*TF*k4;
    dwdx[50] = 1.0*Xa*k28;
    dwdx[51] = 1.0*IXa_VIIIa*k21;
    dwdx[52] = 1.0*TF_VIIa*k9;
    dwdx[53] = 1.0*Va*k28;
    dwdx[54] = 1.0*TFPI*k34;
    dwdx[55] = 1.0*ATIII*k38;
    dwdx[56] = 1.0*VII*k6;
    dwdx[57] = 1.0*TF_VIIa*k12;
    dwdx[58] = 1.0*II*k16;
    dwdx[59] = -1.0*k33;
    dwdx[60] = 1.0*TF_VIIa*k37;
    dwdx[61] = -1.0*k27;
    dwdx[62] = 1.0*II*k30;
    dwdx[63] = 1.0*k32*mIIa;
    dwdx[64] = -1.0*k29;
    dwdx[65] = 1.0*k31;
    dwdx[66] = 1.0*Xa_Va*k32;
    dwdx[67] = 1.0*ATIII*k39;
}