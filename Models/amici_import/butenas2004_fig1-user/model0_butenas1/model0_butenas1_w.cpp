#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_butenas1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*TF*VII*k2 - 1.0*TF_VII*k1;
    w[1] = 1.0*IIa*VIII*k17;
    w[2] = 1.0*IXa*VIIIa*k19 - 1.0*IXa_VIIIa*k18;
    w[3] = 1.0*IXa_VIIIa*X*k21 - 1.0*IXa_VIIIa_X*k20;
    w[4] = 1.0*IXa_VIIIa_X*k22;
    w[5] = 1.0*VIIIa*k24 - 1.0*VIIIa1_L*VIIIa2*k23;
    w[6] = 1.0*IXa_VIIIa_X*k25;
    w[7] = 1.0*IXa_VIIIa*k25;
    w[8] = 1.0*IIa*V*k26;
    w[9] = 1.0*Va*Xa*k28 - 1.0*Xa_Va*k27;
    w[10] = 1.0*II*Xa_Va*k30 - 1.0*Xa_Va_II*k29;
    w[11] = 1.0*Xa_Va_II*k31;
    w[12] = 1.0*Xa_Va*k32*mIIa;
    w[13] = 1.0*TF*VIIa*k4 - 1.0*TF_VIIa*k3;
    w[14] = 1.0*TFPI*Xa*k34 - 1.0*Xa_TFPI*k33;
    w[15] = 1.0*TFPI*TF_VIIa_Xa*k36 - 1.0*TF_VIIa_Xa_TFPI*k35;
    w[16] = 1.0*TF_VIIa*Xa_TFPI*k37;
    w[17] = 1.0*ATIII*Xa*k38;
    w[18] = 1.0*ATIII*k39*mIIa;
    w[19] = 1.0*ATIII*IXa*k40;
    w[20] = 1.0*ATIII*IIa*k41;
    w[21] = 1.0*ATIII*TF_VIIa*k42;
    w[22] = 1.0*IXa*X*k43;
    w[23] = 1.0*V*k44*mIIa;
    w[24] = 1.0*TF_VIIa*VII*k5;
    w[25] = 1.0*VII*Xa*k6;
    w[26] = 1.0*IIa*VII*k7;
    w[27] = 1.0*TF_VIIa*X*k9 - 1.0*TF_VIIa_X*k8;
    w[28] = 1.0*TF_VIIa_X*k10;
    w[29] = 1.0*TF_VIIa*Xa*k12 - 1.0*TF_VIIa_Xa*k11;
    w[30] = 1.0*IX*TF_VIIa*k14 - 1.0*TF_VIIa_IX*k13;
    w[31] = 1.0*TF_VIIa_IX*k15;
    w[32] = 1.0*II*Xa*k16;
}