#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model1_hockin1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = ATIII;
    x_rdata[1] = II;
    x_rdata[2] = IIa;
    x_rdata[3] = IIa_ATIII;
    x_rdata[4] = IX;
    x_rdata[5] = IXa;
    x_rdata[6] = IXa_ATIII;
    x_rdata[7] = IXa_VIIIa;
    x_rdata[8] = IXa_VIIIa_X;
    x_rdata[9] = TF;
    x_rdata[10] = TFPI;
    x_rdata[11] = TF_VII;
    x_rdata[12] = TF_VIIa;
    x_rdata[13] = TF_VIIa_ATIII;
    x_rdata[14] = TF_VIIa_IX;
    x_rdata[15] = TF_VIIa_X;
    x_rdata[16] = TF_VIIa_Xa;
    x_rdata[17] = TF_VIIa_Xa_TFPI;
    x_rdata[18] = V;
    x_rdata[19] = VII;
    x_rdata[20] = VIII;
    x_rdata[21] = VIIIa;
    x_rdata[22] = VIIIa1_L;
    x_rdata[23] = VIIIa2;
    x_rdata[24] = VIIa;
    x_rdata[25] = Va;
    x_rdata[26] = X;
    x_rdata[27] = Xa;
    x_rdata[28] = Xa_ATIII;
    x_rdata[29] = Xa_TFPI;
    x_rdata[30] = Xa_Va;
    x_rdata[31] = Xa_Va_II;
    x_rdata[32] = mIIa;
    x_rdata[33] = mIIa_ATIII;
}