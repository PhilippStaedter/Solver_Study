#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_bungay2(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = APC_PS_Va_l;
    x_rdata[1] = APC_PS_l;
    x_rdata[2] = APC_f;
    x_rdata[3] = APC_l;
    x_rdata[4] = AT_f;
    x_rdata[5] = II_f;
    x_rdata[6] = II_l;
    x_rdata[7] = IIa_AT_f;
    x_rdata[8] = IIa_TM_PC_l;
    x_rdata[9] = IIa_TM_l;
    x_rdata[10] = IIa_f;
    x_rdata[11] = LIPID;
    x_rdata[12] = PC_f;
    x_rdata[13] = PC_l;
    x_rdata[14] = PS_f;
    x_rdata[15] = PS_l;
    x_rdata[16] = TFPI_Xa_TF_VIIa_l;
    x_rdata[17] = TFPI_Xa_l;
    x_rdata[18] = TFPI_f;
    x_rdata[19] = TF_VII_Xa_l;
    x_rdata[20] = TF_VII_l;
    x_rdata[21] = TF_VIIa_X_l;
    x_rdata[22] = TF_VIIa_Xa_l;
    x_rdata[23] = TF_VIIa_l;
    x_rdata[24] = TF_l;
    x_rdata[25] = TM_l;
    x_rdata[26] = VII_Xa_l;
    x_rdata[27] = VII_f;
    x_rdata[28] = VII_l;
    x_rdata[29] = VIIa_f;
    x_rdata[30] = VIIa_l;
    x_rdata[31] = V_IIa_l;
    x_rdata[32] = V_Xa_l;
    x_rdata[33] = V_f;
    x_rdata[34] = V_l;
    x_rdata[35] = V_mIIa_l;
    x_rdata[36] = Va_f;
    x_rdata[37] = Va_l;
    x_rdata[38] = Vai_f;
    x_rdata[39] = Vai_l;
    x_rdata[40] = X_f;
    x_rdata[41] = X_l;
    x_rdata[42] = Xa_AT_f;
    x_rdata[43] = Xa_Va_II_l;
    x_rdata[44] = Xa_Va_l;
    x_rdata[45] = Xa_Va_mIIa_l;
    x_rdata[46] = Xa_f;
    x_rdata[47] = Xa_l;
    x_rdata[48] = alpha2M_IIa_l;
    x_rdata[49] = alpha2M_Xa_l;
    x_rdata[50] = alpha2M_l;
    x_rdata[51] = mIIa_AT_l;
    x_rdata[52] = mIIa_f;
    x_rdata[53] = mIIa_l;
}