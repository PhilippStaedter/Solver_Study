#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Bungay2006a(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = II_f;
    x_rdata[1] = II_l;
    x_rdata[2] = mIIa_f;
    x_rdata[3] = mIIa_l;
    x_rdata[4] = V_f;
    x_rdata[5] = V_l;
    x_rdata[6] = Va_f;
    x_rdata[7] = Va_l;
    x_rdata[8] = VII_f;
    x_rdata[9] = VII_l;
    x_rdata[10] = VIIa_f;
    x_rdata[11] = VIIa_l;
    x_rdata[12] = X_f;
    x_rdata[13] = X_l;
    x_rdata[14] = Xa_f;
    x_rdata[15] = Xa_l;
    x_rdata[16] = APC_f;
    x_rdata[17] = APC_l;
    x_rdata[18] = PS_f;
    x_rdata[19] = PS_l;
    x_rdata[20] = Vai_f;
    x_rdata[21] = Vai_l;
    x_rdata[22] = PC_f;
    x_rdata[23] = PC_l;
    x_rdata[24] = TF_l;
    x_rdata[25] = TF_VIIa_l;
    x_rdata[26] = TF_VII_l;
    x_rdata[27] = TF_VIIa_X_l;
    x_rdata[28] = TF_VIIa_Xa_l;
    x_rdata[29] = TF_VII_Xa_l;
    x_rdata[30] = Xa_Va_l;
    x_rdata[31] = V_Xa_l;
    x_rdata[32] = IIa_f;
    x_rdata[33] = V_IIa_l;
    x_rdata[34] = Xa_Va_II_l;
    x_rdata[35] = Xa_Va_mIIa_l;
    x_rdata[36] = APC_PS_l;
    x_rdata[37] = TFPI_f;
    x_rdata[38] = AT_f;
    x_rdata[39] = IIa_AT_f;
    x_rdata[40] = TFPI_Xa_l;
    x_rdata[41] = TFPI_Xa_TF_VIIa_l;
    x_rdata[42] = APC_PS_Va_l;
    x_rdata[43] = Xa_AT_f;
    x_rdata[44] = VII_Xa_l;
    x_rdata[45] = V_mIIa_l;
    x_rdata[46] = TM_l;
    x_rdata[47] = IIa_TM_l;
    x_rdata[48] = IIa_TM_PC_l;
    x_rdata[49] = mIIa_AT_l;
    x_rdata[50] = LIPID;
    x_rdata[51] = alpha2M_l;
    x_rdata[52] = alpha2M_IIa_l;
    x_rdata[53] = alpha2M_Xa_l;
}