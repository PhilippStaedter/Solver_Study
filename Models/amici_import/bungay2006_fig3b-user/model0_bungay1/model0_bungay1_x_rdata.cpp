#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_bungay1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = APC_PS_VIIIa_l;
    x_rdata[1] = APC_PS_Va_l;
    x_rdata[2] = APC_PS_l;
    x_rdata[3] = APC_f;
    x_rdata[4] = APC_l;
    x_rdata[5] = AT_XIa_l;
    x_rdata[6] = AT_f;
    x_rdata[7] = II_f;
    x_rdata[8] = II_l;
    x_rdata[9] = IIa_AT_f;
    x_rdata[10] = IIa_TM_PC_l;
    x_rdata[11] = IIa_TM_l;
    x_rdata[12] = IIa_f;
    x_rdata[13] = IX_f;
    x_rdata[14] = IX_l;
    x_rdata[15] = IXa_AT_f;
    x_rdata[16] = IXa_VIIIa_X_l;
    x_rdata[17] = IXa_VIIIa_l;
    x_rdata[18] = IXa_f;
    x_rdata[19] = IXa_l;
    x_rdata[20] = LIPID;
    x_rdata[21] = PC_f;
    x_rdata[22] = PC_l;
    x_rdata[23] = PS_f;
    x_rdata[24] = PS_l;
    x_rdata[25] = TFPI_Xa_TF_VIIa_l;
    x_rdata[26] = TFPI_Xa_l;
    x_rdata[27] = TFPI_f;
    x_rdata[28] = TF_VII_Xa_l;
    x_rdata[29] = TF_VII_l;
    x_rdata[30] = TF_VIIa_IX_l;
    x_rdata[31] = TF_VIIa_IXa_l;
    x_rdata[32] = TF_VIIa_X_l;
    x_rdata[33] = TF_VIIa_Xa_l;
    x_rdata[34] = TF_VIIa_l;
    x_rdata[35] = TF_l;
    x_rdata[36] = TM_l;
    x_rdata[37] = VIII_IIa_l;
    x_rdata[38] = VIII_Xa_l;
    x_rdata[39] = VIII_f;
    x_rdata[40] = VIII_l;
    x_rdata[41] = VIII_mIIa_l;
    x_rdata[42] = VIIIa_f;
    x_rdata[43] = VIIIa_l;
    x_rdata[44] = VIIIai_f;
    x_rdata[45] = VIIIai_l;
    x_rdata[46] = VII_Xa_l;
    x_rdata[47] = VII_f;
    x_rdata[48] = VII_l;
    x_rdata[49] = VIIa_f;
    x_rdata[50] = VIIa_l;
    x_rdata[51] = V_IIa_l;
    x_rdata[52] = V_Xa_l;
    x_rdata[53] = V_f;
    x_rdata[54] = V_l;
    x_rdata[55] = V_mIIa_l;
    x_rdata[56] = Va_f;
    x_rdata[57] = Va_l;
    x_rdata[58] = Vai_f;
    x_rdata[59] = Vai_l;
    x_rdata[60] = XI_IIa_l;
    x_rdata[61] = XI_f;
    x_rdata[62] = XIa_IX_l;
    x_rdata[63] = XIa_l;
    x_rdata[64] = X_f;
    x_rdata[65] = X_l;
    x_rdata[66] = Xa_AT_f;
    x_rdata[67] = Xa_Va_II_l;
    x_rdata[68] = Xa_Va_l;
    x_rdata[69] = Xa_Va_mIIa_l;
    x_rdata[70] = Xa_f;
    x_rdata[71] = Xa_l;
    x_rdata[72] = alpha2M_IIa_l;
    x_rdata[73] = alpha2M_Xa_l;
    x_rdata[74] = alpha2M_l;
    x_rdata[75] = mIIa_AT_l;
    x_rdata[76] = mIIa_f;
    x_rdata[77] = mIIa_l;
}