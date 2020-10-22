#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Bungay2003(realtype *x_rdata, const realtype *x, const realtype *tcl){
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
    x_rdata[12] = VIII_f;
    x_rdata[13] = VIII_l;
    x_rdata[14] = VIIIa_f;
    x_rdata[15] = VIIIa_l;
    x_rdata[16] = IX_f;
    x_rdata[17] = IX_l;
    x_rdata[18] = IXa_f;
    x_rdata[19] = IXa_l;
    x_rdata[20] = X_f;
    x_rdata[21] = X_l;
    x_rdata[22] = Xa_f;
    x_rdata[23] = Xa_l;
    x_rdata[24] = APC_f;
    x_rdata[25] = APC_l;
    x_rdata[26] = PS_f;
    x_rdata[27] = PS_l;
    x_rdata[28] = VIIIai_f;
    x_rdata[29] = VIIIai_l;
    x_rdata[30] = Vai_f;
    x_rdata[31] = Vai_l;
    x_rdata[32] = PC_f;
    x_rdata[33] = PC_l;
    x_rdata[34] = TF_l;
    x_rdata[35] = TF_VIIa_l;
    x_rdata[36] = TF_VII_l;
    x_rdata[37] = TF_VIIa_IX_l;
    x_rdata[38] = TF_VIIa_IXa_l;
    x_rdata[39] = TF_VIIa_X_l;
    x_rdata[40] = TF_VIIa_Xa_l;
    x_rdata[41] = TF_VII_Xa_l;
    x_rdata[42] = IXa_VIIIa_l;
    x_rdata[43] = Xa_Va_l;
    x_rdata[44] = IXa_VIIIa_X_l;
    x_rdata[45] = V_Xa_l;
    x_rdata[46] = VIII_Xa_l;
    x_rdata[47] = IIa_f;
    x_rdata[48] = V_IIa_l;
    x_rdata[49] = VIII_IIa_l;
    x_rdata[50] = Xa_Va_II_l;
    x_rdata[51] = Xa_Va_mIIa_l;
    x_rdata[52] = XI_f;
    x_rdata[53] = XI_IIa_l;
    x_rdata[54] = XIa_l;
    x_rdata[55] = APC_PS_l;
    x_rdata[56] = APC_PS_VIIIa_l;
    x_rdata[57] = TFPI_f;
    x_rdata[58] = AT_f;
    x_rdata[59] = IIa_AT_f;
    x_rdata[60] = TFPI_Xa_l;
    x_rdata[61] = TFPI_Xa_TF_VIIa_l;
    x_rdata[62] = APC_PS_Va_l;
    x_rdata[63] = IXa_AT_f;
    x_rdata[64] = Xa_AT_f;
    x_rdata[65] = VII_Xa_l;
    x_rdata[66] = V_mIIa_l;
    x_rdata[67] = VIII_mIIa_l;
    x_rdata[68] = TM_l;
    x_rdata[69] = IIa_TM_l;
    x_rdata[70] = IIa_TM_PC_l;
    x_rdata[71] = mIIa_AT_l;
    x_rdata[72] = XIa_IX_l;
    x_rdata[73] = LIPID;
}