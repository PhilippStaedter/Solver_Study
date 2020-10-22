#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model2_bungay3(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = APC_PS_VIIIa_l;
    x_rdata[1] = APC_PS_Va_l;
    x_rdata[2] = APC_PS_l;
    x_rdata[3] = APC_f;
    x_rdata[4] = APC_l;
    x_rdata[5] = AT_f;
    x_rdata[6] = II_f;
    x_rdata[7] = II_l;
    x_rdata[8] = IIa_AT_f;
    x_rdata[9] = IIa_TM_PC_l;
    x_rdata[10] = IIa_TM_l;
    x_rdata[11] = IIa_f;
    x_rdata[12] = IX_f;
    x_rdata[13] = IX_l;
    x_rdata[14] = IXa_AT_f;
    x_rdata[15] = IXa_VIIIa_X_l;
    x_rdata[16] = IXa_VIIIa_l;
    x_rdata[17] = IXa_f;
    x_rdata[18] = IXa_l;
    x_rdata[19] = LIPID;
    x_rdata[20] = PC_f;
    x_rdata[21] = PC_l;
    x_rdata[22] = PS_f;
    x_rdata[23] = PS_l;
    x_rdata[24] = TFPI_Xa_TF_VIIa_l;
    x_rdata[25] = TFPI_Xa_l;
    x_rdata[26] = TFPI_f;
    x_rdata[27] = TF_VII_Xa_l;
    x_rdata[28] = TF_VII_l;
    x_rdata[29] = TF_VIIa_IX_l;
    x_rdata[30] = TF_VIIa_IXa_l;
    x_rdata[31] = TF_VIIa_X_l;
    x_rdata[32] = TF_VIIa_Xa_l;
    x_rdata[33] = TF_VIIa_l;
    x_rdata[34] = TF_l;
    x_rdata[35] = TM_l;
    x_rdata[36] = VIII_IIa_l;
    x_rdata[37] = VIII_Xa_l;
    x_rdata[38] = VIII_f;
    x_rdata[39] = VIII_l;
    x_rdata[40] = VIII_mIIa_l;
    x_rdata[41] = VIIIa_f;
    x_rdata[42] = VIIIa_l;
    x_rdata[43] = VIIIai_f;
    x_rdata[44] = VIIIai_l;
    x_rdata[45] = VII_Xa_l;
    x_rdata[46] = VII_f;
    x_rdata[47] = VII_l;
    x_rdata[48] = VIIa_f;
    x_rdata[49] = VIIa_l;
    x_rdata[50] = V_IIa_l;
    x_rdata[51] = V_Xa_l;
    x_rdata[52] = V_f;
    x_rdata[53] = V_l;
    x_rdata[54] = V_mIIa_l;
    x_rdata[55] = Va_f;
    x_rdata[56] = Va_l;
    x_rdata[57] = Vai_f;
    x_rdata[58] = Vai_l;
    x_rdata[59] = XI_IIa_l;
    x_rdata[60] = XI_f;
    x_rdata[61] = XIa_IX_l;
    x_rdata[62] = XIa_l;
    x_rdata[63] = X_f;
    x_rdata[64] = X_l;
    x_rdata[65] = Xa_AT_f;
    x_rdata[66] = Xa_Va_II_l;
    x_rdata[67] = Xa_Va_l;
    x_rdata[68] = Xa_Va_mIIa_l;
    x_rdata[69] = Xa_f;
    x_rdata[70] = Xa_l;
    x_rdata[71] = mIIa_AT_l;
    x_rdata[72] = mIIa_f;
    x_rdata[73] = mIIa_l;
}