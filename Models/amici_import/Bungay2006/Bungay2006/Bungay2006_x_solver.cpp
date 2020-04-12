#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Bungay2006(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = II_f;
    x_solver[1] = II_l;
    x_solver[2] = mIIa_f;
    x_solver[3] = mIIa_l;
    x_solver[4] = V_f;
    x_solver[5] = V_l;
    x_solver[6] = Va_f;
    x_solver[7] = Va_l;
    x_solver[8] = VII_f;
    x_solver[9] = VII_l;
    x_solver[10] = VIIa_f;
    x_solver[11] = VIIa_l;
    x_solver[12] = VIII_f;
    x_solver[13] = VIII_l;
    x_solver[14] = VIIIa_f;
    x_solver[15] = VIIIa_l;
    x_solver[16] = IX_f;
    x_solver[17] = IX_l;
    x_solver[18] = IXa_f;
    x_solver[19] = IXa_l;
    x_solver[20] = X_f;
    x_solver[21] = X_l;
    x_solver[22] = Xa_f;
    x_solver[23] = Xa_l;
    x_solver[24] = APC_f;
    x_solver[25] = APC_l;
    x_solver[26] = PS_f;
    x_solver[27] = PS_l;
    x_solver[28] = VIIIai_f;
    x_solver[29] = VIIIai_l;
    x_solver[30] = Vai_f;
    x_solver[31] = Vai_l;
    x_solver[32] = PC_f;
    x_solver[33] = PC_l;
    x_solver[34] = TF_l;
    x_solver[35] = TF_VIIa_l;
    x_solver[36] = TF_VII_l;
    x_solver[37] = TF_VIIa_IX_l;
    x_solver[38] = TF_VIIa_IXa_l;
    x_solver[39] = TF_VIIa_X_l;
    x_solver[40] = TF_VIIa_Xa_l;
    x_solver[41] = TF_VII_Xa_l;
    x_solver[42] = IXa_VIIIa_l;
    x_solver[43] = Xa_Va_l;
    x_solver[44] = IXa_VIIIa_X_l;
    x_solver[45] = V_Xa_l;
    x_solver[46] = VIII_Xa_l;
    x_solver[47] = IIa_f;
    x_solver[48] = V_IIa_l;
    x_solver[49] = VIII_IIa_l;
    x_solver[50] = Xa_Va_II_l;
    x_solver[51] = Xa_Va_mIIa_l;
    x_solver[52] = XI_f;
    x_solver[53] = XI_IIa_l;
    x_solver[54] = XIa_l;
    x_solver[55] = APC_PS_l;
    x_solver[56] = APC_PS_VIIIa_l;
    x_solver[57] = TFPI_f;
    x_solver[58] = AT_f;
    x_solver[59] = IIa_AT_f;
    x_solver[60] = TFPI_Xa_l;
    x_solver[61] = TFPI_Xa_TF_VIIa_l;
    x_solver[62] = APC_PS_Va_l;
    x_solver[63] = IXa_AT_f;
    x_solver[64] = Xa_AT_f;
    x_solver[65] = VII_Xa_l;
    x_solver[66] = V_mIIa_l;
    x_solver[67] = VIII_mIIa_l;
    x_solver[68] = TM_l;
    x_solver[69] = IIa_TM_l;
    x_solver[70] = IIa_TM_PC_l;
    x_solver[71] = mIIa_AT_l;
    x_solver[72] = XIa_IX_l;
    x_solver[73] = LIPID;
    x_solver[74] = alpha2M_l;
    x_solver[75] = alpha2M_IIa_l;
    x_solver[76] = alpha2M_Xa_l;
    x_solver[77] = AT_XIa_l;
}