#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_bungay1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = APC_PS_VIIIa_l;
    x_solver[1] = APC_PS_Va_l;
    x_solver[2] = APC_PS_l;
    x_solver[3] = APC_f;
    x_solver[4] = APC_l;
    x_solver[5] = AT_XIa_l;
    x_solver[6] = AT_f;
    x_solver[7] = II_f;
    x_solver[8] = II_l;
    x_solver[9] = IIa_AT_f;
    x_solver[10] = IIa_TM_PC_l;
    x_solver[11] = IIa_TM_l;
    x_solver[12] = IIa_f;
    x_solver[13] = IX_f;
    x_solver[14] = IX_l;
    x_solver[15] = IXa_AT_f;
    x_solver[16] = IXa_VIIIa_X_l;
    x_solver[17] = IXa_VIIIa_l;
    x_solver[18] = IXa_f;
    x_solver[19] = IXa_l;
    x_solver[20] = LIPID;
    x_solver[21] = PC_f;
    x_solver[22] = PC_l;
    x_solver[23] = PS_f;
    x_solver[24] = PS_l;
    x_solver[25] = TFPI_Xa_TF_VIIa_l;
    x_solver[26] = TFPI_Xa_l;
    x_solver[27] = TFPI_f;
    x_solver[28] = TF_VII_Xa_l;
    x_solver[29] = TF_VII_l;
    x_solver[30] = TF_VIIa_IX_l;
    x_solver[31] = TF_VIIa_IXa_l;
    x_solver[32] = TF_VIIa_X_l;
    x_solver[33] = TF_VIIa_Xa_l;
    x_solver[34] = TF_VIIa_l;
    x_solver[35] = TF_l;
    x_solver[36] = TM_l;
    x_solver[37] = VIII_IIa_l;
    x_solver[38] = VIII_Xa_l;
    x_solver[39] = VIII_f;
    x_solver[40] = VIII_l;
    x_solver[41] = VIII_mIIa_l;
    x_solver[42] = VIIIa_f;
    x_solver[43] = VIIIa_l;
    x_solver[44] = VIIIai_f;
    x_solver[45] = VIIIai_l;
    x_solver[46] = VII_Xa_l;
    x_solver[47] = VII_f;
    x_solver[48] = VII_l;
    x_solver[49] = VIIa_f;
    x_solver[50] = VIIa_l;
    x_solver[51] = V_IIa_l;
    x_solver[52] = V_Xa_l;
    x_solver[53] = V_f;
    x_solver[54] = V_l;
    x_solver[55] = V_mIIa_l;
    x_solver[56] = Va_f;
    x_solver[57] = Va_l;
    x_solver[58] = Vai_f;
    x_solver[59] = Vai_l;
    x_solver[60] = XI_IIa_l;
    x_solver[61] = XI_f;
    x_solver[62] = XIa_IX_l;
    x_solver[63] = XIa_l;
    x_solver[64] = X_f;
    x_solver[65] = X_l;
    x_solver[66] = Xa_AT_f;
    x_solver[67] = Xa_Va_II_l;
    x_solver[68] = Xa_Va_l;
    x_solver[69] = Xa_Va_mIIa_l;
    x_solver[70] = Xa_f;
    x_solver[71] = Xa_l;
    x_solver[72] = alpha2M_IIa_l;
    x_solver[73] = alpha2M_Xa_l;
    x_solver[74] = alpha2M_l;
    x_solver[75] = mIIa_AT_l;
    x_solver[76] = mIIa_f;
    x_solver[77] = mIIa_l;
}