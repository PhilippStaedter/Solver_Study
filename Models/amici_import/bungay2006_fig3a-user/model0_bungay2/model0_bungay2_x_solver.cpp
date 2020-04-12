#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_bungay2(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = APC_PS_Va_l;
    x_solver[1] = APC_PS_l;
    x_solver[2] = APC_f;
    x_solver[3] = APC_l;
    x_solver[4] = AT_f;
    x_solver[5] = II_f;
    x_solver[6] = II_l;
    x_solver[7] = IIa_AT_f;
    x_solver[8] = IIa_TM_PC_l;
    x_solver[9] = IIa_TM_l;
    x_solver[10] = IIa_f;
    x_solver[11] = LIPID;
    x_solver[12] = PC_f;
    x_solver[13] = PC_l;
    x_solver[14] = PS_f;
    x_solver[15] = PS_l;
    x_solver[16] = TFPI_Xa_TF_VIIa_l;
    x_solver[17] = TFPI_Xa_l;
    x_solver[18] = TFPI_f;
    x_solver[19] = TF_VII_Xa_l;
    x_solver[20] = TF_VII_l;
    x_solver[21] = TF_VIIa_X_l;
    x_solver[22] = TF_VIIa_Xa_l;
    x_solver[23] = TF_VIIa_l;
    x_solver[24] = TF_l;
    x_solver[25] = TM_l;
    x_solver[26] = VII_Xa_l;
    x_solver[27] = VII_f;
    x_solver[28] = VII_l;
    x_solver[29] = VIIa_f;
    x_solver[30] = VIIa_l;
    x_solver[31] = V_IIa_l;
    x_solver[32] = V_Xa_l;
    x_solver[33] = V_f;
    x_solver[34] = V_l;
    x_solver[35] = V_mIIa_l;
    x_solver[36] = Va_f;
    x_solver[37] = Va_l;
    x_solver[38] = Vai_f;
    x_solver[39] = Vai_l;
    x_solver[40] = X_f;
    x_solver[41] = X_l;
    x_solver[42] = Xa_AT_f;
    x_solver[43] = Xa_Va_II_l;
    x_solver[44] = Xa_Va_l;
    x_solver[45] = Xa_Va_mIIa_l;
    x_solver[46] = Xa_f;
    x_solver[47] = Xa_l;
    x_solver[48] = alpha2M_IIa_l;
    x_solver[49] = alpha2M_Xa_l;
    x_solver[50] = alpha2M_l;
    x_solver[51] = mIIa_AT_l;
    x_solver[52] = mIIa_f;
    x_solver[53] = mIIa_l;
}