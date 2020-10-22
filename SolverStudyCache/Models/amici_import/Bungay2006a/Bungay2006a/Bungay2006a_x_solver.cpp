#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Bungay2006a(realtype *x_solver, const realtype *x_rdata){
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
    x_solver[12] = X_f;
    x_solver[13] = X_l;
    x_solver[14] = Xa_f;
    x_solver[15] = Xa_l;
    x_solver[16] = APC_f;
    x_solver[17] = APC_l;
    x_solver[18] = PS_f;
    x_solver[19] = PS_l;
    x_solver[20] = Vai_f;
    x_solver[21] = Vai_l;
    x_solver[22] = PC_f;
    x_solver[23] = PC_l;
    x_solver[24] = TF_l;
    x_solver[25] = TF_VIIa_l;
    x_solver[26] = TF_VII_l;
    x_solver[27] = TF_VIIa_X_l;
    x_solver[28] = TF_VIIa_Xa_l;
    x_solver[29] = TF_VII_Xa_l;
    x_solver[30] = Xa_Va_l;
    x_solver[31] = V_Xa_l;
    x_solver[32] = IIa_f;
    x_solver[33] = V_IIa_l;
    x_solver[34] = Xa_Va_II_l;
    x_solver[35] = Xa_Va_mIIa_l;
    x_solver[36] = APC_PS_l;
    x_solver[37] = TFPI_f;
    x_solver[38] = AT_f;
    x_solver[39] = IIa_AT_f;
    x_solver[40] = TFPI_Xa_l;
    x_solver[41] = TFPI_Xa_TF_VIIa_l;
    x_solver[42] = APC_PS_Va_l;
    x_solver[43] = Xa_AT_f;
    x_solver[44] = VII_Xa_l;
    x_solver[45] = V_mIIa_l;
    x_solver[46] = TM_l;
    x_solver[47] = IIa_TM_l;
    x_solver[48] = IIa_TM_PC_l;
    x_solver[49] = mIIa_AT_l;
    x_solver[50] = LIPID;
    x_solver[51] = alpha2M_l;
    x_solver[52] = alpha2M_IIa_l;
    x_solver[53] = alpha2M_Xa_l;
}