#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_bungay2(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = APC_PS_Va_l;
    y[1] = APC_PS_l;
    y[2] = APC_f;
    y[3] = APC_l;
    y[4] = AT_f;
    y[5] = II_f;
    y[6] = II_l;
    y[7] = IIa_AT_f;
    y[8] = IIa_TM_PC_l;
    y[9] = IIa_TM_l;
    y[10] = IIa_f;
    y[11] = LIPID;
    y[12] = PC_f;
    y[13] = PC_l;
    y[14] = PS_f;
    y[15] = PS_l;
    y[16] = TFPI_Xa_TF_VIIa_l;
    y[17] = TFPI_Xa_l;
    y[18] = TFPI_f;
    y[19] = TF_VII_Xa_l;
    y[20] = TF_VII_l;
    y[21] = TF_VIIa_X_l;
    y[22] = TF_VIIa_Xa_l;
    y[23] = TF_VIIa_l;
    y[24] = TF_l;
    y[25] = TM_l;
    y[26] = VII_Xa_l;
    y[27] = VII_f;
    y[28] = VII_l;
    y[29] = VIIa_f;
    y[30] = VIIa_l;
    y[31] = V_IIa_l;
    y[32] = V_Xa_l;
    y[33] = V_f;
    y[34] = V_l;
    y[35] = V_mIIa_l;
    y[36] = Va_f;
    y[37] = Va_l;
    y[38] = Vai_f;
    y[39] = Vai_l;
    y[40] = X_f;
    y[41] = X_l;
    y[42] = Xa_AT_f;
    y[43] = Xa_Va_II_l;
    y[44] = Xa_Va_l;
    y[45] = Xa_Va_mIIa_l;
    y[46] = Xa_f;
    y[47] = Xa_l;
    y[48] = alpha2M_IIa_l;
    y[49] = alpha2M_Xa_l;
    y[50] = alpha2M_l;
    y[51] = mIIa_AT_l;
    y[52] = mIIa_f;
    y[53] = mIIa_l;
}