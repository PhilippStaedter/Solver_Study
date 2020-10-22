#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Bungay2006a(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = II_f;
    y[1] = II_l;
    y[2] = mIIa_f;
    y[3] = mIIa_l;
    y[4] = V_f;
    y[5] = V_l;
    y[6] = Va_f;
    y[7] = Va_l;
    y[8] = VII_f;
    y[9] = VII_l;
    y[10] = VIIa_f;
    y[11] = VIIa_l;
    y[12] = X_f;
    y[13] = X_l;
    y[14] = Xa_f;
    y[15] = Xa_l;
    y[16] = APC_f;
    y[17] = APC_l;
    y[18] = PS_f;
    y[19] = PS_l;
    y[20] = Vai_f;
    y[21] = Vai_l;
    y[22] = PC_f;
    y[23] = PC_l;
    y[24] = TF_l;
    y[25] = TF_VIIa_l;
    y[26] = TF_VII_l;
    y[27] = TF_VIIa_X_l;
    y[28] = TF_VIIa_Xa_l;
    y[29] = TF_VII_Xa_l;
    y[30] = Xa_Va_l;
    y[31] = V_Xa_l;
    y[32] = IIa_f;
    y[33] = V_IIa_l;
    y[34] = Xa_Va_II_l;
    y[35] = Xa_Va_mIIa_l;
    y[36] = APC_PS_l;
    y[37] = TFPI_f;
    y[38] = AT_f;
    y[39] = IIa_AT_f;
    y[40] = TFPI_Xa_l;
    y[41] = TFPI_Xa_TF_VIIa_l;
    y[42] = APC_PS_Va_l;
    y[43] = Xa_AT_f;
    y[44] = VII_Xa_l;
    y[45] = V_mIIa_l;
    y[46] = TM_l;
    y[47] = IIa_TM_l;
    y[48] = IIa_TM_PC_l;
    y[49] = mIIa_AT_l;
    y[50] = LIPID;
    y[51] = alpha2M_l;
    y[52] = alpha2M_IIa_l;
    y[53] = alpha2M_Xa_l;
}