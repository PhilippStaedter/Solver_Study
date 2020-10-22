#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_bungay1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = APC_PS_VIIIa_l;
    y[1] = APC_PS_Va_l;
    y[2] = APC_PS_l;
    y[3] = APC_f;
    y[4] = APC_l;
    y[5] = AT_XIa_l;
    y[6] = AT_f;
    y[7] = II_f;
    y[8] = II_l;
    y[9] = IIa_AT_f;
    y[10] = IIa_TM_PC_l;
    y[11] = IIa_TM_l;
    y[12] = IIa_f;
    y[13] = IX_f;
    y[14] = IX_l;
    y[15] = IXa_AT_f;
    y[16] = IXa_VIIIa_X_l;
    y[17] = IXa_VIIIa_l;
    y[18] = IXa_f;
    y[19] = IXa_l;
    y[20] = LIPID;
    y[21] = PC_f;
    y[22] = PC_l;
    y[23] = PS_f;
    y[24] = PS_l;
    y[25] = TFPI_Xa_TF_VIIa_l;
    y[26] = TFPI_Xa_l;
    y[27] = TFPI_f;
    y[28] = TF_VII_Xa_l;
    y[29] = TF_VII_l;
    y[30] = TF_VIIa_IX_l;
    y[31] = TF_VIIa_IXa_l;
    y[32] = TF_VIIa_X_l;
    y[33] = TF_VIIa_Xa_l;
    y[34] = TF_VIIa_l;
    y[35] = TF_l;
    y[36] = TM_l;
    y[37] = VIII_IIa_l;
    y[38] = VIII_Xa_l;
    y[39] = VIII_f;
    y[40] = VIII_l;
    y[41] = VIII_mIIa_l;
    y[42] = VIIIa_f;
    y[43] = VIIIa_l;
    y[44] = VIIIai_f;
    y[45] = VIIIai_l;
    y[46] = VII_Xa_l;
    y[47] = VII_f;
    y[48] = VII_l;
    y[49] = VIIa_f;
    y[50] = VIIa_l;
    y[51] = V_IIa_l;
    y[52] = V_Xa_l;
    y[53] = V_f;
    y[54] = V_l;
    y[55] = V_mIIa_l;
    y[56] = Va_f;
    y[57] = Va_l;
    y[58] = Vai_f;
    y[59] = Vai_l;
    y[60] = XI_IIa_l;
    y[61] = XI_f;
    y[62] = XIa_IX_l;
    y[63] = XIa_l;
    y[64] = X_f;
    y[65] = X_l;
    y[66] = Xa_AT_f;
    y[67] = Xa_Va_II_l;
    y[68] = Xa_Va_l;
    y[69] = Xa_Va_mIIa_l;
    y[70] = Xa_f;
    y[71] = Xa_l;
    y[72] = alpha2M_IIa_l;
    y[73] = alpha2M_Xa_l;
    y[74] = alpha2M_l;
    y[75] = mIIa_AT_l;
    y[76] = mIIa_f;
    y[77] = mIIa_l;
}