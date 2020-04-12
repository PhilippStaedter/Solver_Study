#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Bungay2003(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
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
    y[12] = VIII_f;
    y[13] = VIII_l;
    y[14] = VIIIa_f;
    y[15] = VIIIa_l;
    y[16] = IX_f;
    y[17] = IX_l;
    y[18] = IXa_f;
    y[19] = IXa_l;
    y[20] = X_f;
    y[21] = X_l;
    y[22] = Xa_f;
    y[23] = Xa_l;
    y[24] = APC_f;
    y[25] = APC_l;
    y[26] = PS_f;
    y[27] = PS_l;
    y[28] = VIIIai_f;
    y[29] = VIIIai_l;
    y[30] = Vai_f;
    y[31] = Vai_l;
    y[32] = PC_f;
    y[33] = PC_l;
    y[34] = TF_l;
    y[35] = TF_VIIa_l;
    y[36] = TF_VII_l;
    y[37] = TF_VIIa_IX_l;
    y[38] = TF_VIIa_IXa_l;
    y[39] = TF_VIIa_X_l;
    y[40] = TF_VIIa_Xa_l;
    y[41] = TF_VII_Xa_l;
    y[42] = IXa_VIIIa_l;
    y[43] = Xa_Va_l;
    y[44] = IXa_VIIIa_X_l;
    y[45] = V_Xa_l;
    y[46] = VIII_Xa_l;
    y[47] = IIa_f;
    y[48] = V_IIa_l;
    y[49] = VIII_IIa_l;
    y[50] = Xa_Va_II_l;
    y[51] = Xa_Va_mIIa_l;
    y[52] = XI_f;
    y[53] = XI_IIa_l;
    y[54] = XIa_l;
    y[55] = APC_PS_l;
    y[56] = APC_PS_VIIIa_l;
    y[57] = TFPI_f;
    y[58] = AT_f;
    y[59] = IIa_AT_f;
    y[60] = TFPI_Xa_l;
    y[61] = TFPI_Xa_TF_VIIa_l;
    y[62] = APC_PS_Va_l;
    y[63] = IXa_AT_f;
    y[64] = Xa_AT_f;
    y[65] = VII_Xa_l;
    y[66] = V_mIIa_l;
    y[67] = VIII_mIIa_l;
    y[68] = TM_l;
    y[69] = IIa_TM_l;
    y[70] = IIa_TM_PC_l;
    y[71] = mIIa_AT_l;
    y[72] = XIa_IX_l;
    y[73] = LIPID;
}