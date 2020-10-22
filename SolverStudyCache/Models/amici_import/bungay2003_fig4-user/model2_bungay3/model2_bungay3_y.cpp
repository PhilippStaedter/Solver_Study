#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model2_bungay3(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = APC_PS_VIIIa_l;
    y[1] = APC_PS_Va_l;
    y[2] = APC_PS_l;
    y[3] = APC_f;
    y[4] = APC_l;
    y[5] = AT_f;
    y[6] = II_f;
    y[7] = II_l;
    y[8] = IIa_AT_f;
    y[9] = IIa_TM_PC_l;
    y[10] = IIa_TM_l;
    y[11] = IIa_f;
    y[12] = IX_f;
    y[13] = IX_l;
    y[14] = IXa_AT_f;
    y[15] = IXa_VIIIa_X_l;
    y[16] = IXa_VIIIa_l;
    y[17] = IXa_f;
    y[18] = IXa_l;
    y[19] = LIPID;
    y[20] = PC_f;
    y[21] = PC_l;
    y[22] = PS_f;
    y[23] = PS_l;
    y[24] = TFPI_Xa_TF_VIIa_l;
    y[25] = TFPI_Xa_l;
    y[26] = TFPI_f;
    y[27] = TF_VII_Xa_l;
    y[28] = TF_VII_l;
    y[29] = TF_VIIa_IX_l;
    y[30] = TF_VIIa_IXa_l;
    y[31] = TF_VIIa_X_l;
    y[32] = TF_VIIa_Xa_l;
    y[33] = TF_VIIa_l;
    y[34] = TF_l;
    y[35] = TM_l;
    y[36] = VIII_IIa_l;
    y[37] = VIII_Xa_l;
    y[38] = VIII_f;
    y[39] = VIII_l;
    y[40] = VIII_mIIa_l;
    y[41] = VIIIa_f;
    y[42] = VIIIa_l;
    y[43] = VIIIai_f;
    y[44] = VIIIai_l;
    y[45] = VII_Xa_l;
    y[46] = VII_f;
    y[47] = VII_l;
    y[48] = VIIa_f;
    y[49] = VIIa_l;
    y[50] = V_IIa_l;
    y[51] = V_Xa_l;
    y[52] = V_f;
    y[53] = V_l;
    y[54] = V_mIIa_l;
    y[55] = Va_f;
    y[56] = Va_l;
    y[57] = Vai_f;
    y[58] = Vai_l;
    y[59] = XI_IIa_l;
    y[60] = XI_f;
    y[61] = XIa_IX_l;
    y[62] = XIa_l;
    y[63] = X_f;
    y[64] = X_l;
    y[65] = Xa_AT_f;
    y[66] = Xa_Va_II_l;
    y[67] = Xa_Va_l;
    y[68] = Xa_Va_mIIa_l;
    y[69] = Xa_f;
    y[70] = Xa_l;
    y[71] = mIIa_AT_l;
    y[72] = mIIa_f;
    y[73] = mIIa_l;
}