#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model4_hockin1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = ATIII;
    y[1] = II;
    y[2] = IIa;
    y[3] = IIa_ATIII;
    y[4] = IX;
    y[5] = IXa;
    y[6] = IXa_ATIII;
    y[7] = IXa_VIIIa;
    y[8] = IXa_VIIIa_X;
    y[9] = TF;
    y[10] = TFPI;
    y[11] = TF_VII;
    y[12] = TF_VIIa;
    y[13] = TF_VIIa_ATIII;
    y[14] = TF_VIIa_IX;
    y[15] = TF_VIIa_X;
    y[16] = TF_VIIa_Xa;
    y[17] = TF_VIIa_Xa_TFPI;
    y[18] = V;
    y[19] = VII;
    y[20] = VIII;
    y[21] = VIIIa;
    y[22] = VIIIa1_L;
    y[23] = VIIIa2;
    y[24] = VIIa;
    y[25] = Va;
    y[26] = X;
    y[27] = Xa;
    y[28] = Xa_ATIII;
    y[29] = Xa_TFPI;
    y[30] = Xa_Va;
    y[31] = Xa_Va_II;
    y[32] = mIIa;
    y[33] = mIIa_ATIII;
}