#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_hockin1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = ATIII;
    x_solver[1] = II;
    x_solver[2] = IIa;
    x_solver[3] = IIa_ATIII;
    x_solver[4] = IX;
    x_solver[5] = IXa;
    x_solver[6] = IXa_ATIII;
    x_solver[7] = IXa_VIIIa;
    x_solver[8] = IXa_VIIIa_X;
    x_solver[9] = TF;
    x_solver[10] = TFPI;
    x_solver[11] = TF_VII;
    x_solver[12] = TF_VIIa;
    x_solver[13] = TF_VIIa_ATIII;
    x_solver[14] = TF_VIIa_IX;
    x_solver[15] = TF_VIIa_X;
    x_solver[16] = TF_VIIa_Xa;
    x_solver[17] = TF_VIIa_Xa_TFPI;
    x_solver[18] = V;
    x_solver[19] = VII;
    x_solver[20] = VIII;
    x_solver[21] = VIIIa;
    x_solver[22] = VIIIa1_L;
    x_solver[23] = VIIIa2;
    x_solver[24] = VIIa;
    x_solver[25] = Va;
    x_solver[26] = X;
    x_solver[27] = Xa;
    x_solver[28] = Xa_ATIII;
    x_solver[29] = Xa_TFPI;
    x_solver[30] = Xa_Va;
    x_solver[31] = Xa_Va_II;
    x_solver[32] = mIIa;
    x_solver[33] = mIIa_ATIII;
}