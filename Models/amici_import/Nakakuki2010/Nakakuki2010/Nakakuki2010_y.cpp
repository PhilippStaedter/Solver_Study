#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Nakakuki2010(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = EGF;
    y[1] = HRG;
    y[2] = A1;
    y[3] = A1_2;
    y[4] = A2;
    y[5] = A2_2;
    y[6] = A3;
    y[7] = A3_2;
    y[8] = DUSPmRNA;
    y[9] = ERK_c;
    y[10] = pERK_c;
    y[11] = ppERK_c;
    y[12] = F;
    y[13] = c_FOS_c;
    y[14] = pc_FOS_c;
    y[15] = c_FOSmRNA;
    y[16] = FmRNA;
    y[17] = Kin;
    y[18] = Kin_2;
    y[19] = pMEK;
    y[20] = MEK;
    y[21] = DUSP_c;
    y[22] = pDUSP_c;
    y[23] = RSK_c;
    y[24] = pRSK_c;
    y[25] = RsD;
    y[26] = RsT;
    y[27] = CREB_n;
    y[28] = pCREB_n;
    y[29] = ERK_n;
    y[30] = pERK_n;
    y[31] = ppERK_n;
    y[32] = Elk1_n;
    y[33] = pElk1_n;
    y[34] = FOSn;
    y[35] = FOSn_2;
    y[36] = Fn;
    y[37] = DUSP_n;
    y[38] = pDUSP_n;
    y[39] = pDUSP_n_ERK_n;
    y[40] = pDUSP_n_pERK_n;
    y[41] = pDUSP_n_ppERK_n;
    y[42] = DUSP_n_ERK_n;
    y[43] = DUSP_n_pERK_n;
    y[44] = DUSP_n_ppERK_n;
    y[45] = PreDUSPmRNA;
    y[46] = PreFOSmRNA;
    y[47] = PreFmRNA;
    y[48] = pRSK_n;
}