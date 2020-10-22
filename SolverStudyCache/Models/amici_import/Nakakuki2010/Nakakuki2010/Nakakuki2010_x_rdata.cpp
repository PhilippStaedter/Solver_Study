#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Nakakuki2010(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = EGF;
    x_rdata[1] = HRG;
    x_rdata[2] = A1;
    x_rdata[3] = A1_2;
    x_rdata[4] = A2;
    x_rdata[5] = A2_2;
    x_rdata[6] = A3;
    x_rdata[7] = A3_2;
    x_rdata[8] = DUSPmRNA;
    x_rdata[9] = ERK_c;
    x_rdata[10] = pERK_c;
    x_rdata[11] = ppERK_c;
    x_rdata[12] = F;
    x_rdata[13] = c_FOS_c;
    x_rdata[14] = pc_FOS_c;
    x_rdata[15] = c_FOSmRNA;
    x_rdata[16] = FmRNA;
    x_rdata[17] = Kin;
    x_rdata[18] = Kin_2;
    x_rdata[19] = pMEK;
    x_rdata[20] = MEK;
    x_rdata[21] = DUSP_c;
    x_rdata[22] = pDUSP_c;
    x_rdata[23] = RSK_c;
    x_rdata[24] = pRSK_c;
    x_rdata[25] = RsD;
    x_rdata[26] = RsT;
    x_rdata[27] = CREB_n;
    x_rdata[28] = pCREB_n;
    x_rdata[29] = ERK_n;
    x_rdata[30] = pERK_n;
    x_rdata[31] = ppERK_n;
    x_rdata[32] = Elk1_n;
    x_rdata[33] = pElk1_n;
    x_rdata[34] = FOSn;
    x_rdata[35] = FOSn_2;
    x_rdata[36] = Fn;
    x_rdata[37] = DUSP_n;
    x_rdata[38] = pDUSP_n;
    x_rdata[39] = pDUSP_n_ERK_n;
    x_rdata[40] = pDUSP_n_pERK_n;
    x_rdata[41] = pDUSP_n_ppERK_n;
    x_rdata[42] = DUSP_n_ERK_n;
    x_rdata[43] = DUSP_n_pERK_n;
    x_rdata[44] = DUSP_n_ppERK_n;
    x_rdata[45] = PreDUSPmRNA;
    x_rdata[46] = PreFOSmRNA;
    x_rdata[47] = PreFmRNA;
    x_rdata[48] = pRSK_n;
}