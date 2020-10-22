#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Nakakuki2010(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = EGF;
    x_solver[1] = HRG;
    x_solver[2] = A1;
    x_solver[3] = A1_2;
    x_solver[4] = A2;
    x_solver[5] = A2_2;
    x_solver[6] = A3;
    x_solver[7] = A3_2;
    x_solver[8] = DUSPmRNA;
    x_solver[9] = ERK_c;
    x_solver[10] = pERK_c;
    x_solver[11] = ppERK_c;
    x_solver[12] = F;
    x_solver[13] = c_FOS_c;
    x_solver[14] = pc_FOS_c;
    x_solver[15] = c_FOSmRNA;
    x_solver[16] = FmRNA;
    x_solver[17] = Kin;
    x_solver[18] = Kin_2;
    x_solver[19] = pMEK;
    x_solver[20] = MEK;
    x_solver[21] = DUSP_c;
    x_solver[22] = pDUSP_c;
    x_solver[23] = RSK_c;
    x_solver[24] = pRSK_c;
    x_solver[25] = RsD;
    x_solver[26] = RsT;
    x_solver[27] = CREB_n;
    x_solver[28] = pCREB_n;
    x_solver[29] = ERK_n;
    x_solver[30] = pERK_n;
    x_solver[31] = ppERK_n;
    x_solver[32] = Elk1_n;
    x_solver[33] = pElk1_n;
    x_solver[34] = FOSn;
    x_solver[35] = FOSn_2;
    x_solver[36] = Fn;
    x_solver[37] = DUSP_n;
    x_solver[38] = pDUSP_n;
    x_solver[39] = pDUSP_n_ERK_n;
    x_solver[40] = pDUSP_n_pERK_n;
    x_solver[41] = pDUSP_n_ppERK_n;
    x_solver[42] = DUSP_n_ERK_n;
    x_solver[43] = DUSP_n_pERK_n;
    x_solver[44] = DUSP_n_ppERK_n;
    x_solver[45] = PreDUSPmRNA;
    x_solver[46] = PreFOSmRNA;
    x_solver[47] = PreFmRNA;
    x_solver[48] = pRSK_n;
}