#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_Eungdamrong2007(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = RasGTP_Golgi_GM;
    x_solver[1] = EGF_EC;
    x_solver[2] = CAPRI_cyt;
    x_solver[3] = serca;
    x_solver[4] = PIP_PM;
    x_solver[5] = PI_PM;
    x_solver[6] = Shc_PM;
    x_solver[7] = CaCAPRI_PM_PM;
    x_solver[8] = RactCa;
    x_solver[9] = Shc_star_PM;
    x_solver[10] = EGFR_PM;
    x_solver[11] = PLC_act_PM;
    x_solver[12] = RasGTP_pal_cyt;
    x_solver[13] = PLC_PM;
    x_solver[14] = PIP2_PM;
    x_solver[15] = Activated_EGFR_PM;
    x_solver[16] = Ca;
    x_solver[17] = Ract;
    x_solver[18] = Rinh;
    x_solver[19] = RinhCa;
    x_solver[20] = IP3;
    x_solver[21] = RasGDP_Golgi_GM;
    x_solver[22] = Ca_RasGRP_GM_GM;
    x_solver[23] = DAG_GM_GM;
    x_solver[24] = RasGRP_DAG_GM;
    x_solver[25] = CaCAPRI_cyt;
    x_solver[26] = DAG_PM;
    x_solver[27] = RasGTP_depal_cyt;
    x_solver[28] = RasGDP_depal_cyt;
    x_solver[29] = RasGDP_pal_cyt;
    x_solver[30] = Ca_PLCe_cyt;
    x_solver[31] = Ras_CaPLCe_GM;
    x_solver[32] = PIP2_GM_GM;
    x_solver[33] = ER_erMembrane;
    x_solver[34] = Ca_ER;
    x_solver[35] = Sos_cyt;
    x_solver[36] = Grb2_cyt;
    x_solver[37] = PLCe_cyt;
    x_solver[38] = buffer_cyt;
    x_solver[39] = ca_buffer_cyt;
    x_solver[40] = SosGrb2_cyt;
    x_solver[41] = SGS_PM;
    x_solver[42] = RasGTP_PM;
    x_solver[43] = RasGDP_PM;
    x_solver[44] = RasGRP_cyt;
    x_solver[45] = CaRasGRP1_cyt;
}