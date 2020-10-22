#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_Eungdamrong2007(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = RasGTP_Golgi_GM;
    x_rdata[1] = EGF_EC;
    x_rdata[2] = CAPRI_cyt;
    x_rdata[3] = serca;
    x_rdata[4] = PIP_PM;
    x_rdata[5] = PI_PM;
    x_rdata[6] = Shc_PM;
    x_rdata[7] = CaCAPRI_PM_PM;
    x_rdata[8] = RactCa;
    x_rdata[9] = Shc_star_PM;
    x_rdata[10] = EGFR_PM;
    x_rdata[11] = PLC_act_PM;
    x_rdata[12] = RasGTP_pal_cyt;
    x_rdata[13] = PLC_PM;
    x_rdata[14] = PIP2_PM;
    x_rdata[15] = Activated_EGFR_PM;
    x_rdata[16] = Ca;
    x_rdata[17] = Ract;
    x_rdata[18] = Rinh;
    x_rdata[19] = RinhCa;
    x_rdata[20] = IP3;
    x_rdata[21] = RasGDP_Golgi_GM;
    x_rdata[22] = Ca_RasGRP_GM_GM;
    x_rdata[23] = DAG_GM_GM;
    x_rdata[24] = RasGRP_DAG_GM;
    x_rdata[25] = CaCAPRI_cyt;
    x_rdata[26] = DAG_PM;
    x_rdata[27] = RasGTP_depal_cyt;
    x_rdata[28] = RasGDP_depal_cyt;
    x_rdata[29] = RasGDP_pal_cyt;
    x_rdata[30] = Ca_PLCe_cyt;
    x_rdata[31] = Ras_CaPLCe_GM;
    x_rdata[32] = PIP2_GM_GM;
    x_rdata[33] = ER_erMembrane;
    x_rdata[34] = Ca_ER;
    x_rdata[35] = Sos_cyt;
    x_rdata[36] = Grb2_cyt;
    x_rdata[37] = PLCe_cyt;
    x_rdata[38] = buffer_cyt;
    x_rdata[39] = ca_buffer_cyt;
    x_rdata[40] = SosGrb2_cyt;
    x_rdata[41] = SGS_PM;
    x_rdata[42] = RasGTP_PM;
    x_rdata[43] = RasGDP_PM;
    x_rdata[44] = RasGRP_cyt;
    x_rdata[45] = CaRasGRP1_cyt;
}