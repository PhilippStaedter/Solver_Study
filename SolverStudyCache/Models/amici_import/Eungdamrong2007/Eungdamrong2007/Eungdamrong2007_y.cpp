#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_Eungdamrong2007(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = RasGTP_Golgi_GM;
    y[1] = EGF_EC;
    y[2] = CAPRI_cyt;
    y[3] = serca;
    y[4] = PIP_PM;
    y[5] = PI_PM;
    y[6] = Shc_PM;
    y[7] = CaCAPRI_PM_PM;
    y[8] = RactCa;
    y[9] = Shc_star_PM;
    y[10] = EGFR_PM;
    y[11] = PLC_act_PM;
    y[12] = RasGTP_pal_cyt;
    y[13] = PLC_PM;
    y[14] = PIP2_PM;
    y[15] = Activated_EGFR_PM;
    y[16] = Ca;
    y[17] = Ract;
    y[18] = Rinh;
    y[19] = RinhCa;
    y[20] = IP3;
    y[21] = RasGDP_Golgi_GM;
    y[22] = Ca_RasGRP_GM_GM;
    y[23] = DAG_GM_GM;
    y[24] = RasGRP_DAG_GM;
    y[25] = CaCAPRI_cyt;
    y[26] = DAG_PM;
    y[27] = RasGTP_depal_cyt;
    y[28] = RasGDP_depal_cyt;
    y[29] = RasGDP_pal_cyt;
    y[30] = Ca_PLCe_cyt;
    y[31] = Ras_CaPLCe_GM;
    y[32] = PIP2_GM_GM;
    y[33] = ER_erMembrane;
    y[34] = Ca_ER;
    y[35] = Sos_cyt;
    y[36] = Grb2_cyt;
    y[37] = PLCe_cyt;
    y[38] = buffer_cyt;
    y[39] = ca_buffer_cyt;
    y[40] = SosGrb2_cyt;
    y[41] = SGS_PM;
    y[42] = RasGTP_PM;
    y[43] = RasGDP_PM;
    y[44] = RasGRP_cyt;
    y[45] = CaRasGRP1_cyt;
}