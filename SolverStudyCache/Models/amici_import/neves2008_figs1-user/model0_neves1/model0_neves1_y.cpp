#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "p.h"
#include "k.h"
#include "x.h"

void y_model0_neves1(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = AC_PKA_cyto_mem;
    y[1] = AC_active_cyto_mem;
    y[2] = AC_cyto_mem;
    y[3] = AMP_cyto;
    y[4] = ATP_cyto;
    y[5] = BAR_G_cyto_mem;
    y[6] = BAR_cyto_mem;
    y[7] = B_Raf_active_cyto;
    y[8] = B_Raf_cyto;
    y[9] = GRK_bg_cyto;
    y[10] = GRK_cyto;
    y[11] = G_GDP_cyto;
    y[12] = G_a_s_cyto;
    y[13] = G_protein_cyto;
    y[14] = MAPK_active_cyto;
    y[15] = MAPK_cyto;
    y[16] = MEK_active_cyto;
    y[17] = MEK_cyto;
    y[18] = PDE4_P_cyto;
    y[19] = PDE4_cyto;
    y[20] = PDE_high_km_cyto;
    y[21] = PKA_cyto;
    y[22] = PP2A_cyto;
    y[23] = PP_PDE_cyto;
    y[24] = PTP_PKA_cyto;
    y[25] = PTP_PP_cyto;
    y[26] = PTP_cyto;
    y[27] = R2C2_cyto;
    y[28] = bg_cyto;
    y[29] = c2_R2C2_cyto;
    y[30] = c3_R2C2_cyto;
    y[31] = cAMP_cyto;
    y[32] = c_R2C2_cyto;
    y[33] = iso_BAR_G_cyto_mem;
    y[34] = iso_BAR_cyto_mem;
    y[35] = iso_BAR_p_cyto_mem;
    y[36] = iso_extra;
}