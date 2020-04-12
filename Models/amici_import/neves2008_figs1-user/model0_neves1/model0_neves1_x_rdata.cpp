#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "x.h"

void x_rdata_model0_neves1(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = AC_PKA_cyto_mem;
    x_rdata[1] = AC_active_cyto_mem;
    x_rdata[2] = AC_cyto_mem;
    x_rdata[3] = AMP_cyto;
    x_rdata[4] = ATP_cyto;
    x_rdata[5] = BAR_G_cyto_mem;
    x_rdata[6] = BAR_cyto_mem;
    x_rdata[7] = B_Raf_active_cyto;
    x_rdata[8] = B_Raf_cyto;
    x_rdata[9] = GRK_bg_cyto;
    x_rdata[10] = GRK_cyto;
    x_rdata[11] = G_GDP_cyto;
    x_rdata[12] = G_a_s_cyto;
    x_rdata[13] = G_protein_cyto;
    x_rdata[14] = MAPK_active_cyto;
    x_rdata[15] = MAPK_cyto;
    x_rdata[16] = MEK_active_cyto;
    x_rdata[17] = MEK_cyto;
    x_rdata[18] = PDE4_P_cyto;
    x_rdata[19] = PDE4_cyto;
    x_rdata[20] = PDE_high_km_cyto;
    x_rdata[21] = PKA_cyto;
    x_rdata[22] = PP2A_cyto;
    x_rdata[23] = PP_PDE_cyto;
    x_rdata[24] = PTP_PKA_cyto;
    x_rdata[25] = PTP_PP_cyto;
    x_rdata[26] = PTP_cyto;
    x_rdata[27] = R2C2_cyto;
    x_rdata[28] = bg_cyto;
    x_rdata[29] = c2_R2C2_cyto;
    x_rdata[30] = c3_R2C2_cyto;
    x_rdata[31] = cAMP_cyto;
    x_rdata[32] = c_R2C2_cyto;
    x_rdata[33] = iso_BAR_G_cyto_mem;
    x_rdata[34] = iso_BAR_cyto_mem;
    x_rdata[35] = iso_BAR_p_cyto_mem;
    x_rdata[36] = iso_extra;
}