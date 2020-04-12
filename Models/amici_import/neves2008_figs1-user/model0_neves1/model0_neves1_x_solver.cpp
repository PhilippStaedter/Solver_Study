#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "x_rdata.h"

void x_solver_model0_neves1(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = AC_PKA_cyto_mem;
    x_solver[1] = AC_active_cyto_mem;
    x_solver[2] = AC_cyto_mem;
    x_solver[3] = AMP_cyto;
    x_solver[4] = ATP_cyto;
    x_solver[5] = BAR_G_cyto_mem;
    x_solver[6] = BAR_cyto_mem;
    x_solver[7] = B_Raf_active_cyto;
    x_solver[8] = B_Raf_cyto;
    x_solver[9] = GRK_bg_cyto;
    x_solver[10] = GRK_cyto;
    x_solver[11] = G_GDP_cyto;
    x_solver[12] = G_a_s_cyto;
    x_solver[13] = G_protein_cyto;
    x_solver[14] = MAPK_active_cyto;
    x_solver[15] = MAPK_cyto;
    x_solver[16] = MEK_active_cyto;
    x_solver[17] = MEK_cyto;
    x_solver[18] = PDE4_P_cyto;
    x_solver[19] = PDE4_cyto;
    x_solver[20] = PDE_high_km_cyto;
    x_solver[21] = PKA_cyto;
    x_solver[22] = PP2A_cyto;
    x_solver[23] = PP_PDE_cyto;
    x_solver[24] = PTP_PKA_cyto;
    x_solver[25] = PTP_PP_cyto;
    x_solver[26] = PTP_cyto;
    x_solver[27] = R2C2_cyto;
    x_solver[28] = bg_cyto;
    x_solver[29] = c2_R2C2_cyto;
    x_solver[30] = c3_R2C2_cyto;
    x_solver[31] = cAMP_cyto;
    x_solver[32] = c_R2C2_cyto;
    x_solver[33] = iso_BAR_G_cyto_mem;
    x_solver[34] = iso_BAR_cyto_mem;
    x_solver[35] = iso_BAR_p_cyto_mem;
    x_solver[36] = iso_extra;
}