#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void w_model0_neves1(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    w[0] = 1.0*(2.7593514420370556e-6*A1_Kf*c2_R2C2_cyto*cAMP_cyto - 0.0016611295681063099*A1_Kr*c3_R2C2_cyto)/KMOLE;
    w[1] = 1.0*(2.7593514420370556e-6*A2_Kf*c3_R2C2_cyto*cAMP_cyto - 0.0016611295681063099*A2_Kr*PKA_cyto)/KMOLE;
    w[2] = 0.00033222591362126199*AC_activation_Kf_AC_activation*AC_cyto_mem*G_a_s_cyto - 0.20000000000000001*AC_activation_Kr_AC_activation*AC_active_cyto_mem;
    w[3] = 0.002823920265780727*AC_active_cyto_mem*ATP_cyto/(AC_active_Km_AC_active + 0.0016611295681063099*ATP_cyto);
    w[4] = 6.6445182724252401e-5*AC_cyto_mem*ATP_cyto/(AC_basal_Km_AC_basal + 0.0016611295681063099*ATP_cyto);
    w[5] = 1.0*(2.7593514420370556e-6*B1_Kf*R2C2_cyto*cAMP_cyto - 0.0016611295681063099*B1_Kr*c_R2C2_cyto)/KMOLE;
    w[6] = 1.0*(2.7593514420370556e-6*B2_Kf*cAMP_cyto*c_R2C2_cyto - 0.0016611295681063099*B2_Kr*c2_R2C2_cyto)/KMOLE;
    w[7] = 3.4551495016611247e-5*GRK_cyto*iso_BAR_cyto_mem/(GRK_Km_grk + iso_BAR_cyto_mem);
    w[8] = 0.00044518272425249112*GRK_bg_cyto*iso_BAR_cyto_mem/(GRK_bg_Km_GRK_bg + iso_BAR_cyto_mem);
    w[9] = 1.0*(0.0016611295681063099*GTPase_Kf_GTPase*G_a_s_cyto - 0.0016611295681063099*GTPase_Kr_GTPase*G_GDP_cyto)/KMOLE;
    w[10] = -0.20000000000000001*BAR_G_cyto_mem*G_binds_BAR_Kr_G_binds_BAR + 0.00033222591362126199*BAR_cyto_mem*G_binds_BAR_Kf_G_binds_BAR*G_protein_cyto;
    w[11] = 0.00033222591362126199*G_binds_iso_BAR_Kf_G_binds_iso_BAR*G_protein_cyto*iso_BAR_cyto_mem - 0.20000000000000001*G_binds_iso_BAR_Kr_G_binds_iso_BAR*iso_BAR_G_cyto_mem;
    w[12] = 4.1390271630555837e-7*MAPK_cyto*MEK_active_cyto/(KMOLE*(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km));
    w[13] = 2.2074811536296445e-5*PDE4_cyto*cAMP_cyto/(KMOLE*(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto));
    w[14] = 2.759351442037056e-5*PDE4_cyto*PKA_cyto/(KMOLE*(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km));
    w[15] = 2.7593514420370556e-6*PKA_cyto*PTP_cyto*kcat_PKA_P_PTP/(KMOLE*(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto));
    w[16] = 2.7593514420370556e-6*B_Raf_cyto*PKA_cyto*kcat_PKA_activates_Raf/(KMOLE*(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km));
    w[17] = 2.7593514420370556e-6*MAPK_active_cyto*PP2A_cyto*kcat_PPase_MAPK/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km));
    w[18] = 2.7593514420370556e-6*B_Raf_active_cyto*PP2A_cyto*kcat_PPase_Raf/(KMOLE*(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km));
    w[19] = 2.7593514420370556e-6*MEK_active_cyto*PP2A_cyto*kcat_PPase_mek/(KMOLE*(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km));
    w[20] = 2.7593514420370556e-6*MAPK_active_cyto*PTP_cyto*kcat_PTP/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_Km));
    w[21] = 2.7593514420370556e-6*MAPK_active_cyto*PTP_PKA_cyto*kcat_PTP_PKA/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km));
    w[22] = 2.8973190141389085e-7*B_Raf_active_cyto*MEK_cyto/(KMOLE*(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km));
    w[23] = -5.5187028840741116e-7*G_a_s_cyto*activate_Gs_Kr_activate_Gs*bg_cyto*iso_BAR_cyto_mem + 0.20000000000000001*activate_Gs_Kf_activate_Gs*iso_BAR_G_cyto_mem;
    w[24] = 1.0*(-0.0016611295681063099*GRK_bg_cyto*bg_binds_GRK_Kr_bg_binds_GRK + 2.7593514420370556e-6*GRK_cyto*bg_binds_GRK_Kf_bg_binds_GRK*bg_cyto)/KMOLE;
    w[25] = 2.2074811536296445e-5*PDE_high_km_cyto*cAMP_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km));
    w[26] = 0.00033222591362126199*BAR_cyto_mem*iso_binds_BAR_Kf*iso_extra - 0.20000000000000001*iso_BAR_cyto_mem*iso_binds_BAR_Kr;
    w[27] = 0.00033222591362126199*BAR_G_cyto_mem*iso_binds_BAR_g_Kf*iso_extra - 0.20000000000000001*iso_BAR_G_cyto_mem*iso_binds_BAR_g_Kr;
    w[28] = 5.5187028840741121e-5*PDE4_P_cyto*cAMP_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p));
    w[29] = 1.379675721018528e-5*PDE4_P_cyto*PP_PDE_cyto/(KMOLE*(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4));
    w[30] = 2.7593514420370556e-6*PTP_PKA_cyto*PTP_PP_cyto*kcat_pp_ptp_pp_ptp/(KMOLE*(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km));
    w[31] = 1.0*(2.7593514420370556e-6*G_GDP_cyto*bg_cyto*trimer_Kf_trimer - 0.0016611295681063099*G_protein_cyto*trimer_Kr_trimer)/KMOLE;
}