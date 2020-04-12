#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"
#include <cmath>


#include "w.h"
#include "tcl.h"
#include "p.h"
#include "k.h"
#include "x.h"

void dwdx_model0_neves1(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dwdx[0] = -0.20000000000000001*AC_activation_Kr_AC_activation;
    dwdx[1] = 0.002823920265780727*ATP_cyto/(AC_active_Km_AC_active + 0.0016611295681063099*ATP_cyto);
    dwdx[2] = 0.00033222591362126199*AC_activation_Kf_AC_activation*G_a_s_cyto;
    dwdx[3] = 6.6445182724252401e-5*ATP_cyto/(AC_basal_Km_AC_basal + 0.0016611295681063099*ATP_cyto);
    dwdx[4] = -4.6908974514629954e-6*AC_active_cyto_mem*ATP_cyto/pow(AC_active_Km_AC_active + 0.0016611295681063099*ATP_cyto, 2) + 0.002823920265780727*AC_active_cyto_mem/(AC_active_Km_AC_active + 0.0016611295681063099*ATP_cyto);
    dwdx[5] = -1.1037405768148224e-7*AC_cyto_mem*ATP_cyto/pow(AC_basal_Km_AC_basal + 0.0016611295681063099*ATP_cyto, 2) + 6.6445182724252401e-5*AC_cyto_mem/(AC_basal_Km_AC_basal + 0.0016611295681063099*ATP_cyto);
    dwdx[6] = -0.20000000000000001*G_binds_BAR_Kr_G_binds_BAR;
    dwdx[7] = 0.00033222591362126199*iso_binds_BAR_g_Kf*iso_extra;
    dwdx[8] = 0.00033222591362126199*G_binds_BAR_Kf_G_binds_BAR*G_protein_cyto;
    dwdx[9] = 0.00033222591362126199*iso_binds_BAR_Kf*iso_extra;
    dwdx[10] = -4.5836402691645375e-9*B_Raf_active_cyto*PP2A_cyto*kcat_PPase_Raf/(KMOLE*pow(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km, 2)) + 2.7593514420370556e-6*PP2A_cyto*kcat_PPase_Raf/(KMOLE*(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km));
    dwdx[11] = 2.8973190141389085e-7*MEK_cyto/(KMOLE*(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km));
    dwdx[12] = -4.5836402691645375e-9*B_Raf_cyto*PKA_cyto*kcat_PKA_activates_Raf/(KMOLE*pow(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km, 2)) + 2.7593514420370556e-6*PKA_cyto*kcat_PKA_activates_Raf/(KMOLE*(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km));
    dwdx[13] = 0.00044518272425249112*iso_BAR_cyto_mem/(GRK_bg_Km_GRK_bg + iso_BAR_cyto_mem);
    dwdx[14] = -0.0016611295681063099*bg_binds_GRK_Kr_bg_binds_GRK/KMOLE;
    dwdx[15] = 3.4551495016611247e-5*iso_BAR_cyto_mem/(GRK_Km_grk + iso_BAR_cyto_mem);
    dwdx[16] = 2.7593514420370556e-6*bg_binds_GRK_Kf_bg_binds_GRK*bg_cyto/KMOLE;
    dwdx[17] = -0.0016611295681063099*GTPase_Kr_GTPase/KMOLE;
    dwdx[18] = 2.7593514420370556e-6*bg_cyto*trimer_Kf_trimer/KMOLE;
    dwdx[19] = 0.00033222591362126199*AC_activation_Kf_AC_activation*AC_cyto_mem;
    dwdx[20] = 0.0016611295681063099*GTPase_Kf_GTPase/KMOLE;
    dwdx[21] = -5.5187028840741116e-7*activate_Gs_Kr_activate_Gs*bg_cyto*iso_BAR_cyto_mem;
    dwdx[22] = 0.00033222591362126199*BAR_cyto_mem*G_binds_BAR_Kf_G_binds_BAR;
    dwdx[23] = 0.00033222591362126199*G_binds_iso_BAR_Kf_G_binds_iso_BAR*iso_BAR_cyto_mem;
    dwdx[24] = -0.0016611295681063099*trimer_Kr_trimer/KMOLE;
    dwdx[25] = -4.5836402691645375e-9*MAPK_active_cyto*PP2A_cyto*kcat_PPase_MAPK/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km, 2)) + 2.7593514420370556e-6*PP2A_cyto*kcat_PPase_MAPK/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km));
    dwdx[26] = -4.5836402691645375e-9*MAPK_active_cyto*PTP_cyto*kcat_PTP/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PTP_Km, 2)) + 2.7593514420370556e-6*PTP_cyto*kcat_PTP/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_Km));
    dwdx[27] = -4.5836402691645375e-9*MAPK_active_cyto*PTP_PKA_cyto*kcat_PTP_PKA/(KMOLE*pow(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km, 2)) + 2.7593514420370556e-6*PTP_PKA_cyto*kcat_PTP_PKA/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km));
    dwdx[28] = -6.8754604037468067e-10*MAPK_cyto*MEK_active_cyto/(KMOLE*pow(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km, 2)) + 4.1390271630555837e-7*MEK_active_cyto/(KMOLE*(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km));
    dwdx[29] = 4.1390271630555837e-7*MAPK_cyto/(KMOLE*(0.0016611295681063099*MAPK_cyto + MEK_activates_MAPK_Km));
    dwdx[30] = -4.5836402691645375e-9*MEK_active_cyto*PP2A_cyto*kcat_PPase_mek/(KMOLE*pow(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km, 2)) + 2.7593514420370556e-6*PP2A_cyto*kcat_PPase_mek/(KMOLE*(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km));
    dwdx[31] = -4.8128222826227647e-10*B_Raf_active_cyto*MEK_cyto/(KMOLE*pow(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km, 2)) + 2.8973190141389085e-7*B_Raf_active_cyto/(KMOLE*(0.0016611295681063099*MEK_cyto + Raf_activates_MEK_Km));
    dwdx[32] = 5.5187028840741121e-5*cAMP_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p));
    dwdx[33] = -2.2918201345822693e-8*PDE4_P_cyto*PP_PDE_cyto/(KMOLE*pow(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4, 2)) + 1.379675721018528e-5*PP_PDE_cyto/(KMOLE*(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4));
    dwdx[34] = 2.2074811536296445e-5*cAMP_cyto/(KMOLE*(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto));
    dwdx[35] = -4.5836402691645385e-8*PDE4_cyto*PKA_cyto/(KMOLE*pow(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km, 2)) + 2.759351442037056e-5*PKA_cyto/(KMOLE*(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km));
    dwdx[36] = 2.2074811536296445e-5*cAMP_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km));
    dwdx[37] = -0.0016611295681063099*A2_Kr/KMOLE;
    dwdx[38] = 2.759351442037056e-5*PDE4_cyto/(KMOLE*(0.0016611295681063099*PDE4_cyto + PKA_P_PDE_Km));
    dwdx[39] = 2.7593514420370556e-6*PTP_cyto*kcat_PKA_P_PTP/(KMOLE*(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto));
    dwdx[40] = 2.7593514420370556e-6*B_Raf_cyto*kcat_PKA_activates_Raf/(KMOLE*(0.0016611295681063099*B_Raf_cyto + PKA_activates_Raf_Km));
    dwdx[41] = 2.7593514420370556e-6*MAPK_active_cyto*kcat_PPase_MAPK/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PPase_MAPK_Km));
    dwdx[42] = 2.7593514420370556e-6*B_Raf_active_cyto*kcat_PPase_Raf/(KMOLE*(0.0016611295681063099*B_Raf_active_cyto + PPase_Raf_Km));
    dwdx[43] = 2.7593514420370556e-6*MEK_active_cyto*kcat_PPase_mek/(KMOLE*(0.0016611295681063099*MEK_active_cyto + PPase_mek_Km));
    dwdx[44] = 1.379675721018528e-5*PDE4_P_cyto/(KMOLE*(0.0016611295681063099*PDE4_P_cyto + pp2a_4_Km_pp2a_4));
    dwdx[45] = 2.7593514420370556e-6*MAPK_active_cyto*kcat_PTP_PKA/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_PKA_Km));
    dwdx[46] = -4.5836402691645375e-9*PTP_PKA_cyto*PTP_PP_cyto*kcat_pp_ptp_pp_ptp/(KMOLE*pow(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km, 2)) + 2.7593514420370556e-6*PTP_PP_cyto*kcat_pp_ptp_pp_ptp/(KMOLE*(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km));
    dwdx[47] = 2.7593514420370556e-6*PTP_PKA_cyto*kcat_pp_ptp_pp_ptp/(KMOLE*(0.0016611295681063099*PTP_PKA_cyto + pp_ptp_Km));
    dwdx[48] = -4.5836402691645375e-9*PKA_cyto*PTP_cyto*kcat_PKA_P_PTP/(KMOLE*pow(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto, 2)) + 2.7593514420370556e-6*PKA_cyto*kcat_PKA_P_PTP/(KMOLE*(PKA_P_PTP_Km + 0.0016611295681063099*PTP_cyto));
    dwdx[49] = 2.7593514420370556e-6*MAPK_active_cyto*kcat_PTP/(KMOLE*(0.0016611295681063099*MAPK_active_cyto + PTP_Km));
    dwdx[50] = 2.7593514420370556e-6*B1_Kf*cAMP_cyto/KMOLE;
    dwdx[51] = -5.5187028840741116e-7*G_a_s_cyto*activate_Gs_Kr_activate_Gs*iso_BAR_cyto_mem;
    dwdx[52] = 2.7593514420370556e-6*GRK_cyto*bg_binds_GRK_Kf_bg_binds_GRK/KMOLE;
    dwdx[53] = 2.7593514420370556e-6*G_GDP_cyto*trimer_Kf_trimer/KMOLE;
    dwdx[54] = 2.7593514420370556e-6*A1_Kf*cAMP_cyto/KMOLE;
    dwdx[55] = -0.0016611295681063099*B2_Kr/KMOLE;
    dwdx[56] = -0.0016611295681063099*A1_Kr/KMOLE;
    dwdx[57] = 2.7593514420370556e-6*A2_Kf*cAMP_cyto/KMOLE;
    dwdx[58] = 2.7593514420370556e-6*A1_Kf*c2_R2C2_cyto/KMOLE;
    dwdx[59] = 2.7593514420370556e-6*A2_Kf*c3_R2C2_cyto/KMOLE;
    dwdx[60] = 2.7593514420370556e-6*B1_Kf*R2C2_cyto/KMOLE;
    dwdx[61] = 2.7593514420370556e-6*B2_Kf*c_R2C2_cyto/KMOLE;
    dwdx[62] = -3.66691221533163e-8*PDE4_cyto*cAMP_cyto/(KMOLE*pow(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto, 2)) + 2.2074811536296445e-5*PDE4_cyto/(KMOLE*(PDE4_Km_PDE4 + 0.0016611295681063099*cAMP_cyto));
    dwdx[63] = -3.66691221533163e-8*PDE_high_km_cyto*cAMP_cyto/(KMOLE*pow(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km, 2)) + 2.2074811536296445e-5*PDE_high_km_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + highKM_PDE_Km));
    dwdx[64] = -9.167280538329077e-8*PDE4_P_cyto*cAMP_cyto/(KMOLE*pow(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p, 2)) + 5.5187028840741121e-5*PDE4_P_cyto/(KMOLE*(0.0016611295681063099*cAMP_cyto + pde4_p_Km_pde4_p));
    dwdx[65] = -0.0016611295681063099*B1_Kr/KMOLE;
    dwdx[66] = 2.7593514420370556e-6*B2_Kf*cAMP_cyto/KMOLE;
    dwdx[67] = -0.20000000000000001*G_binds_iso_BAR_Kr_G_binds_iso_BAR;
    dwdx[68] = 0.20000000000000001*activate_Gs_Kf_activate_Gs;
    dwdx[69] = -0.20000000000000001*iso_binds_BAR_g_Kr;
    dwdx[70] = -3.4551495016611247e-5*GRK_cyto*iso_BAR_cyto_mem/pow(GRK_Km_grk + iso_BAR_cyto_mem, 2) + 3.4551495016611247e-5*GRK_cyto/(GRK_Km_grk + iso_BAR_cyto_mem);
    dwdx[71] = -0.00044518272425249112*GRK_bg_cyto*iso_BAR_cyto_mem/pow(GRK_bg_Km_GRK_bg + iso_BAR_cyto_mem, 2) + 0.00044518272425249112*GRK_bg_cyto/(GRK_bg_Km_GRK_bg + iso_BAR_cyto_mem);
    dwdx[72] = 0.00033222591362126199*G_binds_iso_BAR_Kf_G_binds_iso_BAR*G_protein_cyto;
    dwdx[73] = -5.5187028840741116e-7*G_a_s_cyto*activate_Gs_Kr_activate_Gs*bg_cyto;
    dwdx[74] = -0.20000000000000001*iso_binds_BAR_Kr;
    dwdx[75] = 0.00033222591362126199*BAR_cyto_mem*iso_binds_BAR_Kf;
    dwdx[76] = 0.00033222591362126199*BAR_G_cyto_mem*iso_binds_BAR_g_Kf;
}